#include "palabos3D.h"
#include "palabos3D.hh"
#include <cstdlib>
#include <iostream>

using namespace plb;
using namespace std;

// Use double-precision arithmetics
typedef double T;
// Use a grid which additionally to the f's stores two variables for
//   the external force term.
#define DESCRIPTOR descriptors::ShanChenD3Q19Descriptor

T kinematicViscosity        = 0.;
T averageInletVelocity      = 0.;
plint referenceResolution   = 0.;

T uAveLB                    = 0.;
T fluidDensity              = 0.;
T volume                    = 0.;
T Re                        = 0.;
T G                         = 0.;

T simTime = 0;
plint maxIter = 0;
plint writeInterval = 0;
T epsilon = 0;
bool performOutput = false;

T currentTime = 0;

/// Initial condition: heavy fluid on top, light fluid on bottom.
/** This functional is going to be used as an argument to the function "applyIndexed",
 *  to setup the initial condition. For efficiency reasons, this approach should
 *  always be preferred over explicit space loops in end-user codes.
 */
template<typename T, template<typename U> class Descriptor>
class TwoLayerInitializer : public OneCellIndexedWithRandFunctional3D<T,Descriptor> {
public:
    TwoLayerInitializer(plint nx_, bool topLayer_)
        : nx(nx_),
          topLayer(topLayer_)
    { }
    TwoLayerInitializer<T,Descriptor>* clone() const {
        return new TwoLayerInitializer<T,Descriptor>(*this);
    }
    virtual void execute(plint iX, plint iY, plint iZ, T rand_val, Cell<T,Descriptor>& cell) const {
        T densityFluctuations = 1.e-2;
        T almostNoFluid       = 1.e-4;
        Array<T,3> zeroVelocity (0.,0.,0.);

        T rho = (T)1;
        // Add a random perturbation to the initial condition to instantiate the
        //   instability.
        if ( (topLayer && iX>=nx/2) || (!topLayer && iX <nx/2) ) {
            rho += rand_val * densityFluctuations;
        }
        else {
            rho = almostNoFluid;
        }

        iniCellAtEquilibrium(cell, rho, zeroVelocity);
    }
private:
    plint nx;
    bool topLayer;
};

// InstantiateDynamicsFunctional2D
template<typename T, template<typename U> class Descriptor>
class SetDensity : public BoxProcessingFunctional3D_L<T,Descriptor> {
public:
    SetDensity(plint dx_, plint dy_, plint dz_)
        : dx(dx_), dy(dy_), dz(dz_)
    { }
    SetDensity<T,Descriptor>* clone() const {
        return new SetDensity<T,Descriptor>(*this);
    }
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice) {
        T rho = (T)1;
        for (plint iX=domain.x0; iX<=domain.x1; ++iX){
            for (plint iY=domain.y0; iY<=domain.y1; ++iY){
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ){
                    rho = lattice.get(iX+dx, iY+dy, iZ+dz).computeDensity();
                    lattice.get(iX, iY, iZ).defineDensity(rho);
                }
            }
        }
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const{
        modified[0] = modif::dataStructure;
    }
private:
    plint dx;
    plint dy;
    plint dz;
};

void writeVTK(  MultiBlockLattice3D<T,DESCRIPTOR>& leftFluid,
                MultiBlockLattice3D<T,DESCRIPTOR>& rightFluid,
                T dx, T dt, plint iter)
{
    VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), 1.);
    vtkOut.writeData<float>(*computeDensity(leftFluid), "densityleft", 1.);
    vtkOut.writeData<3,float>(*computeVelocity(leftFluid), "velocityleft", dx/dt);
    vtkOut.writeData<float>(*computeDensity(rightFluid), "densityright", 1.);
    vtkOut.writeData<3,float>(*computeVelocity(rightFluid), "velocityright", dx/dt);
}

void readParameters(XMLreader const& document)
{
    document["geometry"]["averageInletVelocity"].read(averageInletVelocity);

    document["fluid"]["kinematicViscosity"].read(kinematicViscosity);
    document["fluid"]["density"].read(fluidDensity);
    document["fluid"]["volume"].read(volume);
    document["fluid"]["Re"].read(Re);
    document["fluid"]["G"].read(G);

    document["numerics"]["referenceResolution"].read(referenceResolution);
    document["numerics"]["uAveLB"].read(uAveLB);

    document["simulation"]["simTime"].read(simTime);
    document["simulation"]["maxIter"].read(maxIter);
    document["simulation"]["writeInterval"].read(writeInterval);
    document["simulation"]["epsilon"].read(epsilon);
    document["simulation"]["performOutput"].read(performOutput);
}

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");
    srand(global::mpi().getRank());

    string paramXmlFileName;
    try{
        global::argv(1).read(paramXmlFileName);
    }
    catch (PlbIOException& exception) {
        pcout << "Wrong parameter; the syntax is: "
        << (std::string)global::argv(0) << " parameter-input-file.xml" << std::endl;
        return -1;
    }

    // Read the parameter XML input file
    try {
        XMLreader document(paramXmlFileName);
        readParameters(paramXmlFileName);
    }
    catch (PlbIOException& exception) {
        pcout << "Error in input file " << paramXmlFileName
        << ": " << exception.what() << std::endl;
        return -1;
    }

    T dx = (1./referenceResolution)*0.0015;
    T dt = uAveLB / averageInletVelocity *dx;
    T nuLB = (uAveLB*referenceResolution) / Re;
    T tau = (3.*nuLB+0.5);
    T omega = 1./tau;

    if (performOutput) {
         pcout << "dx = " << dx << std::endl;
         pcout << "dt = " << dt << std::endl;
         pcout << "nuLB = " << nuLB << std::endl;
         pcout << "tau = " << tau << std::endl;
         pcout << "G = " <<  G << std::endl;
     }

     if(tau <= 0.5){
         pcout << "tau has invalid value." << std::endl;
         exit(1);
     }

    const plint nx   = referenceResolution/15+2;
    const plint ny   = referenceResolution+2;
    const plint nz   = referenceResolution/15+2;
    plint inletSize = referenceResolution/30;

    // Use regularized BGK dynamics to improve numerical stability (but note that
    //   BGK dynamics works well too).
    MultiBlockLattice3D<T, DESCRIPTOR> rightFluid (
            nx,ny,nz, new ExternalMomentRegularizedBGKdynamics<T, DESCRIPTOR>(omega) );
    MultiBlockLattice3D<T, DESCRIPTOR> leftFluid (
            nx,ny,nz, new ExternalMomentRegularizedBGKdynamics<T, DESCRIPTOR>(omega) );

    // Store a pointer to all lattices (two in the present application) in a vector to
    //   create the Shan/Chen coupling therm. The heavy fluid being at the first place
    //   in the vector, the coupling term is going to be executed at the end of the call
    //   to collideAndStream() or stream() for the heavy fluid.
    vector<MultiBlockLattice3D<T, DESCRIPTOR>* > blockLattices;
    blockLattices.push_back(&rightFluid);
    blockLattices.push_back(&leftFluid);

    // The argument "constOmegaValues" to the Shan/Chen processor is optional,
    //   and is used for efficiency reasons only. It tells the data processor
    //   that the relaxation times are constant, and that their inverse must be
    //   computed only once.
    std::vector<T> constOmegaValues;
    constOmegaValues.push_back(omega);
    constOmegaValues.push_back(omega);
    plint processorLevel = 1;
    integrateProcessingFunctional (
            new ShanChenMultiComponentProcessor3D<T,DESCRIPTOR>(G,constOmegaValues),
            Box3D(0,nx-1,1,ny-2,0,nz-1),
            blockLattices,
            processorLevel );

        // The setup is: periodicity along horizontal direction, bounce-back on top
        // and bottom. The upper half is initially filled with fluid 1 + random noise,
        // and the lower half with fluid 2. Only fluid 1 experiences a forces,
        // directed downwards.

        Box3D LeftWall(0,0, 0, ny-1, 0, nz-1);
        Box3D RightWall(nx-1,nx-1, 0, ny-1, 0, nz-1);
        Box3D FrontWall(1,nx-2,1,ny-2,nz-1,nz-1);
        Box3D BackWall(1,nx-2,1,ny-2,0,0);
        Box3D outletDomain(1,nx-2, 0,0, 0, nz-1);
        Box3D inlet1Domain(1,inletSize, ny-1,ny-1, 0, nz-1);
        Box3D inlet2Domain(inletSize+1, nx-2, ny-1, ny-1, 0, nz-1);

        pcout << "nx: " << nx << ", ny: " << ny << std::endl;

        // Bounce-back on bottom wall.
        defineDynamics(leftFluid, outletDomain, new VelocityBounceBack<T, DESCRIPTOR>((T)1., Array<T,3>((T) 0., -uAveLB, (T) 0.)) );
        defineDynamics(rightFluid, outletDomain, new VelocityBounceBack<T, DESCRIPTOR>((T)1., Array<T,3>((T) 0., -uAveLB, (T) 0.)) );
        // Bounce-back on top wall.
            // On the left inlet
        defineDynamics(leftFluid, inlet1Domain, new VelocityBounceBack<T, DESCRIPTOR>((T)1, Array<T,3>((T) 0., -uAveLB, (T) 0.)) );
        defineDynamics(rightFluid, inlet1Domain, new VelocityBounceBack<T, DESCRIPTOR>((T)1.e-4, Array<T,3>((T) 0., -uAveLB, (T) 0.)) );
            // On the right inlet
        defineDynamics(leftFluid, inlet2Domain, new VelocityBounceBack<T, DESCRIPTOR>((T)1.e-4, Array<T,3>((T) 0., -uAveLB, (T) 0.)) );
        defineDynamics(rightFluid, inlet2Domain, new VelocityBounceBack<T, DESCRIPTOR>((T)1, Array<T,3>((T) 0., -uAveLB, (T) 0.)) );
        // Define Bounceback left wall
        defineDynamics(leftFluid, LeftWall, new VelocityBounceBack<T, DESCRIPTOR>((T)1., Array<T,3>((T) 0., (T) 0., (T) 0.)) );
        defineDynamics(rightFluid, LeftWall, new VelocityBounceBack<T, DESCRIPTOR>((T)1.e-4, Array<T,3>((T) 0., (T) 0., (T) 0.)) );
        // Define Bounceback right wall
        defineDynamics(leftFluid, RightWall, new VelocityBounceBack<T, DESCRIPTOR>((T)1.e-4, Array<T,3>((T) 0., (T) 0., (T) 0.)) );
        defineDynamics(rightFluid, RightWall, new VelocityBounceBack<T, DESCRIPTOR>((T)1., Array<T,3>((T) 0., (T) 0., (T) 0.)) );
        // Define Bounceback front wall
        defineDynamics(leftFluid, FrontWall, new VelocityBounceBack<T, DESCRIPTOR>((T)1., Array<T,3>((T) 0., (T) 0., (T) 0.)) );
        defineDynamics(rightFluid, FrontWall, new VelocityBounceBack<T, DESCRIPTOR>((T)1.e-4, Array<T,3>((T) 0., (T) 0., (T) 0.)) );
        // Define Bounceback back wall
        defineDynamics(leftFluid, BackWall, new VelocityBounceBack<T, DESCRIPTOR>((T)1.e-4, Array<T,3>((T) 0., (T) 0., (T) 0.)) );
        defineDynamics(rightFluid, BackWall, new VelocityBounceBack<T, DESCRIPTOR>((T)1., Array<T,3>((T) 0., (T) 0., (T) 0.)) );


        // Initialize top layer.
        applyIndexed(rightFluid, Box3D(0, nx-1, 0, ny-1, 0, nz-1),
                     new TwoLayerInitializer<T,DESCRIPTOR>(nx, true) );
        // Initialize bottom layer.
        applyIndexed(leftFluid, Box3D(0, nx-1, 0, ny-1, 0, nz-1),
                     new TwoLayerInitializer<T,DESCRIPTOR>(nx, false) );

        leftFluid.initialize();
        rightFluid.initialize();

    T epsilon = 1e-2;
    T time = T();
    util::ValueTracer<T> velocityTracer1((T)1, ny, epsilon);
    util::ValueTracer<T> velocityTracer2((T)1, ny, epsilon);
    global::timer("Simulation").start();
    const plint convergence = 100;
    const plint statIter = 1000;
    plint iT = 0;

    pcout << "Starting simulation" << endl;
    // Main loop over time iterations.
    while(!velocityTracer1.hasConverged() || !velocityTracer2.hasConverged())
    {
        if (iT%writeInterval==0) {
            pcout << "Writing vtk at iteration: " << iT << std::endl;
            writeVTK(leftFluid, rightFluid, (T)1, (T)1, iT);
        }

        if (iT%convergence==0) {
            velocityTracer1.takeValue(computeAverageEnergy(leftFluid));
            velocityTracer2.takeValue(computeAverageEnergy(rightFluid));
        }

        // Execute lattice Boltzmann iteration
        applyProcessingFunctional(new SetDensity<T,DESCRIPTOR>(1,0,0), LeftWall, rightFluid);
        applyProcessingFunctional(new SetDensity<T,DESCRIPTOR>(1,0,0), LeftWall, leftFluid);
        applyProcessingFunctional(new SetDensity<T,DESCRIPTOR>(-1,0,0), RightWall, rightFluid);
        applyProcessingFunctional(new SetDensity<T,DESCRIPTOR>(-1,0,0), RightWall, leftFluid);
        applyProcessingFunctional(new SetDensity<T,DESCRIPTOR>(0,0,-1), FrontWall, rightFluid);
        applyProcessingFunctional(new SetDensity<T,DESCRIPTOR>(0,0,-1), FrontWall, leftFluid);
        applyProcessingFunctional(new SetDensity<T,DESCRIPTOR>(0,0,1), BackWall, rightFluid);
        applyProcessingFunctional(new SetDensity<T,DESCRIPTOR>(0,0,1), BackWall, leftFluid);
        applyProcessingFunctional(new SetDensity<T,DESCRIPTOR>(0,1,0), outletDomain, rightFluid);
        applyProcessingFunctional(new SetDensity<T,DESCRIPTOR>(0,1,0), outletDomain, leftFluid);
        leftFluid.collideAndStream();
        rightFluid.collideAndStream();

        if (iT%statIter==0) {
            pcout << "Average density fluid one = "
                  << getStoredAverageDensity<T>(rightFluid);
            pcout << ", average density fluid two = "
                  << getStoredAverageDensity<T>(leftFluid) << endl;
            pcout << "Convergence: " << velocityTracer1.hasConverged() << ", "
                    << velocityTracer2.hasConverged() << std::endl;
        }
        ++iT;
    }
    time = global::timer("Simulation").stop();
    pcout << "Simulation took: " << time << std::endl;
}
