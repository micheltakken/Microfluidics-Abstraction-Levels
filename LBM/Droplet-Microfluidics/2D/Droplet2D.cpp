#include "palabos2D.h"
#include "palabos2D.hh"
#include "palabos3D.h"
#include "palabos3D.hh"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::D2Q9Descriptor
#define DESCRIPTORSHAN descriptors::ShanChenD2Q9Descriptor

plint extraLayer            = 0.;   // Make the bounding box larger; for visualization purposes
                                    // only. For the simulation, it is ok to have extraLayer=0.
const plint blockSize       = 20;   // Zero means: no sparse representation.
const plint envelopeWidth   = 1;    // For standard BGK Dynamics
const plint extendedEnvelopeWidth = 2;  // Because the Guo off lattice boundary condition
                                        // needs 2-cell neighbor access.

T kinematicViscosity        = 0.;
T averageInletVelocity      = 0.;
plint referenceResolution   = 0.;

plint referenceDirection    = 0;
plint openingSortDirection  = 0;

T uAveLB                    = 0.;
T fluidDensity              = 0.;
T volume                    = 0.;
T Re                        = 0.;
T G                         = 0.;

T simTime = 0;
plint maxIter = 0;
plint writeInterval = 0;
plint maxLevel = 0;
T epsilon = 0;
bool performOutput = false;

TriangleSet<T>* triangleSet = 0;
T currentTime = 0;

// Structure which defines an 'opening'. The provided STL file contains holes
// that need to be defined as Inlets and Outlets.
template<typename T>
struct Opening {
    bool inlet;
    Array<T,3> point1;
    Array<T,3> point2;
};

std::vector<Opening<T>> openings;

// Initial condition: Water droplet immersed in oil in the inlet channel.
// This functional is going to be used as an argument to the function "applyIndexed",
// to setup the initial condition.
template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
class CopyVelocity : public BoxProcessingFunctional2D_LL<T1,Descriptor1,T2,Descriptor2>
{
public:
    CopyVelocity(bool droplet_)
        : droplet(droplet_)
        { }
    CopyVelocity<T1,Descriptor1,T2,Descriptor2>* clone() const {
        return new CopyVelocity<T1,Descriptor1,T2,Descriptor2>(*this);
    }
    virtual void process(Box2D domain, BlockLattice2D<T1,Descriptor1>& latticeFrom,
                        BlockLattice2D<T2,Descriptor2>& latticeTo) {
        Dot2D relativePosition = latticeFrom.getLocation();
        Dot2D offset = computeRelativeDisplacement(latticeFrom, latticeTo);
        Array<T,2> tempVel (0.,0.);
        Array<plint,2> center (40,103);
        plint radius2 = 100;
        T almostNoFluid = 1.e-4;
        T rho = 1.;
        T loc;
        for (plint iX=domain.x0; iX<=domain.x1; ++iX){
            for (plint iY=domain.y0; iY<=domain.y1; ++iY){
                loc = (iX+relativePosition.x-center[0])*(iX+relativePosition.x-center[0])
                    + (iY+relativePosition.y-center[1])*(iY+relativePosition.y-center[1]);
                latticeFrom.get(iX,iY).computeVelocity(tempVel);
                if (!droplet && (std::abs(center[0]-iX) >= 40 || std::abs(center[1]-iY) >= 18)){
                    rho = latticeFrom.get(iX,iY).computeDensity();
                } else if (droplet && std::abs(center[0]-iX) < 40 && std::abs(center[1]-iY) < 18){
                    rho = latticeFrom.get(iX,iY).computeDensity();
                } else {
                    rho = almostNoFluid;
                }
                //rho = latticeFrom.get(iX,iY).computeDensity();
                iniCellAtEquilibrium(latticeTo.get(iX+offset.x,iY+offset.y), rho, tempVel);
            }
        }
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const{
        for (pluint iBlock=0; iBlock<modified.size(); ++iBlock) {
            modified[iBlock] = modif::staticVariables;
        }
    }
private:
    bool droplet;
};

// Initiate 2D field from STL file
template<typename T>
class domainInitializer3D : public BoxProcessingFunctional3D_S<T> {
public:
    domainInitializer3D(ScalarField2D<T>* write_)
        : write(write_)
    { }
    domainInitializer3D<T>* clone() const {
        return new domainInitializer3D<T>(*this);
    }
    virtual void process(Box3D domain, ScalarField3D<T>& from) {
        PLB_ASSERT(domain.z0 == domain.z1);
        Dot3D rel = from.getLocation();
        for (plint iY=domain.y0; iY<=domain.y1; ++iY){
            for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
                write->get(iX+rel.x, iY+rel.y) = from.get(iX,iY,domain.z0);
            }
        }
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const{
        modified[0] = modif::dataStructure;
    }
private:
    ScalarField2D<T>* write;
};

// Initiate 2D field from STL file
template<typename T>
class domainInitializer2D : public BoxProcessingFunctional2D_S<T> {
public:
    domainInitializer2D(ScalarField2D<T>* write_)
        : write(write_)
    { }
    domainInitializer2D<T>* clone() const {
        return new domainInitializer2D<T>(*this);
    }
    virtual void process(Box2D domain, ScalarField2D<T>& to) {
        Dot2D rel = to.getLocation();
        for (plint iY=domain.y0; iY<=domain.y1; ++iY){
            for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
                to.get(iX,iY) = write->get(iX+rel.x, iY+rel.y);
                //pcout << to.get(iX,iY);
            }
            //pcout << std::endl;
        }
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const{
        modified[0] = modif::dataStructure;
    }
private:
    ScalarField2D<T>* write;
};

// InstantiateDynamicsFunctional2D
template<typename T, template<typename U> class Descriptor>
class VelBCInitializer : public BoxProcessingFunctional2D_L<T,Descriptor> {
public:
    VelBCInitializer(bool droplet_, Array<T,2> velocity_)
        : droplet(droplet_),
        u(velocity_)
    { }
    VelBCInitializer<T,Descriptor>* clone() const {
        return new VelBCInitializer<T,Descriptor>(*this);
    }
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice) {
        T almostNoFluid = 1.e-4;
        T rho = (T)1;

        if ( droplet ) {
            rho = almostNoFluid;
            for (plint iX=domain.x0; iX<=domain.x1; ++iX){
                for (plint iY=domain.y0; iY<=domain.y1; ++iY){
                    lattice.attributeDynamics(iX,iY, new VelocityBounceBack<T, Descriptor>(rho, u));
                }
            }
        }
        if ( !droplet ){
            for (plint iX=domain.x0; iX<=domain.x1; ++iX){
                for (plint iY=domain.y0; iY<=domain.y1; ++iY){
                    lattice.attributeDynamics(iX,iY, new VelocityBounceBack<T, Descriptor>(rho, u));
                }
            }
        }
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const{
        modified[0] = modif::dataStructure;
    }
private:
    bool droplet;
    Array<T,2> u;
};

template<typename T, template<typename U> class Descriptor>
class PreBCInitializer : public BoxProcessingFunctional2D_L<T,Descriptor> {
public:
    PreBCInitializer(bool droplet_, T pressure_)
        : droplet(droplet_),
        pressure(pressure_)
    { }
    PreBCInitializer<T,Descriptor>* clone() const {
        return new PreBCInitializer<T,Descriptor>(*this);
    }
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice) {
        T rho = (T)1;
        T almostNoFluid = 1.e-4;

        if ( droplet ){
            rho = almostNoFluid;
            for (plint iX=domain.x0; iX<=domain.x1; ++iX){
                for (plint iY=domain.y0; iY<=domain.y1; ++iY){
                    lattice.attributeDynamics(iX,iY, new AntiBounceBack<T, Descriptor>(rho, Array<T,2>((T) 0., (T) 0.)));
                }
            }
        }
        if ( !droplet ){
            rho = pressure*DESCRIPTOR<T>::invCs2 + (T)1;
            for (plint iX=domain.x0; iX<=domain.x1; ++iX){
                for (plint iY=domain.y0; iY<=domain.y1; ++iY){
                    lattice.attributeDynamics(iX,iY, new AntiBounceBack<T, Descriptor>(rho, Array<T,2>((T) 0., (T) 0.)));
                }
            }
        }
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const{
        modified[0] = modif::dataStructure;
    }
private:
    bool droplet;
    T pressure;
};

void caseSetup (
        MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
        MultiScalarField2D<int>& voxelizedDomain)
{
        // Create the initial condition
        initializeAtEquilibrium(
            lattice, lattice.getBoundingBox(), (T) 1., Array<T,2>((T) 0., (T) 0.));
}

void caseSetup (
        MultiBlockLattice2D<T,DESCRIPTORSHAN>& lightFluid,
        MultiBlockLattice2D<T,DESCRIPTORSHAN>& heavyFluid,
        MultiScalarField2D<int>& voxelizedDomain)
{
        // Create the initial condition
        initializeAtEquilibrium(
            lightFluid, lightFluid.getBoundingBox(), (T) 1., Array<T,2>((T) 0., (T) 0.));
        initializeAtEquilibrium(
            heavyFluid, heavyFluid.getBoundingBox(), (T) 1., Array<T,2>((T) 0., (T) 0.));
}

/// Produce a GIF snapshot of the velocity-norm.
void writeGif(MultiBlockLattice2D<T,DESCRIPTOR>& lattice, plint iter)
{
    const plint imSize = 600;

    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif(createFileName("u", iter, 6),
                               *computeVelocityNorm(lattice),
                               imSize, imSize );
}

/// Write the full velocity and the velocity-norm into a VTK file.
void writeVTK(  MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
                T dx, T dt, plint iter)
{
    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), 1.);
    vtkOut.writeData<float>(*computeDensity(lattice), "densityLight", 1.);
    vtkOut.writeData<2,float>(*computeVelocity(lattice), "velocityLight", dx/dt);
}

void writeVTK(MultiBlockLattice2D<T,DESCRIPTORSHAN>& lightFluid,
                MultiBlockLattice2D<T,DESCRIPTORSHAN>& heavyFluid,
                T dx, T dt, plint iter)
{
    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), 1.);
    vtkOut.writeData<float>(*computeDensity(lightFluid), "densityLight", 1.);
    vtkOut.writeData<float>(*computeDensity(heavyFluid), "densityHeavy", 1.);
    vtkOut.writeData<2,float>(*computeVelocity(lightFluid), "velocityLight", dx/dt);
    vtkOut.writeData<2,float>(*computeVelocity(heavyFluid), "velocityHeavy", dx/dt);
}

void readParameters(XMLreader const& document)
{
    std::string meshFileName;
    std::vector<std::string> openingType;
    document["geometry"]["mesh"].read(meshFileName);
    document["geometry"]["averageInletVelocity"].read(averageInletVelocity);
    document["geometry"]["openings"]["sortDirection"].read(openingSortDirection);
    document["geometry"]["openings"]["type"].read(openingType);

    document["fluid"]["kinematicViscosity"].read(kinematicViscosity);
    document["fluid"]["density"].read(fluidDensity);
    document["fluid"]["volume"].read(volume);
    document["fluid"]["Re"].read(Re);
    document["fluid"]["G"].read(G);

    document["numerics"]["referenceDirection"].read(referenceDirection);
    document["numerics"]["referenceResolution"].read(referenceResolution);
    document["numerics"]["uAveLB"].read(uAveLB);

    document["simulation"]["simTime"].read(simTime);
    document["simulation"]["maxIter"].read(maxIter);
    document["simulation"]["writeInterval"].read(writeInterval);
    document["simulation"]["maxLevel"].read(maxLevel);
    document["simulation"]["epsilon"].read(epsilon);
    document["simulation"]["performOutput"].read(performOutput);

    triangleSet = new TriangleSet<T>(meshFileName, DBL);
    openings.resize(openingType.size());
}

std::unique_ptr<MultiBlockLattice2D<T,DESCRIPTOR> > run (
    plint level, MultiBlockLattice2D<T,DESCRIPTOR>* iniVal=0 )
{
    plint margin = 1; // Extra margin of allocated cells around the obstacle.
    plint borderWidth = 1; // Because the Guo boundary condition acts in a one-cell
                           // layer. Requirement: margin >= borderwidth

   // The resolution is doubled at each coordinate direction with the increase of the
   //   resolution level by one. The parameter ``referenceResolution'' is by definition
   //   the resolution at grid refinement level 0.
   plint resolution = referenceResolution * util::twoToThePower(level);

   DEFscaledMesh<T>* defMesh =
       new DEFscaledMesh<T>(*triangleSet, referenceResolution, referenceDirection, margin, extraLayer);
   TriangleBoundary3D<T> boundary(*defMesh);
   delete defMesh;
   boundary.getMesh().inflate();

   T dx = boundary.getDx();
   T dt = uAveLB / averageInletVelocity *dx;
   T nuLB = (uAveLB*resolution) / Re;
   T tau = (3.*nuLB+0.5);
   T omega = 1./tau;

   Array<T,3>location(boundary.getPhysicalLocation());

   if (performOutput) {
        pcout << "dx = " << dx << std::endl;
        pcout << "dt = " << dt << std::endl;
        pcout << "nuLB = " << nuLB << std::endl;
        pcout << "tau = " << tau << std::endl;
    }

    if(tau <= 0.5){
        pcout << "tau has invalid value." << std::endl;
        exit(1);
    }

    const int flowType = voxelFlag::inside;
    VoxelizedDomain3D<T> voxelizedDomain (
        boundary, flowType, extraLayer, borderWidth, extendedEnvelopeWidth, blockSize);

    MultiScalarField3D<int> flagMatrix3D((MultiBlock3D&)voxelizedDomain.getVoxelMatrix());
    setToConstant(flagMatrix3D, voxelizedDomain.getVoxelMatrix(),
                    voxelFlag::inside, flagMatrix3D.getBoundingBox(), 1);
    setToConstant(flagMatrix3D, voxelizedDomain.getVoxelMatrix(),
                    voxelFlag::innerBorder, flagMatrix3D.getBoundingBox(), 1);
    setToConstant(flagMatrix3D, voxelizedDomain.getVoxelMatrix(),
                    voxelFlag::outerBorder, flagMatrix3D.getBoundingBox(), 2);

    ScalarField2D<int> temp(flagMatrix3D.getBoundingBox().getNx(),
                            flagMatrix3D.getBoundingBox().getNy(),
                            0);

    MultiScalarField2D<int> flagMatrix2D(
        flagMatrix3D.getBoundingBox().getNx(),
        flagMatrix3D.getBoundingBox().getNy(),
        0);

    Box3D Domain3D( 0, flagMatrix3D.getBoundingBox().getNx(),
                        0, flagMatrix3D.getBoundingBox().getNy(),
                        util::roundToInt(flagMatrix3D.getBoundingBox().getNz()/2),
                        util::roundToInt(flagMatrix3D.getBoundingBox().getNz()/2));
    Box2D Domain2D( 0, flagMatrix2D.getBoundingBox().getNx(),
                        0, flagMatrix2D.getBoundingBox().getNy());

    Dynamics<T,DESCRIPTOR>* dynamics = new NoDynamics<T, DESCRIPTOR>;
    std::unique_ptr<MultiBlockLattice2D<T,DESCRIPTOR>> lattice
        = generateMultiBlockLattice(Domain2D, dynamics, envelopeWidth);
    lattice->toggleInternalStatistics(false);

    openings[0].point1 = Array<T,3> (0.00025, -0.0002, 0.);
    openings[0].point2 = Array<T,3> (0.00035, -0.0002, 0.0001);
    openings[1].point1 = Array<T,3> (0., 0., 0.);
    openings[1].point2 = Array<T,3> (0., 0.0001, 0.0001);
    openings[2].point1 = Array<T,3> (0.00025, 0.0003, 0.);
    openings[2].point2 = Array<T,3> (0.00035, 0.0003, 0.0001);

    Array<T,3> inletPos1((openings[1].point1 - location)/dx);
    Array<T,3> inletPos2((openings[1].point2 - location)/dx);
    Array<T,3> outlet1Pos1((openings[0].point1 - location)/dx);
    Array<T,3> outlet1Pos2((openings[0].point2 - location)/dx);
    Array<T,3> outlet2Pos1((openings[2].point1 - location)/dx);
    Array<T,3> outlet2Pos2((openings[2].point2 - location)/dx);

    Box2D inletDomain(  util::roundToInt(inletPos1[0])-1, util::roundToInt(inletPos2[0])-1,
                        util::roundToInt(inletPos1[1]), util::roundToInt(inletPos2[1]));
    Box2D outlet1Domain(  util::roundToInt(outlet1Pos1[0]), util::roundToInt(outlet1Pos2[0]),
                        util::roundToInt(outlet1Pos1[1])-1, util::roundToInt(outlet1Pos2[1])-1);
    Box2D outlet2Domain(  util::roundToInt(outlet2Pos1[0]), util::roundToInt(outlet2Pos2[0]),
                        util::roundToInt(outlet2Pos1[1])+1, util::roundToInt(outlet2Pos2[1])+1);

    applyProcessingFunctional(new domainInitializer3D<int>(&temp), Domain3D, flagMatrix3D);
    setToConstant(temp, inletDomain, 3);
    setToConstant(temp, outlet1Domain, 3);
    setToConstant(temp, outlet2Domain, 3);
    applyProcessingFunctional(new domainInitializer2D<int>(&temp), Domain2D, flagMatrix2D);

    // Set initial conditions and boundary conditions
    defineDynamics(*lattice, flagMatrix2D, lattice->getBoundingBox(), new BGKdynamics<T, DESCRIPTOR>(omega), 1);
    defineDynamics(*lattice, flagMatrix2D, lattice->getBoundingBox(), new BounceBack<T,DESCRIPTOR>(1.), 2);

    applyProcessingFunctional(new VelBCInitializer<T,DESCRIPTOR>(false, Array<T,2>(uAveLB, (T) 0.)), inletDomain, *lattice);
    applyProcessingFunctional(new PreBCInitializer<T,DESCRIPTOR>(false, -1.0e-3), outlet1Domain, *lattice);
    applyProcessingFunctional(new PreBCInitializer<T,DESCRIPTOR>(false, -2.0e-3), outlet2Domain, *lattice);

    caseSetup(*lattice, flagMatrix2D);
    lattice->initialize();
    if(iniVal) {
        Box2D toDomain(lattice->getBoundingBox());
        Box2D fromDomain(toDomain.shift(margin,margin));
        copy(*iniVal, fromDomain, *lattice, toDomain, modif::staticVariables);
    }

    // The ValueTracer is needed to check when a chosen quantity (in our case the average energy)
    //   has converged, so to conclude that steady state has been reached for the specific grid
    //   refinement level and stop the simulation.
    plint convergenceIter=20;
    util::ValueTracer<T> velocityTracer(0.05*convergenceIter, resolution, epsilon);
    global::timer("iteration").restart();
    plint i = util::roundToInt(currentTime/dt);
    lattice->resetTime(i);
    bool checkForErrors = true;

    // Collision and streaming iterations for lattice.
    while(!velocityTracer.hasConverged() && currentTime<simTime)
    {
        if (i%100==0 && performOutput) {
            pcout << "Iteration " << i << "; \t"
                << "T= " << currentTime << "; \t"
                << "Average energy: "
                << computeAverageEnergy(*lattice)*util::sqr(dx/dt) << std::endl;
                //writeVTK(*lattice, dx, dt, i);
        }
        if (i%convergenceIter==0) {
            velocityTracer.takeValue(computeAverageEnergy(*lattice));
        }

        lattice->collideAndStream();
        if (checkForErrors) {
            abortIfErrorsOccurred();
            checkForErrors = false;
        }

        ++i;
        currentTime = i*dt;
        if (i>1000){
            //exit(0);
        }
    }

    // Setup of Shan-Chen droplet problem
    if (level == maxLevel) {
        if (performOutput) {
            pcout << "dx=" << dx << std::endl;
            pcout << "dt=" << dt << std::endl;
            pcout << getMultiBlockInfo(voxelizedDomain.getVoxelMatrix()) << std::endl;
        }
        std::unique_ptr<MultiBlockLattice2D<T,DESCRIPTORSHAN>> lightFluid
                = generateMultiBlockLattice<T,DESCRIPTORSHAN> (
                        Domain2D, new NoDynamics<T, DESCRIPTORSHAN>, envelopeWidth);

        std::unique_ptr<MultiBlockLattice2D<T,DESCRIPTORSHAN>> heavyFluid
                = generateMultiBlockLattice<T,DESCRIPTORSHAN> (
                        Domain2D, new NoDynamics<T, DESCRIPTORSHAN>, envelopeWidth);

        // Set background conditions and domain
        // Set initial conditions and boundary conditions
        defineDynamics(*lightFluid, flagMatrix2D, lattice->getBoundingBox(), new ExternalMomentRegularizedBGKdynamics<T,DESCRIPTORSHAN>(omega), 1);
        defineDynamics(*lightFluid, flagMatrix2D, lattice->getBoundingBox(), new BounceBack<T,DESCRIPTORSHAN>(1.), 2);
        defineDynamics(*heavyFluid, flagMatrix2D, lattice->getBoundingBox(), new ExternalMomentRegularizedBGKdynamics<T,DESCRIPTORSHAN>(omega), 1);
        defineDynamics(*heavyFluid, flagMatrix2D, lattice->getBoundingBox(), new BounceBack<T,DESCRIPTORSHAN>(0.), 2);

        applyProcessingFunctional(new VelBCInitializer<T,DESCRIPTORSHAN>(false, Array<T,2>(uAveLB, (T) 0.)), inletDomain, *lightFluid);
        applyProcessingFunctional(new PreBCInitializer<T,DESCRIPTORSHAN>(false, -1.0e-3), outlet1Domain, *lightFluid);
        applyProcessingFunctional(new PreBCInitializer<T,DESCRIPTORSHAN>(false, -2.0e-3), outlet2Domain, *lightFluid);
        applyProcessingFunctional(new VelBCInitializer<T,DESCRIPTORSHAN>(true, Array<T,2>(uAveLB, (T) 0.)), inletDomain, *heavyFluid);
        applyProcessingFunctional(new PreBCInitializer<T,DESCRIPTORSHAN>(true, -1.0e-3), outlet1Domain, *heavyFluid);
        applyProcessingFunctional(new PreBCInitializer<T,DESCRIPTORSHAN>(true, -2.0e-3), outlet2Domain, *heavyFluid);

        vector<MultiBlockLattice2D<T,DESCRIPTORSHAN>* > blockLattices;
        blockLattices.push_back(heavyFluid.get());
        blockLattices.push_back(lightFluid.get());
        plint processorLevel = 1;

        integrateProcessingFunctional (
            new ShanChenMultiComponentProcessor2D<T,DESCRIPTORSHAN>(G),
            Domain2D, blockLattices, processorLevel );

        //caseSetup
        applyProcessingFunctional(new CopyVelocity<T,DESCRIPTOR,T,DESCRIPTORSHAN>(false), lattice->getBoundingBox(), *lattice, *lightFluid);
        applyProcessingFunctional(new CopyVelocity<T,DESCRIPTOR,T,DESCRIPTORSHAN>(true), lattice->getBoundingBox(), *lattice, *heavyFluid);

        lightFluid->initialize();
        heavyFluid->initialize();

        pcout << std::endl << "Starting main sequence for droplet simulation." << std::endl;
        // Main loop for Droplet flow
        for (plint iT=0; iT<maxIter; ++iT) {
            // Write vtk file every writeInterval timeseps
            if (iT%writeInterval==0){
                pcout << "Writing vtk file at iT = " << iT << endl;
                writeVTK(*lightFluid, *heavyFluid, dx, dt, iT);
            }
            // Execute lattice Boltzmann iteration
            lightFluid->collideAndStream();
            heavyFluid->collideAndStream();
        }
        exit(0);
    }
    return lattice;
}

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);

    global::directories().setOutputDir("./tmp/");
    global::IOpolicy().activateParallelIO(true);

    string paramXmlFileName;
    try{
        global::argv(1).read(paramXmlFileName);
    }
    catch (PlbIOException& exception) {
        pcout << "Wrong parameter; the syntax is: "
        << (std::string)global::argv(0) << " parameter-input-file.xml" << std::endl;
        return -1;
    }

    try {
        XMLreader document(paramXmlFileName);
        readParameters(paramXmlFileName);
    }
    catch (PlbIOException& exception) {
        pcout << "Error in input file " << paramXmlFileName
        << ": " << exception.what() << std::endl;
        return -1;
    }

    global::timer("global").start();
    plint iniLevel=0;
    std::unique_ptr<MultiBlockLattice2D<T,DESCRIPTOR> > iniConditionLattice(nullptr);
    try{
        // Loop for initial velocity and density conditions
        for (plint level=iniLevel; level<=maxLevel; ++level) {
            pcout << std::endl << "Running new simulation at level " << level << std::endl;
            std::unique_ptr<MultiBlockLattice2D<T,DESCRIPTOR>> convergedLattice (
                run(level, iniConditionLattice.get()) );
            if (level != maxLevel) {
                plint dxScale = -1;
                plint dtScale = -2;
                iniConditionLattice = std::unique_ptr<MultiBlockLattice2D<T,DESCRIPTOR> > (
                    refine(*convergedLattice, dxScale, dtScale, new BGKdynamics<T,DESCRIPTOR>(1.)) );
            }
        }
    }
    catch(PlbException& exception) {
        pcout << exception.what() << std::endl;
        return -1;
    }
}
