#include "palabos2D.h"
#include "palabos2D.hh"   // include full template code

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "functions.h"
#include "newtonRaphson.h"
#include "trapeziumIntegration.h"

using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
#define DESCRIPTOR D2Q9Descriptor


/// Velocity on the parabolic Poiseuille profile
T poiseuilleVelocity(plint iY, CarreauFlowParam<T> const& parameters) {
    T y = (T)iY / parameters.getResolution();
    return 4.*parameters.getLatticeU() * (y-y*y);
}

/// Linearly decreasing pressure profile
T poiseuillePressure(plint iX, CarreauFlowParam<T> const& parameters) {
    T Lx = parameters.getNx()-1;
    T Ly = parameters.getNy()-1;
    return 8.*parameters.getLatticeNu0()*parameters.getLatticeU() / (Ly*Ly) * (Lx/(T)2-(T)iX);
}

/// Convert pressure to density according to ideal gas law
T poiseuilleDensity(plint iX, CarreauFlowParam<T> const& parameters) {
    return poiseuillePressure(iX,parameters)*DESCRIPTOR<T>::invCs2 + (T)1;
}

/// A functional, used to initialize the velocity for the boundary conditions
template<typename T>
class PoiseuilleVelocity {
public:
    PoiseuilleVelocity(CarreauFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    void operator()(plint iX, plint iY, Array<T,2>& u) const {
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = T();
    }
private:
    CarreauFlowParam<T> parameters;
};

/// A functional, used to initialize the density for the boundary conditions
template<typename T>
class PoiseuilleDensity {
public:
    PoiseuilleDensity(CarreauFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    T operator()(plint iX, plint iY) const {
        return poiseuilleDensity(iX,parameters);
    }
private:
    CarreauFlowParam<T> parameters;
};

// InstantiateDynamicsFunctional2D
template<typename T, template<typename U> class Descriptor>
class VelBCInitializer : public BoxProcessingFunctional2D_L<T,Descriptor> {
public:
    VelBCInitializer(CarreauFlowParam<T> const& parameters_)
        : parameters(parameters_)
    { }
    VelBCInitializer<T,Descriptor>* clone() const {
        return new VelBCInitializer<T,Descriptor>(*this);
    }
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice) {
        Array<T,2> u;
        T rho = (T)1;
        u[1] = T();

        for (plint iX=domain.x0; iX<=domain.x1; ++iX){
            rho = poiseuilleDensity(iX, parameters);
            for (plint iY=domain.y0; iY<=domain.y1; ++iY){
                //pcout << iY << std::endl;
                u[0] = poiseuilleVelocity(iY-1, parameters);
                pcout << "loc: (" <<iX <<", "<<iY<<")   rho: " << rho << ", speed: (" << u[0]<<", " << u[1] << ")" << endl;
                lattice.attributeDynamics(iX,iY, new VelocityBounceBack<T, DESCRIPTOR>(rho, u));
                iniCellAtEquilibrium(lattice.get(iX,iY), rho, u);
            }
        }
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const{
        modified[0] = modif::dataStructure;
    }
private:
    CarreauFlowParam<T> const& parameters;
};

/// A functional, used to create an initial condition for with zero velocity,
///   and linearly decreasing pressure.
template<typename T>
class PoiseuilleDensityAndVelocity {
public:
    PoiseuilleDensityAndVelocity(CarreauFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    void operator()(plint iX, plint iY, T& rho, Array<T,2>& u) const {
        rho = poiseuilleDensity(iX, parameters);
        //rho = (T) 1;
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = T();
        //pcout << "loc: (" <<iX <<", "<<iY<<")   rho: " << rho << ", pressure: " << poiseuillePressure(iX, parameters) << endl;
    }
private:
    CarreauFlowParam<T> parameters;
};

void definePoiseuilleGeometry( MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
                               CarreauFlowParam<T> const& parameters,
                               OnLatticeBoundaryCondition2D<T,DESCRIPTOR>& boundaryCondition,
                               T tol, plint maxIter)
{
    setCompositeDynamics (
            lattice,
            lattice.getBoundingBox(),
            new CarreauDynamics<T,DESCRIPTOR,1>(new NoDynamics<T,DESCRIPTOR>) );

    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();

    Box2D inlet(0, 0, 1, ny-2);
    Box2D outlet(nx-1, nx-1, 1, ny-2);
    Box2D bottomWall(0, nx-1, 0, 0);
    Box2D topWall(0, nx-1, ny-1, ny-1);

    // Bounce-back boundary conditions on the bottom and top walls
    defineDynamics(lattice, Box2D(0, nx-1, 0, 0), new BounceBack<T, DESCRIPTOR>((T) 1.) );
    defineDynamics(lattice, Box2D(0, nx-1, ny-1, ny-1), new BounceBack<T, DESCRIPTOR>((T) 1.) );

    // Initialize inlet velocity conditions
    defineDynamics(lattice, Box2D(0, 0, 1, ny-2), new VelocityBounceBack<T, DESCRIPTOR>((T)1., Array<T,2>(parameters.getLatticeU(), (T)0.)));
    // Initialize outlet velocity conditions
    defineDynamics(lattice, Box2D(nx-1, nx-1, 1, ny-2), new AntiBounceBack<T, DESCRIPTOR>((T)1., Array<T,2>((T)0., (T)0.)));

    // Initialize all cells at an equilibrium distribution, with a velocity and density
    //   value of the analytical Poiseuille solution.
    initializeAtEquilibrium (
           lattice, lattice.getBoundingBox(), (T) 1., Array<T,2>((T)0.,(T)0.) );

    pcout << "lattice totally initialized." << endl;
}

T computeRMSerror ( MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
                    MultiTensorField2D<T,2> &analyticalVelocity,
                    CarreauFlowParam<T> const& parameters )
{
    MultiTensorField2D<T,2> numericalVelocity(lattice);
    computeVelocity(lattice, numericalVelocity, lattice.getBoundingBox());

    // Divide by lattice velocity to normalize the error
    return 1./parameters.getLatticeU() *
    // Compute RMS difference between analytical and numerical solution
    std::sqrt( computeAverage( *computeNormSqr(
        *subtract(analyticalVelocity, numericalVelocity)
    ) ) );
}

void writeVTK(MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
            T dx, T dt, plint iter)
{
    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<2,float>(*computeVelocity(lattice), "vel", dx/dt);
    vtkOut.writeData<float>(*computeOmega(lattice), "omega", (T)1);
    vtkOut.writeData<float>(*computeDensity(lattice), "density", (T)1);
}

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");
    global::timer("simTime").start();

    if (argc != 2)
    {
        pcout << "Error N must be specified." << std::endl;
        exit(1);
    }

    const plint Nref = 100;
    const T uMaxRef = 0.01;
    const plint N = atoi(argv[1]);

    const T uMax = uMaxRef * (T)Nref / (T)N;

    CarreauFlowParam<T> parameters(
            (T) 0.01,      // physical uMax
            (T) 0.001,    // lattice uMax
            (T) 0.0641,    // Re
            (T) 0.0001,    // physical Length
            (T) 11, // Cu
            (T) 1.56e-6,   // NuInf (Only the case nuInf = 0 has been implemented for the velocity inlet bc)
            (T) 0.392,      // n
            (T) 0.644,      // a
             N,          // N
             5.,         // lattice lx
             1.          // lattice ly
    );

    const plint maxT     = 10000000;
    T LatticeU = parameters.getLatticeU();
    pcout << uMax << std::endl;
    pcout << LatticeU << std::endl;

    writeLogFile(parameters, "Carreau Poseuille Flow");

    global::CarreauParameters().setNu0(parameters.getLatticeNu0());
    global::CarreauParameters().setNuInf(parameters.getLatticeNuInf());
    global::CarreauParameters().setLambda(parameters.getLatticeLambda());
    global::CarreauParameters().setExponent(parameters.getExponent());
    global::CarreauParameters().setA(parameters.getA());

    pcout << "dt = " << parameters.getDeltaT() << ", dx = " << parameters.getDeltaX() << std::endl;
    pcout << "nu0 = " << parameters.getLatticeNu0() << ", nuInf = " << parameters.getLatticeNuInf() <<
    ", Lambda = " << parameters.getLatticeLambda() <<  ", exponent = " << parameters.getExponent() <<
    ", a = " << parameters.getA() << std::endl;
    pcout << "u = " << parameters.getLatticeU() << std::endl;

    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
            parameters.getNx(), parameters.getNy(),
            new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega0()) );

    defineDynamics(lattice, lattice.getBoundingBox(), new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega0()) );

    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();

    T tolAll = 1.0e-11;
    plint maxIterAll = 100;

    definePoiseuilleGeometry(lattice, parameters, *boundaryCondition, tolAll, maxIterAll);

    util::ValueTracer<T> converge(uMax,Nref,1.0e-3);

    T tIni = global::timer("simTime").stop();
    pcout << "time elapsed for iniGeometry:" << tIni << endl;
    global::timer("simTime").start();

    plint iT = 0;
    for (iT=0; iT<maxT; ++iT)
    {
        converge.takeValue(getStoredAverageEnergy(lattice),true);
        if (iT % 1000 == 0)
        {
            pcout << iT << " : Writing image." << endl;
            writeVTK(lattice, parameters.getDeltaX(), parameters.getDeltaT() ,iT);
        }

        if (converge.hasConverged())
        {
            pcout << "Simulation converged." << endl;
            break;
        }

        lattice.collideAndStream();
    }

    T tEnd = global::timer("simTime").stop();

    T totalTime = tEnd-tIni;
    T N1000 = lattice.getNx()/(T)1000;
    pcout << "N=" << N << endl;
    pcout << "number of processors: " << global::mpi().getSize() << endl;
    pcout << "simulation time: " << totalTime << endl;
    pcout << "total time: " << tEnd << endl;
    pcout << "total iterations: " << iT << endl;
    pcout << "Msus: " << N1000*N1000*(T)iT/totalTime << endl;

    delete boundaryCondition;
}
