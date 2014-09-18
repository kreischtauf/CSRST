/***************************************************************************
                               main.cpp
                         -------------------
    begin                : Wed Sep 30 2009
    copyright            : (C) 2009 by Christof Kraus
    email                : christof.kraus-csrst@my.mail.de
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <math.h>
#include <float.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>

#include <boost/mpi/environment.hpp>

#include "defs.hh"
#include "utils.hh"

#include "Commands.hh"
#include "Physics.hh"
#include "Distribution.hh"
#include "PartBunch.hh"
#include "TETMSolver.hh"
#include "PoissonSolver.hh"
#include "SCBoundaryCondition.hh"
#include "Bend.hh"
#include "Communicator.hh"

std::ofstream dbg;
bool dbg_output;
#define DBGOUT dbg << "main.cpp: " << __LINE__ << "\t"

// Create an ArgTraits for the 2D vector type that declares it to be
// of string like type
namespace TCLAP {
    template<>
    struct ArgTraits<Vect2D> {
        typedef StringLike ValueCategory;
    };
}

void copyFields(VField_Edge_t & toE,
                VField_Cell_t & toH,
                VField_Edge_t & toJ,
                const VField_Edge_t & fromE,
                const VField_Cell_t & fromH,
                const VField_Edge_t & fromJ,
                const Vektor<int, DIM> & start,
                const bool & all2all,
                const bool & copyCurrent);
void exchangeBoundaries(VField_Edge_t & surroundingE,
                        VField_Cell_t & surroundingH,
                        VField_Edge_t & embeddedE,
                        VField_Cell_t & embeddedH,
                        const Vektor<int, DIM> & start);
void moveMesh(VField_Edge_t & miniMeshE,
              VField_Cell_t & miniMeshH,
              VField_Edge_t & miniMeshJ,
              VField_Edge_t & EFD,
              VField_Cell_t & HFD,
              NDIndex<DIM> & ddim,
              const NDIndex<DIM> & newGlobalPos);

int main(int argc, char *argv[]){
    Ippl ippl(argc, argv);
    boost::mpi::environment env(argc, argv);
    Communicator::initialize();

    boost::mpi::communicator world;

    Inform msg(argv[0]);
    Inform msg2all(argv[0], INFORM_ALL_NODES);

    msg << "\n"
        << CSRST_PACKAGE_NAME << " " << VERSION_MAJOR << "." << VERSION_MINOR << " (" << VERSION_PATCH << ")"
        << "\n";

    const std::string defaultFileName("Output.h5part");
    bool inStraightSectionBeforeBend = true;
    bool inBend = false;
    const size_t maxNumPoissonIterations = 500;
    const double relativeConvergenceTolerance = 1e-12;

    // declartion of the electric and the magnetic field. The memory is allocated in a separated command
    VField_Edge_t *EFD = 0; /* Pointer to the electric field */
    VField_Cell_t *HFD = 0; /* Pointer to the magnetic field */
    VField_Edge_t *JFD = 0; /* Pointer to the charge current field */
    VField_t *AlphaE = 0; /* Pointer to the electric components of the parameters for PML */
    VField_t *AlphaH = 0; /* Pointer to the magnetic components of the parameters for PML */

    double currentSimulatedTime = 0.0;

    Timings::TimerRef ipplToyFDTDTimer = Timings::getTimer("ipplToyFDTDTimer");
    Timings::startTimer(ipplToyFDTDTimer);
    Timings::TimerRef miniMeshTimer = Timings::getTimer("move mini mesh");
    Timings::TimerRef othersTimer = Timings::getTimer("other stuff");
    Timings::startTimer(othersTimer);
    // command line arguments parsing
    TCLAP::CmdLine cmdl("MaxwADIs: a 2D Maxwell simulation tool using an implicit ADI scheme", ' ', VERSION_STRING);
    TCLAP::ValueArg<string> restart   ("r", "restart", "file from which to restart", false, "", "String", cmdl);
    TCLAP::ValueArg<int> startstep    ("n", "startstep", "step to start from", false, 0, "Integer", cmdl);
    TCLAP::ValueArg<int> numsteps     ("N", "numsteps", "number of steps to perform", false, 10, "Integer", cmdl);
    TCLAP::ValueArg<int> numgpx       ("x", "numgridpx", "number of gridpoints in x-direction", false, 100, "Integer", cmdl);
    TCLAP::ValueArg<int> numgpy       ("y", "numgridpy", "number of gridpoints in y-direction", false, 100, "Integer", cmdl);
    TCLAP::ValueArg<int> plotModulus  ("m", "plotmodulus", "plot fields every pm steps", false, 100000, "Integer", cmdl);
    // TCLAP::ValueArg<int> saveModulus  ("M", "savemodulus", "save fields and particles every M steps", false, 100000, "Integer", cmdl);
    TCLAP::ValueArg<int> statModulus  ("", "statmodulus", "save bunch statistics every M steps", false, 100000, "Integer", cmdl);
    TCLAP::ValueArg<int> lineModulus  ("", "linemodulus", "save EM fields along line every M steps", false, 100000, "Integer", cmdl);
    TCLAP::ValueArg<int> rebModulus   ("", "reballance", "reballance the particles every B steps", false, 100000, "Integer", cmdl);
    TCLAP::ValueArg<double> phi       ("P","phi","angle of deflection", false, 0.0, "Double", cmdl);
    TCLAP::ValueArg<double> width     ("w","width","width of vacuum chamber", false, 0.10, "Double", cmdl);
    TCLAP::ValueArg<double> lenAfter  ("l","length","length of drift after bend", false, 0.10, "Double", cmdl);
    TCLAP::ValueArg<double> lenBefore ("","lengthbefore","length of drift before bend", false, 0.10, "Double", cmdl);
    TCLAP::ValueArg<double> ekin      ("E","Ekin","kinetic energy [MeV]", false, 250.0, "Double", cmdl);
    TCLAP::ValueArg<double> bzext     ("B","Bz","transverse magnetic field [T]", false, 0.10, "Double", cmdl);
    TCLAP::ValueArg<double> tcharge   ("Q","Qtotal","total charge [nC]", false, -1.0, "Double", cmdl);
    TCLAP::ValueArg<double> mean_Rx   ("", "mean_Rx", "initial x-position of bunch", false, 0.05, "Double", cmdl);
    TCLAP::ValueArg<double> mean_Ry   ("", "mean_Ry", "initial y-position of bunch", false, 0.0, "Double", cmdl);
    TCLAP::ValueArg<double> sigma_x   ("s", "sigma_x", "standard deviation of initial distribution", false, 0.01, "Double", cmdl);
    TCLAP::ValueArg<double> sigma_p   ("S", "sigma_p", "standard deviation of initial momenta", false, 1.0, "Double", cmdl);
    TCLAP::ValueArg<int> numparts     ("", "Np", "number of particles", false, 50000, "Integer", cmdl);
    TCLAP::ValueArg<string> filename  ("", "h5file", "read/write data from/to file", false, defaultFileName.c_str(), "String", cmdl);
    TCLAP::SwitchArg without_initials ("I", "no_initials", "don't calculate the initial conditions of the em fields", cmdl, false);
    TCLAP::SwitchArg debugoutput      ("d", "debug", "debug output", cmdl, false);
    try {
	cmdl.parse(argc, argv);
    } catch(std::exception &e) {
	std::cout << e.what() << std::endl;
	return EXIT_FAILURE;
    }

    std::string effectiveFileName = filename.getValue();
    // h5_int32_t saveFlag = H5_O_RDWR;

    LoadH5Part *lhpc = 0;
    if (restart.getValue() != "") {
        effectiveFileName = restart.getValue();
        // saveFlag = H5_O_APPEND;

        lhpc = new LoadH5Part(restart.getValue(), cmdl);
    }

    const int startAtStep = startstep.getValue();
    const int maximumIterations = numsteps.getValue();
    const int Nx = numgpx.getValue();
    const int Ny = numgpy.getValue();
    dbg_output = debugoutput.getValue();

    if (dbg_output) {
        // one debug file per core
        char buffer[18];
        sprintf(buffer, "debug_%05d.txt", Ippl::myNode());
        dbg.open(buffer);
    }

    AllocateFields afc(Nx, Ny);
    afc.execute(EFD, HFD, JFD, AlphaE, AlphaH, msg);

    Bend geo(lenBefore.getValue(),
             lenAfter.getValue(),
             width.getValue(),
             phi.getValue(),
             ekin.getValue(),
             bzext.getValue());

    geo.Init(*AlphaE,
             *AlphaH);

    // space and time differentials (dimensions of a cell) in meters / seconds
    Vector_t hr((*EFD).get_mesh().get_meshSpacing(0),
                (*EFD).get_mesh().get_meshSpacing(1));
    (*HFD).get_mesh().set_meshSpacing(&(hr(0)));
    (*HFD).get_mesh().set_origin((*EFD).get_mesh().get_origin());

    const double dt = hr(0) / Physics::c;

    // print out some simulation parameters:
    msg << "plot modulus = " <<  plotModulus.getValue() << "\n"
        << "meshing parameters: " << "\n"
        << "cells " <<  Nx << "  " <<  Ny << "\n"
        << "dx: " <<  hr(0) << ", dy: " <<  hr(1) << ", dt: " << dt << "\n"
        << "meter^2 simulation region " << Nx * hr(0) << " x " <<  Ny * hr(1) << "\n"
        << endl;

    SCBoundaryCondition<Bend, TETMSolver> bc(2, geo);

    Distribution distro(cmdl, *EFD, lhpc);
    NDIndex<DIM> ddim = distro.getDomain(*EFD);
    for (unsigned int d = 0; d < DIM; ++ d) {
        ddim[d] = Index(ddim[d].min() - 5, ddim[d].max() + 5);
    }
    Vector_t origin = (*EFD).get_mesh().get_origin();
    Vector_t miniMeshOrigin(origin(0) + ddim[0].min() * hr(0),
                            origin(1) + ddim[1].min() * hr(1));
    Index mmII(0, ddim[0].length());
    Index mmJJ(0, ddim[1].length());
    Mesh_t miniMesh(Index(ddim[0].length()), Index(ddim[1].length()), &hr[0], miniMeshOrigin);
    Mesh_t miniMeshCell(Index(ddim[0].length() + 1), Index(ddim[1].length() + 1), &hr[0], miniMeshOrigin);
    e_dim_tag decomp[] = {PARALLEL, PARALLEL, SERIAL};
    FieldLayout_Edge_t miniMeshEdgeFL(miniMesh, decomp);
    FieldLayout_Cell_t miniMeshCellFL(miniMeshCell, decomp);

    VField_Edge_t miniMeshE(miniMesh, miniMeshEdgeFL, GuardCellSizes<DIM>(GUARDCELLSIZE));
    VField_Cell_t miniMeshH(miniMeshCell, miniMeshCellFL, GuardCellSizes<DIM>(GUARDCELLSIZE));
    VField_Edge_t miniMeshJ(miniMesh, miniMeshEdgeFL, GuardCellSizes<DIM>(GUARDCELLSIZE));

    Timings::stopTimer(othersTimer);
    PartBunch bunch(miniMeshEdgeFL,
                    miniMeshCellFL,
                    miniMesh,
                    distro,
                    geo,
                    -1e-9 / Physics::m_e);

    if (restart.getValue() != "") {
        lhpc->restoreFields(*EFD, *HFD, *JFD, currentSimulatedTime);
    } else if (!without_initials.getValue()){
        bunch.commNumLocalParticles();

        PoissonSolver poisson(bunch,
                              geo,
                              bc.getBoundaryCells(),
                              EFD->get_mesh(),
                              EFD->getLayout(),
                              PBC,
                              PEC);

        poisson.computeField(*EFD,
                             *HFD,
                             *JFD,
                             dt,
                             relativeConvergenceTolerance,
                             maxNumPoissonIterations);

        bc.resetOutsideDomain(*EFD, geo);
    } else {
        bunch.makeReady(dt/2);
    }
    NDIndex<DIM> oldGlobalPos = bunch.getGlobalPDomainInclOld(miniMesh);

    {
        std::vector<NDIndex<DIM> > lSDoms, lEDoms;
        Utils::getLocalDomains(EFD->getLayout(), lSDoms);
        Utils::getLocalDomains(miniMeshEdgeFL, lEDoms);
        unsigned int totalSSize = EFD->getLayout().getDomain().size();
        unsigned int totalESize = miniMeshEdgeFL.getDomain().size();
        double idealSSize = totalSSize * (1.0 / Ippl::getNodes());
        double idealESize = totalESize * (1.0 / Ippl::getNodes());
        if (Ippl::myNode() == 0) {
            msg << " ----------------------------------------------------\n"
                << "    L O C A L  D O M A I N E S ... \n"
                << " ----------------------------------------------------\n";
            msg << " SURROUNDING MESH \n";
            for (int i = 0; i < Ippl::getNodes(); ++ i) {
                int size = lSDoms[i].size();
                double diff = (size - idealSSize) / idealSSize;
                msg << "Node" << std::setw(5) << i << ": {"
                    << "[" << std::setw(5) << lSDoms[i][0].first() << ":" << std::setw(5) << lSDoms[i][0].last() << "],"
                    << "[" << std::setw(5) << lSDoms[i][1].first() << ":" << std::setw(5) << lSDoms[i][1].last() << "]"
                    << "}, size: " << std::setw(7) << size
                    << std::setw(8) << std::setiosflags(std::ios::showpos)
                    << std::setprecision(3) << diff * 100 << " %   "
                    << std::resetiosflags(std::ios::showpos)
                    << "\n";
            }
            msg << "\n EMBEDDED MESH \n";
            msg << "shift: " << ddim[0].first() << ", " << ddim[1].first() << "\n";
            for (int i = 0; i < Ippl::getNodes(); ++ i) {
                int size = lEDoms[i].size();
                double diff = (size - idealESize) / idealESize;
                msg << "Node" << std::setw(5) << i << ": {"
                    << "[" << std::setw(5) << lEDoms[i][0].first() << ":" << std::setw(5) << lEDoms[i][0].last() << "],"
                    << "[" << std::setw(5) << lEDoms[i][1].first() << ":" << std::setw(5) << lEDoms[i][1].last() << "]"
                    << "}, size: " << std::setw(7) << size
                    << std::setw(8) << std::setiosflags(std::ios::showpos)
                    << std::setprecision(3) << diff * 100 << " %   "
                    << std::resetiosflags(std::ios::showpos)
                    << "\n";
            }
            msg << " ----------------------------------------------------\n"
                << " ----------------------------------------------------\n"
                << endl;
        }

        const std::vector<size_t> & localNums = bunch.getNumLocalParticles();
        // const std::vector<unsigned int>  localNums(Ippl::getNodes(), 0);
        // unsigned int myLocalNum = bunch.getLocalNum();
        // boost::mpi::gather(world, myLocalNum, localNums, 0);
        unsigned int totalNum = bunch.getTotalNum();
        double idealNum = totalNum * (1.0 / Ippl::getNodes());
        msg << " ----------------------------------------------------\n"
            << "    L O C A L  P A R T I C L E S ... \n"
            << " ----------------------------------------------------\n";
        for (unsigned int k = 0; k < localNums.size(); ++ k) {
            double diff = (localNums[k] - idealNum) / idealNum;
            msg << "Node "
                << std::setw(5) << k << ": "
                << std::setw(8) << localNums[k]
                << std::setw(8) << std::setiosflags(std::ios::showpos)
                << std::setprecision(3) << diff * 100 << " %   "
                << std::resetiosflags(std::ios::showpos)
                << "\n";
        }
        msg << endl;

        NDIndex<DIM> lPDom = bunch.getLocalPDomain(miniMesh);
        const std::vector<NDIndex<DIM> > localPDomains = Communicator::getAllLocalDomains(lPDom);
        msg << " ----------------------------------------------------\n"
            << "    L O C A L  P A R T I C L E  D O M A I N S ... \n"
            << " ----------------------------------------------------\n";
        for (unsigned int k = 0; k < localNums.size(); ++ k) {
            msg << "Node "
                << std::setw(5) << k << ": "
                << "[" << std::setw(5) << localPDomains[k][0].first() << ":" << std::setw(5) << localPDomains[k][0].last() << "],"
                << "[" << std::setw(5) << localPDomains[k][1].first() << ":" << std::setw(5) << localPDomains[k][1].last() << "]"
                << "}" << "\n";
        }
        msg << endl;


    }

    copyFields(miniMeshE, miniMeshH, miniMeshJ,
               *EFD, *HFD, *JFD,
               Vektor<int, DIM>(ddim[0].min(), ddim[1].min()), true, false);

    Timings::startTimer(othersTimer);
    exchangeBoundaries(*EFD, *HFD,
                       miniMeshE, miniMeshH,
                       Vektor<int, DIM>(ddim[0].min(), ddim[1].min()));

    bunch.balanceParticles();

    miniMeshE.fillGuardCells();
    miniMeshH.fillGuardCells();

    TETMSolver fieldSolver(*AlphaE,
                           *AlphaH,
                           bunch,
                           msg,
                           *EFD,
                           dt);

    // CalcEnergyCommand cec(*EFD);

    PlotBinaryVTKCommand pavc(*EFD, *HFD, *JFD);
    PlotBinaryVTKCommand mmpavc(miniMeshE, miniMeshH, miniMeshJ, "mm");

    // SaveH5Part h5file(cmdl, effectiveFileName, saveFlag);
    if (restart.getValue() == "") {
    //     // output initial state
        pavc.execute(0.0, bunch, msg);
        mmpavc.execute(0.0, bunch, msg);
    //     h5file.write(*EFD, *HFD, *JFD, bunch, currentSimulatedTime, msg);
    }

    SaveEMFieldsCommand semfc;
    IntegrateElongCommand ielc(miniMeshE, statModulus.getValue() * dt);
    DataSinkCommand dsc(*EFD);
    PartBunchState bunchState;
    dsc.execute(bunch,
                currentSimulatedTime + dt,
                *EFD,
                *HFD,
                miniMeshE,
                miniMeshH,
                Vektor<int,DIM>(ddim[0].min(), ddim[1].min()),
                bunchState,
                msg);
    // dsc.execute(bunch, currentSimulatedTime, bunchState, msg);

    semfc.execute(*EFD, *HFD, *JFD, bunchState, msg);

    Timings::stopTimer(othersTimer);

    (*JFD) = 0.0;
    // begin main loop /////////////////////////////////////////////////////////////////////////
    for(int iteration = startAtStep; iteration < startAtStep + maximumIterations; iteration++) {
        msg << "Iteration " << iteration
            << ", current simulation time: " <<  currentSimulatedTime << endl;
        DBGOUT << "new timestep; iteration =  " << iteration << std::endl;

        fieldSolver.updateField_firstStep(*EFD,
                                          *HFD,
                                          *JFD,
                                          miniMeshE,
                                          miniMeshH,
                                          miniMeshJ,
                                          ddim,
                                          bc,
                                          geo);

        fieldSolver.updateField_secondStep(*EFD,
                                           *HFD,
                                           *JFD,
                                           miniMeshE,
                                           miniMeshH,
                                           miniMeshJ,
                                           ddim,
                                           bc,
                                           geo);

        // begin visual output /////////////////////////////////////////////////
        if (iteration % plotModulus.getValue() == plotModulus.getValue() - 1) {
            copyFields(*EFD, *HFD, *JFD,
                       miniMeshE, miniMeshH, miniMeshJ,
                       Vektor<int, DIM>(-ddim[0].min(), -ddim[1].min()), true, true);

            pavc.execute(currentSimulatedTime + dt, bunch, msg);
            mmpavc.execute(currentSimulatedTime + dt, bunch, msg);
        }

        // if (iteration % saveModulus.getValue() == saveModulus.getValue() - 1) {
        //     copyFields(*EFD, *HFD, *JFD,
        //                miniMeshE, miniMeshH, miniMeshJ,
        //                Vektor<int, DIM>(-ddim[0].min(), -ddim[1].min()), true, true);
        //     h5file.write(*EFD, *HFD, *JFD, bunch, currentSimulatedTime + dt, msg);
        // }

        // end visual output ///////////////////////////////////////////////////

        // begin statistic output //////////////////////////////////////////////
        if (iteration % statModulus.getValue() == statModulus.getValue() - 1) {
            dsc.execute(bunch,
                        currentSimulatedTime + dt,
                        *EFD,
                        *HFD,
                        miniMeshE,
                        miniMeshH,
                        Vektor<int,DIM>(ddim[0].min(), ddim[1].min()),
                        bunchState,
                        msg);
            if (inStraightSectionBeforeBend && geo.isInsideBend(bunchState.get_rmax())) {

                ielc.write("Data/energyChangeBefore.dat",
                           miniMeshE,
                           bunchState);
                ielc.reset();
                inStraightSectionBeforeBend = false;
                inBend = true;
            }
            if (inBend &&
                !geo.isInsideBend(bunchState.get_rmax()) &&
                !geo.isInsideBend(bunchState.get_rmin())) {

                ielc.write("Data/energyChangeInBend.dat",
                           miniMeshE,
                           bunchState);
                ielc.reset();
                inBend = false;
            }

            ielc.execute(miniMeshE, bunchState, msg);
        }

        if (iteration % lineModulus.getValue() == lineModulus.getValue() - 1) {
            if (iteration % statModulus.getValue() != statModulus.getValue() - 1) {
                bunch.calcBeamParameters(bunchState);
            }
            if (iteration % plotModulus.getValue() != plotModulus.getValue() - 1) {
                copyFields(*EFD, *HFD, *JFD,
                           miniMeshE, miniMeshH, miniMeshJ,
                           Vektor<int, DIM>(-ddim[0].min(), -ddim[1].min()), true, false);
            }
            semfc.execute(*EFD, *HFD, *JFD, bunchState, msg);
        }

        // cec.execute(*EFD,
        //             *HFD,
        //             miniMeshE,
        //             miniMeshH,
        //             Vektor<int,DIM>(ddim[0].min(), ddim[1].min()),
        //             currentSimulatedTime + dt,
        //             msg);

        // end statistic output ////////////////////////////////////////////////

        // if (iteration % rebModulus.getValue() == rebModulus.getValue() - 1) {
        //     bunch.rebalance();
        // }

        (*JFD) = 0.0;

        Timings::startTimer(miniMeshTimer);
        moveMesh(miniMeshE, miniMeshH, miniMeshJ,
                 *EFD, *HFD,
                 ddim,
                 bunch.getGlobalBounds());
        Timings::stopTimer(miniMeshTimer);

        // time in simulated seconds that the simulation has progressed:
        currentSimulatedTime += dt;

    }
    // end main loop ///////////////////////////////////////////////////////////////////////////

    msg << "\n"
        << "\n"
        << "Meshing parameters: \n"
        << "cells " <<  Nx << " " <<  Ny << "\n"
        << "mesh spacing " << hr(0) << " " << hr(1)
        << endl;

    bunch.calcBeamParameters(bunchState);
    ielc.write("Data/energyChangeAfter.dat",
               miniMeshE,
               bunchState);

    pavc.closeFile();
    mmpavc.closeFile();

    // free fields
    delete lhpc;

    if (dbg_output) {
        dbg.close();
    }
    Timings::stopTimer(ipplToyFDTDTimer);
    IpplTimings::print();
}

void copyFields(VField_Edge_t & toE,
                VField_Cell_t & toH,
                VField_Edge_t & toJ,
                const VField_Edge_t & fromE,
                const VField_Cell_t & fromH,
                const VField_Edge_t & fromJ,
                const Vektor<int, DIM> & start,
                const bool & all2all,
                const bool & copyCurrent)
{
    Communicator::copy(fromE, toE, start, all2all);
    Communicator::copy(fromH, toH, start, all2all);
    if (copyCurrent)
        Communicator::copy(fromJ, toJ, start, all2all);
}

void exchangeBoundaries(VField_Edge_t & surroundingE,
                        VField_Cell_t & surroundingH,
                        VField_Edge_t & embeddedE,
                        VField_Cell_t & embeddedH,
                        const Vektor<int, DIM> & start)
{
    Communicator::exchangeBoundaries(surroundingH, embeddedH, start, HFieldToDouble());
    Communicator::exchangeBoundaries(surroundingE, embeddedE, start, ExFieldToDouble());
    Communicator::exchangeBoundaries(surroundingE, embeddedE, start, EyFieldToDouble());
}


void moveMesh(VField_Edge_t & miniMeshE,
              VField_Cell_t & miniMeshH,
              VField_Edge_t & miniMeshJ,
              VField_Edge_t & EFD,
              VField_Cell_t & HFD,
              NDIndex<DIM> & ddim,
              const NDIndex<DIM> & newGlobalPos)
{
    NDIndex<DIM> idealPos;
    bool fits = true;
    for (unsigned int d = 0; d < DIM; ++ d) {
        if (newGlobalPos[d].length() >= ddim[d].length()) {
            fits = false;
            break;
        } else {
            unsigned int diff = ddim[d].length() - newGlobalPos[d].length();
            unsigned int upperMargin =
                static_cast<unsigned int>(std::floor(diff / 2.0));
            unsigned int lowerMargin = diff - upperMargin;
            idealPos[d] = Index(lowerMargin - 1, ddim[d].length() - upperMargin - 2);
        }
    }
    if (!fits) {
        DBGOUT << "YOU HAVE A MAJOR PROBLEM" << std::endl;
        PAssert(fits);
    }

    int moveI = newGlobalPos[0].last() - idealPos[0].last() - 1;
    int maxJ = std::max(0, newGlobalPos[1].last() - idealPos[1].last() - 1);
    int minJ = std::min(0, newGlobalPos[1].last() - idealPos[1].last() - 1);

    NDIndex<DIM> ldomMM = miniMeshE.getLayout().getLocalNDIndex();
    Index immII(ldomMM[0].first() - moveI, ldomMM[0].last());
    Index immJJ(ldomMM[1].first() - maxJ, ldomMM[1].last() - minJ);

    if (minJ + maxJ == 0) {
        immJJ = Index(immJJ.first() - 1, immJJ.last() + 1);
    }

    if (minJ == -1) {
        NDIndex<DIM> elem1, elem2;
        for (int j = immJJ.last(); j >= immJJ.first(); -- j) {
            elem1[1] = Index(j, j);
            elem2[1] = Index(j + minJ + maxJ, j + minJ + maxJ);
            for (int i = immII.first(); i <= immII.last(); ++ i) {
                elem1[0] = Index(i,i);
                elem2[0] = Index(i + moveI, i + moveI);
                miniMeshE.localElement(elem1) = miniMeshE.localElement(elem2);
                miniMeshH.localElement(elem1) = miniMeshH.localElement(elem2);
                miniMeshJ.localElement(elem1) = miniMeshJ.localElement(elem2);
            }
        }
    } else {
        NDIndex<DIM> elem1, elem2;
        for (int j = immJJ.first(); j <= immJJ.last(); ++ j) {
            elem1[1] = Index(j, j);
            elem2[1] = Index(j + minJ + maxJ, j + minJ + maxJ);
            for (int i = immII.first(); i <= immII.last(); ++ i) {
                elem1[0] = Index(i,i);
                elem2[0] = Index(i + moveI, i + moveI);
                miniMeshE.localElement(elem1) = miniMeshE.localElement(elem2);
                miniMeshH.localElement(elem1) = miniMeshH.localElement(elem2);
                miniMeshJ.localElement(elem1) = miniMeshJ.localElement(elem2);
            }
        }
    }

    ddim[0] = ddim[0] + moveI;
    ddim[1] = ddim[1] + minJ + maxJ;

    Vector_t miniMeshOrigin = miniMeshE.get_mesh().get_origin();
    Vector_t hr(miniMeshE.get_mesh().get_meshSpacing(0),
                miniMeshE.get_mesh().get_meshSpacing(1));
    miniMeshOrigin(0) += moveI * hr(0);
    miniMeshOrigin(1) += (minJ + maxJ) * hr(1);

    miniMeshE.get_mesh().set_origin(miniMeshOrigin);
    miniMeshH.get_mesh().set_origin(miniMeshOrigin);

    Communicator::exchangeBoundariesMovingMesh(EFD,
                                               HFD,
                                               miniMeshE,
                                               miniMeshH,
                                               Vektor<int,DIM>(ddim[0].first(), ddim[1].first()),
                                               Vektor<int,DIM>(moveI, minJ + maxJ));

    miniMeshE.fillGuardCells();
    miniMeshH.fillGuardCells();
    EFD.fillGuardCells();
    HFD.fillGuardCells();
}
