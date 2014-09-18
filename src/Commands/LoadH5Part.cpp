/***************************************************************************
                            LoadH5Part.cpp
                         -------------------
    begin                : Mon Oct 4 2010
    copyright            : (C) 2010 by Christof Kraus
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

#include "LoadH5Part.hh"
#include "../Physics.hh"
#include "../PartBunch.hh"

#include <set>

extern std::ostream dbg;

LoadH5Part::LoadH5Part(const std::string & fname,
                       TCLAP::CmdLine & cmdl):
    _H5file(0),
    _fileName(fname)
{
    openFile();
    readH5FileAttributes(cmdl);
    closeFile();
    _inputTimer = Timings::getTimer("h5block_input");
}

LoadH5Part::~LoadH5Part()
{ }

void LoadH5Part::closeFile()
{
    H5CloseFile(_H5file);
    Ippl::Comm->barrier();

    _H5file = 0;
}

void LoadH5Part::openFile()
{
    _H5file = H5OpenFile(_fileName.c_str(), H5_O_RDONLY, MPI_COMM_WORLD);
    if ((char*)_H5file == (char*)H5_ERR) {
        std::cerr << "h5 file open failed! exiting" << std::endl;
        exit(0);
    }

    h5_int64_t lastStep = H5GetNumSteps(_H5file) - 1;
    H5SetStep(_H5file, lastStep);
}

void LoadH5Part::readH5FileAttributes(TCLAP::CmdLine & cmdl)
{
    std::set<std::string> excludedArgs({"numgridpx", "numgridpy", "Np", "phi"});
    std::list<TCLAP::Arg*> callArguments = cmdl.getArgList();
    for (auto arg: callArguments) {
        if (arg->isSet() &&
            excludedArgs.find(arg->getName()) == excludedArgs.end()) {
            continue;
        }

        if (arg->getName() == "startstep") {
            TCLAP::ValueArg<int> & Varg = *dynamic_cast<TCLAP::ValueArg<int>*>(arg);
            h5_int64_t value;
            std::string name = Varg.getName();
            H5ReadFileAttribInt64(_H5file, name.c_str(), &value);
            int startStep = value;

            name = "numsteps";
            H5ReadFileAttribInt64(_H5file, name.c_str(), &value);
            int numSteps = value;

            name = "savemodulus";
            H5ReadFileAttribInt64(_H5file, name.c_str(), &value);
            int saveModulus = value;

            Varg.getValue() = ((startStep + numSteps) / saveModulus) * saveModulus;

            continue;
        }

        if (dynamic_cast<TCLAP::ValueArg<std::string>*>(arg)) {
            TCLAP::ValueArg<std::string> & Varg = *dynamic_cast<TCLAP::ValueArg<std::string>*>(arg);
            char value[512];
            std::string name = Varg.getName();
            H5ReadFileAttribString(_H5file, name.c_str(), value);
            Varg.getValue() = std::string(value);

            continue;
        }

        if (dynamic_cast<TCLAP::ValueArg<double>*>(arg)) {
            TCLAP::ValueArg<double> & Varg = *dynamic_cast<TCLAP::ValueArg<double>*>(arg);
            h5_float64_t value;
            std::string name = Varg.getName();
            H5ReadFileAttribFloat64(_H5file, name.c_str(), &value);
            Varg.getValue() = value;

            continue;
        }

        if (dynamic_cast<TCLAP::ValueArg<int>*>(arg)) {
            TCLAP::ValueArg<int> & Varg = *dynamic_cast<TCLAP::ValueArg<int>*>(arg);
            h5_int64_t value;
            std::string name = Varg.getName();
            H5ReadFileAttribInt64(_H5file, name.c_str(), &value);
            Varg.getValue() = value;

            continue;
        }

        if (dynamic_cast<TCLAP::ValueArg<unsigned>*>(arg)) {
            TCLAP::ValueArg<unsigned> & Varg = *dynamic_cast<TCLAP::ValueArg<unsigned>*>(arg);
            h5_int64_t value;
            std::string name = Varg.getName();
            H5ReadFileAttribInt64(_H5file, name.c_str(), &value);
            Varg.getValue() = value;

            continue;
        }

        if (dynamic_cast<TCLAP::ValueArg<Vect2D>*>(arg)) {
            TCLAP::ValueArg<Vect2D> & Varg = *dynamic_cast<TCLAP::ValueArg<Vect2D>*>(arg);
            h5_float64_t value_f64[2];
            std::string name = Varg.getName();
            H5ReadFileAttribFloat64(_H5file, name.c_str(), value_f64);
            Vect2D & value = Varg.getValue();
            value.v[0] = value_f64[0];
            value.v[1] = value_f64[1];

            continue;
        }

        if (dynamic_cast<TCLAP::SwitchArg*>(arg)) {
            TCLAP::SwitchArg & Sarg = *dynamic_cast<TCLAP::SwitchArg*>(arg);
            h5_int64_t value_i64;
            std::string name = Sarg.getName();
            H5ReadFileAttribInt64(_H5file, name.c_str(), &value_i64);
            Sarg.setValue(value_i64 == 1);

            continue;
        }

        std::cerr << "**** can't save the argument '" << arg->getName()
                  << "' because it's type can't be determined ****\n"
                  << "**** going to exit ****" << std::endl;

        exit(0);

    }
}

void LoadH5Part::restoreFields(VField_Edge_t & EFD,
                               VField_Cell_t & HFD,
                               VField_Edge_t & JFD,
                               double & t)
{
    FieldLayout<DIM> & FL = EFD.getLayout();
    const NDIndex<DIM> & lDom = FL.getLocalNDIndex();
    NDIndex<DIM> elem;

    const int Nx = lDom[0].length();
    const int Ny = lDom[1].length();

    h5_float64_t* x_data = new h5_float64_t[Nx * Ny];
    h5_float64_t* y_data = new h5_float64_t[Nx * Ny];
    h5_float64_t* z_data = new h5_float64_t[Nx * Ny];

    openFile();
    h5_ssize_t lastStep = H5GetNumSteps(_H5file) - 1;
    H5SetStep(_H5file, lastStep);
    H5ReadStepAttribFloat64(_H5file, "TIME", &t);
    t *= 1e-9;
    H5Block3dSetView(_H5file,
                     lDom[0].first(), lDom[0].last(),
                     lDom[1].first(), lDom[1].last(),
                     0, 0);

    long ii = 0;
    H5Block3dReadVector3dFieldFloat64(_H5file,  "Efield", x_data, y_data, z_data);
    for (int j = lDom[1].first(); j <= lDom[1].last(); ++ j) {
        elem[1] = Index(j, j);
        for (int i = lDom[0].first(); i <= lDom[0].last(); ++ i) {
            elem[0] = Index(i, i);
            EFD.localElement(elem)(0) = x_data[ii];
            EFD.localElement(elem)(1) = y_data[ii];
            ++ii;
        }
    }

    ii = 0;
    H5Block3dReadVector3dFieldFloat64(_H5file,  "Hfield", x_data, y_data, z_data);
    for (int j = lDom[1].first(); j <= lDom[1].last(); ++ j) {
        elem[1] = Index(j, j);
        for (int i = lDom[0].first(); i <= lDom[0].last(); ++ i) {
            elem[0] = Index(i, i);
            HFD.localElement(elem)(0) = x_data[ii];
            HFD.localElement(elem)(1) = y_data[ii];
            ++ii;
        }
    }

    ii = 0;
    H5Block3dReadVector3dFieldFloat64(_H5file,  "Jfield", x_data, y_data, z_data);
    for (int j = lDom[1].first(); j <= lDom[1].last(); ++ j) {
        elem[1] = Index(j, j);
        for (int i = lDom[0].first(); i <= lDom[0].last(); ++ i) {
            elem[0] = Index(i, i);
            JFD.localElement(elem)(0) = x_data[ii];
            JFD.localElement(elem)(1) = y_data[ii];
            ++ii;
        }
    }

    h5_float64_t dx;
    h5_float64_t dy;
    h5_float64_t dz;
    H5Block3dGetFieldSpacing(_H5file, "Efield", &dx, &dy, &dz);
    Vector_t hr(dx, dy);
    EFD.get_mesh().set_meshSpacing(&(hr[0]));
    HFD.get_mesh().set_meshSpacing(&(hr[0]));

    H5Block3dGetFieldOrigin(_H5file, "Efield", &dx, &dy, &dz);
    Vector_t origin(dx, dy);
    EFD.get_mesh().set_origin(origin);

    closeFile();

    delete[] x_data;
    delete[] y_data;
    delete[] z_data;
}

void LoadH5Part::restoreBunch(std::vector<Vector_t> & R,
                              std::vector<Vector_t> & oldR,
                              std::vector<Vector_t> & P,
                              std::vector<double> & Q)
{
    openFile();

    h5_int64_t StartView;
    h5_int64_t EndView;

    h5_int64_t nPartTot = H5PartGetNumParticles(_H5file);
    h5_int64_t numLocalParticles = std::floor(static_cast<double>(nPartTot) / Ippl::getNodes());
    h5_int64_t missing = nPartTot - numLocalParticles * Ippl::getNodes();
    if (Ippl::myNode() < missing) {
        numLocalParticles ++;
        StartView = Ippl::myNode() * numLocalParticles;
    } else {
        StartView = Ippl::myNode() * numLocalParticles + missing;
    }
    EndView = StartView + numLocalParticles - 1;

    R.resize(numLocalParticles);
    oldR.resize(numLocalParticles);
    P.resize(numLocalParticles);
    Q.resize(numLocalParticles);

    // StartView = 0;
    // EndView = nPartTot - 1;
    // numLocalParticles = EndView - StartView + 1;

    H5PartSetView(_H5file, StartView, EndView);

    h5_float64_t * part_data = new h5_float64_t[2 * numLocalParticles];

    h5_float64_t * x_data = part_data;
    h5_float64_t * y_data = part_data + numLocalParticles;

    H5PartReadDataFloat64(_H5file, "x", x_data);
    H5PartReadDataFloat64(_H5file, "y", y_data);

    size_t ii = 0;
    for (int i = 0; i < numLocalParticles; ++ i) {
        // if (i % Ippl::getNodes() == Ippl::myNode()) {
        R[ii++] = Vector_t(x_data[i], y_data[i]);
        // }
    }

    H5PartReadDataFloat64(_H5file, "oldx", x_data);
    H5PartReadDataFloat64(_H5file, "oldy", y_data);

    ii = 0;
    for (int i = 0; i < numLocalParticles; ++ i) {
        // if (i % Ippl::getNodes() == Ippl::myNode()) {
        oldR[ii++] = Vector_t(x_data[i], y_data[i]);
        // }
    }

    H5PartReadDataFloat64(_H5file, "px", x_data);
    H5PartReadDataFloat64(_H5file, "py", y_data);

    ii = 0;
    for (int i = 0; i < numLocalParticles; ++ i) {
        // if (i % Ippl::getNodes() == Ippl::myNode()) {
        P[ii++] = Vector_t(x_data[i], y_data[i]);
        // }
    }

    H5PartReadDataFloat64(_H5file, "q", x_data);

    ii = 0;
    for (int i = 0; i < numLocalParticles; ++ i) {
        // if (i % Ippl::getNodes() == Ippl::myNode()) {
        Q[ii++] = x_data[i];
        // }
    }

    delete[] part_data;

    closeFile();
}
