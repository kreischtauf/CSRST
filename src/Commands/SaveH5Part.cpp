/***************************************************************************
                            SaveH5Part.cpp
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

#include "SaveH5Part.hh"
#include "../PartBunch.hh"

extern std::ostream dbg;

SaveH5Part::SaveH5Part(TCLAP::CmdLine & cmdl,
                       const std::string & filename,
                       const h5_int32_t & flag):
    _x_data(0),
    _y_data(0),
    _z_data(0),
    _part_data(0),
    _H5call(0),
    _prev_nloc(0),
    _domain_size(0),
    _fileName(filename)
{
    openFile(flag);
    if (flag != H5_O_APPEND) {
        writeH5FileAttributes(cmdl);
    } else {
        writeReducedH5FileAttributes(cmdl);
        _H5call = H5GetNumSteps(_H5file);
    }
    closeFile();

    _outputTimer = Timings::getTimer("h5block_output");
}

SaveH5Part::~SaveH5Part()
{
    if (_H5call > 0) {
        delete[] _part_data;
        delete[] _x_data;
        delete[] _y_data;
        delete[] _z_data;
    }

    _part_data = 0;
    _x_data = 0;
    _y_data = 0;
    _z_data = 0;
}

void SaveH5Part::closeFile()
{
    H5CloseFile(_H5file);
    Ippl::Comm->barrier();

    _H5file = 0;
}

void SaveH5Part::openFile(const h5_int32_t & flag)
{
    _H5file = H5OpenFile(_fileName.c_str(), flag, MPI_COMM_WORLD);
    if ((char*)_H5file == (char*)H5_ERR) {
        std::cerr << "h5 file open failed! exiting" << std::endl;
        exit(0);
    }
}

void SaveH5Part::write(const VField_Edge_t & EFD,
                       const VField_Cell_t & HFD,
                       const VField_Edge_t & JFD,
                       const PartBunch & bunch,
                       const double & t,
                       Inform & msg)
{
    FieldLayout<DIM> & FL = EFD.getLayout();
    const NDIndex<DIM> & lDom = FL.getLocalNDIndex();
    NDIndex<DIM> elem;

    h5_int64_t *int_part_data;
    const size_t effLocalNumParticles = bunch.getLocalNP();
    const size_t one{1};
    h5_int64_t nLoc = std::max(one, effLocalNumParticles);

    const double time = t * 1e9; // save time in ns

    const int Nx = lDom[0].length();
    const int Ny = lDom[1].length();

    if (Nx * Ny > _domain_size) {
        delete[] _x_data;
        delete[] _y_data;
        delete[] _z_data;

        _x_data = new h5_float64_t[Nx * Ny];
        _y_data = new h5_float64_t[Nx * Ny];
        _z_data = new h5_float64_t[Nx * Ny];
        _domain_size = Nx * Ny;
    }

    if (nLoc > _prev_nloc) {
        delete[] _part_data;
        _part_data = new h5_float64_t[nLoc];
        _prev_nloc = nLoc;
    }
    int_part_data = (h5_int64_t *)_part_data;
    _part_data[0] = 0.0;

    Timings::startTimer(_outputTimer);

    openFile(H5_O_APPEND);
    H5SetStep(_H5file, _H5call);

    H5WriteStepAttribFloat64(_H5file, "TIME", &time, 1);

    H5PartSetNumParticles(_H5file, nLoc);

    for (size_t i = 0; i < effLocalNumParticles; ++ i) {
        _part_data[i] = bunch.R[i](0);
    }
    H5PartWriteDataFloat64(_H5file, "x", _part_data);

    for (size_t i = 0; i < effLocalNumParticles; ++ i) {
        _part_data[i] = bunch.R[i](1);
    }
    H5PartWriteDataFloat64(_H5file, "y", _part_data);

    for (size_t i = 0; i < effLocalNumParticles; ++ i) {
        _part_data[i] = 1.00000001;
    }
    H5PartWriteDataFloat64(_H5file, "z", _part_data);

    for (size_t i = 0; i < effLocalNumParticles; ++ i) {
        _part_data[i] = bunch.oldR[i](0);
    }
    H5PartWriteDataFloat64(_H5file, "oldx", _part_data);

    for (size_t i = 0; i < effLocalNumParticles; ++ i) {
        _part_data[i] = bunch.oldR[i](1);
    }
    H5PartWriteDataFloat64(_H5file, "oldy", _part_data);

    for (size_t i = 0; i < effLocalNumParticles; ++ i) {
        _part_data[i] = bunch.P[i](0);
    }
    H5PartWriteDataFloat64(_H5file, "px", _part_data);

    for (size_t i = 0; i < effLocalNumParticles; ++ i) {
        _part_data[i] = bunch.P[i](1);
    }
    H5PartWriteDataFloat64(_H5file, "py", _part_data);

    for (size_t i = 0; i < effLocalNumParticles; ++ i) {
        _part_data[i] = 0.0;
    }
    H5PartWriteDataFloat64(_H5file, "pz", _part_data);

    for (size_t i = 0; i < effLocalNumParticles; ++ i) {
        _part_data[i] = bunch.Q[i];
    }
    H5PartWriteDataFloat64(_H5file, "q", _part_data);

    int_part_data[0] = -1;
    for(size_t i = 0; i < effLocalNumParticles; ++ i) {
        int_part_data[i] =  bunch.ID[i];
    }
    H5PartWriteDataInt64(_H5file, "id", int_part_data);


    H5Block3dSetView(_H5file,
                     lDom[0].first(), lDom[0].last(),
                     lDom[1].first(), lDom[1].last(),
                     0, 0);

    long ii = 0;
    // h5block uses the fortran convention of storing data:
    // INTEGER, DIMENSION(2,3) :: a
    // => {a(1,1), a(2,1), a(1,2), a(2,2), a(1,3), a(2,3)}
    for (int j = lDom[1].first(); j <= lDom[1].last(); ++ j) {
        elem[1] = Index(j, j);
        for (int i = lDom[0].first(); i <= lDom[0].last(); ++ i) {
            elem[0] = Index(i, i);
            _x_data[ii] = EFD.localElement(elem)(0);
            _y_data[ii] = EFD.localElement(elem)(1);
            _z_data[ii] = 0.0;
            ++ii;
        }
    }
    H5Block3dWriteVector3dFieldFloat64(_H5file,  "Efield", _x_data, _y_data, _z_data);

    h5_float64_t dx = EFD.get_mesh().get_meshSpacing(0);
    h5_float64_t dy = EFD.get_mesh().get_meshSpacing(1);
    H5Block3dSetFieldSpacing(_H5file,
                             "Efield",
                             dx,
                             dy,
                             h5_float64_t(1.0));

    Vector_t origin = EFD.get_mesh().get_origin();
    H5Block3dSetFieldOrigin(_H5file,
                            "Efield",
                            h5_float64_t(origin(0)),
                            h5_float64_t(origin(1)),
                            h5_float64_t(0.0));

    ii = 0;
    for (int j = lDom[1].first(); j <= lDom[1].last(); ++ j) {
        elem[1] = Index(j, j);
        for (int i = lDom[0].first(); i <= lDom[0].last(); ++ i) {
            elem[0] = Index(i, i);
            _x_data[ii] = HFD.localElement(elem)(0);
            _y_data[ii] = HFD.localElement(elem)(1);
            _z_data[ii] = HFD.localElement(elem)(0) + HFD.localElement(elem)(1);
            ++ii;
        }
    }
    H5Block3dWriteVector3dFieldFloat64(_H5file,  "Hfield", _x_data, _y_data, _z_data);

    H5Block3dSetFieldSpacing(_H5file,
                             "Hfield",
                             dx,
                             dy,
                             h5_float64_t(1.0));

    H5Block3dSetFieldOrigin(_H5file,
                            "Hfield",
                            h5_float64_t(origin(0) + 0.5 * dx),
                            h5_float64_t(origin(1) + 0.5 * dy),
                            h5_float64_t(0.0));

    ii = 0;
    for (int j = lDom[1].first(); j <= lDom[1].last(); ++ j) {
        elem[1] = Index(j, j);
        for (int i = lDom[0].first(); i <= lDom[0].last(); ++ i) {
            elem[0] = Index(i, i);
            _x_data[ii] = JFD.localElement(elem)(0);
            _y_data[ii] = JFD.localElement(elem)(1);
            _z_data[ii] = 0.0;
            ++ii;
        }
    }

    H5Block3dWriteVector3dFieldFloat64(_H5file,  "Jfield", _x_data, _y_data, _z_data);
    H5Block3dSetFieldSpacing(_H5file,
                             "Jfield",
                             dx,
                             dy,
                             h5_float64_t(1.0));
    H5Block3dSetFieldOrigin(_H5file,
                            "Jfield",
                            h5_float64_t(origin(0)),
                            h5_float64_t(origin(1)),
                            h5_float64_t(0.0));

    Timings::stopTimer(_outputTimer);

    ++_H5call;
    closeFile();
}


void SaveH5Part::writeH5FileAttributes(TCLAP::CmdLine & cmdl)
{
    H5WriteFileAttribString(_H5file, "tUnit", "s");
    H5WriteFileAttribString(_H5file, "xUnit", "m");
    H5WriteFileAttribString(_H5file, "yUnit", "m");
    H5WriteFileAttribString(_H5file, "zUnit", "m");
    H5WriteFileAttribString(_H5file, "pxUnit", "#beta#gamma");
    H5WriteFileAttribString(_H5file, "pyUnit", "#beta#gamma");
    H5WriteFileAttribString(_H5file, "pzUnit", "#beta#gamma");
    H5WriteFileAttribString(_H5file, "qUnit", "Cb");

    H5WriteFileAttribString(_H5file, "idUnit", "1");

    H5WriteFileAttribString(_H5file, "SPOSUnit", "m");
    H5WriteFileAttribString(_H5file, "TIMEUnit", "ns");

    std::list<TCLAP::Arg*> callArguments = cmdl.getArgList();
    for (auto arg: callArguments) {
        if (dynamic_cast<TCLAP::ValueArg<std::string>*>(arg)) {
            TCLAP::ValueArg<std::string> & Varg = *dynamic_cast<TCLAP::ValueArg<std::string>*>(arg);
            std::string value = Varg.getValue();
            std::string name = Varg.getName();
            H5WriteFileAttribString(_H5file, name.c_str(), value.c_str());

            continue;
        }

        if (dynamic_cast<TCLAP::ValueArg<double>*>(arg)) {
            TCLAP::ValueArg<double> & Varg = *dynamic_cast<TCLAP::ValueArg<double>*>(arg);
            h5_float64_t value = Varg.getValue();
            std::string name = Varg.getName();
            H5WriteFileAttribFloat64(_H5file, name.c_str(), &value, 1);

            continue;
        }

        if (dynamic_cast<TCLAP::ValueArg<int>*>(arg)) {
            TCLAP::ValueArg<int> & Varg = *dynamic_cast<TCLAP::ValueArg<int>*>(arg);
            h5_int64_t value = Varg.getValue();
            std::string name = Varg.getName();
            H5WriteFileAttribInt64(_H5file, name.c_str(), &value, 1);

            continue;
        }

        if (dynamic_cast<TCLAP::ValueArg<unsigned>*>(arg)) {
            TCLAP::ValueArg<unsigned> & Varg = *dynamic_cast<TCLAP::ValueArg<unsigned>*>(arg);
            h5_int64_t value = Varg.getValue();
            std::string name = Varg.getName();
            H5WriteFileAttribInt64(_H5file, name.c_str(), &value, 1);

            continue;
        }

        if (dynamic_cast<TCLAP::ValueArg<Vect2D>*>(arg)) {
            TCLAP::ValueArg<Vect2D> & Varg = *dynamic_cast<TCLAP::ValueArg<Vect2D>*>(arg);
            Vect2D value = Varg.getValue();
            h5_float64_t value_f64[] = {value.v[0], value.v[1]};
            std::string name = Varg.getName();
            H5WriteFileAttribFloat64(_H5file, name.c_str(), value_f64, 2);

            continue;
        }

        if (dynamic_cast<TCLAP::SwitchArg*>(arg)) {
            TCLAP::SwitchArg & Sarg = *dynamic_cast<TCLAP::SwitchArg*>(arg);
            bool value = Sarg.getValue();
            h5_int64_t value_i64 = value? 1: 0;
            std::string name = Sarg.getName();
            H5WriteFileAttribInt64(_H5file, name.c_str(), &value_i64, 1);

            continue;
        }

        std::cerr << "**** can't save the argument '" << arg->getName()
                  << "' because it's type can't be determined ****\n"
                  << "**** going to exit ****" << std::endl;

        exit(0);

    }
}

void SaveH5Part::writeReducedH5FileAttributes(TCLAP::CmdLine & cmdl)
{
    std::set<std::string> includedArgs({"startstep", "numsteps", "savemodulus", "order", "debug"});
    std::list<TCLAP::Arg*> callArguments = cmdl.getArgList();
    for (auto arg: callArguments) {
        if (includedArgs.find(arg->getName()) == includedArgs.end()) continue;
        /*
        if (dynamic_cast<TCLAP::ValueArg<std::string>*>(arg)) {
            TCLAP::ValueArg<std::string> & Varg = *dynamic_cast<TCLAP::ValueArg<std::string>*>(arg);
            std::string value = Varg.getValue();
            std::string name = Varg.getName();
            H5ChangeFileAttribString(_H5file, name.c_str(), value.c_str());

            continue;
        }

        if (dynamic_cast<TCLAP::ValueArg<double>*>(arg)) {
            TCLAP::ValueArg<double> & Varg = *dynamic_cast<TCLAP::ValueArg<double>*>(arg);
            h5_float64_t value = Varg.getValue();
            std::string name = Varg.getName();
            H5ChangeFileAttribFloat64(_H5file, name.c_str(), &value, 1);

            continue;
        }
        */

        if (dynamic_cast<TCLAP::ValueArg<int>*>(arg)) {
            TCLAP::ValueArg<int> & Varg = *dynamic_cast<TCLAP::ValueArg<int>*>(arg);
            h5_int64_t value = Varg.getValue();
            std::string name = Varg.getName();
            H5ChangeFileAttribInt64(_H5file, name.c_str(), &value, 1);

            continue;
        }

        if (dynamic_cast<TCLAP::ValueArg<unsigned>*>(arg)) {
            TCLAP::ValueArg<unsigned> & Varg = *dynamic_cast<TCLAP::ValueArg<unsigned>*>(arg);
            h5_int64_t value = Varg.getValue();
            std::string name = Varg.getName();
            H5ChangeFileAttribInt64(_H5file, name.c_str(), &value, 1);

            continue;
        }

        /*
        if (dynamic_cast<TCLAP::ValueArg<Vect2D>*>(arg)) {
            TCLAP::ValueArg<Vect2D> & Varg = *dynamic_cast<TCLAP::ValueArg<Vect2D>*>(arg);
            Vect2D value = Varg.getValue();
            h5_float64_t value_f64[] = {value.v[0], value.v[1]};
            std::string name = Varg.getName();
            H5ChangeFileAttribFloat64(_H5file, name.c_str(), value_f64, 2);

            continue;
        }
        */

        if (dynamic_cast<TCLAP::SwitchArg*>(arg)) {
            TCLAP::SwitchArg & Sarg = *dynamic_cast<TCLAP::SwitchArg*>(arg);
            bool value = Sarg.getValue();
            h5_int64_t value_i64 = value? 1: 0;
            std::string name = Sarg.getName();
            H5ChangeFileAttribInt64(_H5file, name.c_str(), &value_i64, 1);

            continue;
        }

        std::cerr << "**** can't save the argument '" << arg->getName()
                  << "' because it's type can't be determined ****\n"
                  << "**** going to exit ****" << std::endl;

        exit(0);

    }
}
