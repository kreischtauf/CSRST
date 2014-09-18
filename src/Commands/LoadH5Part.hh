/***************************************************************************
                            LoadH5Part.hh
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

#ifndef LOADH5PARTCOMMAND_HH
#define LOADH5PARTCOMMAND_HH

#include "mpi.h"
#include "hdf5.h"
#include "H5hut.h"
#include "../defs.hh"

class PartBunch;

class LoadH5Part {
public:
    LoadH5Part(const std::string & fname,
               TCLAP::CmdLine & cmdl);

    ~LoadH5Part();

    void restoreFields(VField_Edge_t & EFD,
                       VField_Cell_t & HFD,
                       VField_Edge_t & JFD,
                       double & t);

    void restoreBunch(std::vector<Vector_t> & R,
                      std::vector<Vector_t> & oldR,
                      std::vector<Vector_t> & P,
                      std::vector<double> & Q);
private:

    void readH5FileAttributes(TCLAP::CmdLine & cmdl);

    void closeFile();

    void openFile();

    h5_file* _H5file;

    std::string _fileName;

    Timings::TimerRef _inputTimer;
};

#endif
