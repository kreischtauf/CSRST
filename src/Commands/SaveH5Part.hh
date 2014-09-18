/***************************************************************************
                            SaveH5Part.hh
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

#ifndef SAVEH5PARTCOMMAND_HH
#define SAVEH5PARTCOMMAND_HH

#include "hdf5.h"
#include "H5hut.h"

#include "../defs.hh"

class PartBunch;

class SaveH5Part {
public:
    SaveH5Part(TCLAP::CmdLine & cmdl,
               const std::string & filename,
               const h5_int32_t & flag = H5_O_RDWR);

    ~SaveH5Part();

    void write(const VField_Edge_t & EFD,
               const VField_Cell_t & HFD,
               const VField_Edge_t & JFD,
               const PartBunch & bunch,
               const double & t,
               Inform & msg);

private:
    void writeH5FileAttributes(TCLAP::CmdLine & cmdl);

    void writeReducedH5FileAttributes(TCLAP::CmdLine & cmdl);

    void closeFile();

    void openFile(const h5_int32_t & flag);

    h5_file_t* _H5file;

    h5_float64_t* _x_data;
    h5_float64_t* _y_data;
    h5_float64_t* _z_data;

    h5_float64_t* _part_data;

    h5_int64_t _H5call;

    unsigned int _prev_nloc;
    int _domain_size;

    std::string _fileName;
    Timings::TimerRef _outputTimer;
};

#endif
