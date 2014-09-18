/***************************************************************************
                             DataSink.hh
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

#ifndef DATASINK_HH
#define DATASINK_HH

#include <fstream>
#include "../defs.hh"
#include "CalcEnergy.hh"

class PartBunch;
class PartBunchState;

class DataSinkCommand {
public:
    DataSinkCommand(const VField_Edge_t & EFD);

    ~DataSinkCommand();

    void writeSDDSHeader();

    void execute(const PartBunch & bunch,
                 const double & t,
                 const VField_Edge_t & EFD,
                 const VField_Cell_t & HFD,
                 const VField_Edge_t & mmEFD,
                 const VField_Cell_t & mmHFD,
                 const Vektor<int, DIM> & start,
                 PartBunchState & state,
                 Inform & msg);


private:
    CalcEnergyCommand _cec;
    std::ofstream _sink;
};

#endif
