/***************************************************************************
                            CalcEnergy.hh
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

#ifndef CALCENERGY_HH
#define CALCENERGY_HH

#include <fstream>
#include "../defs.hh"

class CalcEnergyCommand {
public:
    CalcEnergyCommand(const VField_Edge_t & EFD);

    ~CalcEnergyCommand();

    double execute(const VField_Edge_t & EFD,
                   const VField_Cell_t & HFD,
                   const VField_Edge_t & mmEFD,
                   const VField_Cell_t & mmHFD,
                   const Vektor<int, DIM> & start,
                   const double & t,
                   Inform & msg);


private:

    double _energy_EFD_old;
    double _area;
    // std::string _file_name;
    // std::ofstream _energy_out;
    Timings::TimerRef _executeTimer;
};

#endif
