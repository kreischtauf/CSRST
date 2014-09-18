/***************************************************************************
                          IntegrateElong.hh
                         -------------------
    begin                : Thu Jan 30 2014
    copyright            : (C) 2014 by Christof Kraus
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

#ifndef INTEGRATEELONG_HH
#define INTEGRATEELONG_HH

#include <fstream>
#include "../defs.hh"
#include "../utils.hh"

class PartBunch;
class PartBunchState;

class IntegrateElongCommand {
public:
    IntegrateElongCommand(const VField_Edge_t & mmEFD,
                          const double & dt);

    ~IntegrateElongCommand();

    void execute(const VField_Edge_t & mmEFD,
                 const PartBunchState & state,
                 Inform & msg);

    void write(std::string fname,
               const VField_Edge_t & EFD,
               const PartBunchState & state) const;

    void reset();

private:

    static
    void getPartInside(Vector_t & begin,
                       Vector_t & end,
                       const VField_Edge_t & mmEFD);

    size_t _numSamples;
    double _dt;
    std::vector<Utils::KahanAccumulation> _integral;
    Timings::TimerRef _saveTimer;
};

inline
void IntegrateElongCommand::reset()
{
    std::fill(_integral.begin(), _integral.end(), Utils::KahanAccumulation());
}
#endif
