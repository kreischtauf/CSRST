/***************************************************************************
                           SaveEMFields.hh
                         -------------------
    begin                : Tue Jan 14 2014
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

#ifndef SAVEEMFIELDS_HH
#define SAVEEMFIELDS_HH

#include <fstream>
#include "../defs.hh"

class PartBunch;
class PartBunchState;

class SaveEMFieldsCommand {
public:
    SaveEMFieldsCommand();

    ~SaveEMFieldsCommand();

    void execute(const VField_Edge_t & EFD,
                 const VField_Cell_t & HFD,
                 const VField_Edge_t & JFD,
                 const PartBunchState & state,
                 Inform & msg);


private:

    static
    void getPartInside(Vector_t & begin,
                       Vector_t & end,
                       const VField_Edge_t & EFD);

    Timings::TimerRef _saveTimer;
    unsigned int _iteration;
};

#endif
