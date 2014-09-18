/***************************************************************************
                              Timings.hh
                         -------------------
    begin                : Mon Oct 11 2013
    copyright            : (C) 2013 by Christof Kraus
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

#ifndef TIMINGS_HH
#define TIMINGS_HH

#include "Ippl.h"

#include <iomanip>
#define noTIMINGSDEBUGOUTPUT

#ifdef TIMINGSDEBUGOUTPUT
#include <boost/timer/timer.hpp>
#endif

class Timings {
public:
    typedef IpplTimings::TimerRef TimerRef;
    static TimerRef getTimer(const char *);
    static void startTimer(TimerRef);
    static void stopTimer(TimerRef);
#ifdef TIMINGSDEBUGOUTPUT
private:
    static TimerRef _current;
    static boost::timer::cpu_timer _mainTimer;
    static std::map<TimerRef, std::string> _ref2Name;
#endif
};

inline
Timings::TimerRef Timings::getTimer(const char *nm)
{
    TimerRef ret = IpplTimings::getTimer(nm);

#ifdef TIMINGSDEBUGOUTPUT
    std::string s(nm);
    auto loc = _ref2Name.find(ret);
    if (loc == _ref2Name.end()) {
        _ref2Name.insert(std::pair<TimerRef, std::string>(ret, s));
    }
#endif
    return ret;
}

inline
void Timings::startTimer(TimerRef t)
{
    IpplTimings::startTimer(t);

#ifdef TIMINGSDEBUGOUTPUT
    Inform msg("Timings");
    std::string newName = _ref2Name[t];
    if (_current != -1U) {
        std::string currentName = _ref2Name[_current];
        msg << currentName << " and " << newName
            << " are overlaping" << endl;
    } else {
        if (newName != "ipplToyFDTDTimer")
            _current = t;
    }
    msg << "starting " << newName;
    boost::timer::cpu_times const elapsed_times(_mainTimer.elapsed());
    boost::timer::nanosecond_type const elapsed(elapsed_times.system
                                                + elapsed_times.user);
    msg << " at " << std::setprecision(5) << elapsed*1e-9 << endl;
#endif
    return;
}

inline
void Timings::stopTimer(TimerRef t)
{
    IpplTimings::stopTimer(t);

#ifdef TIMINGSDEBUGOUTPUT
    Inform msg("Timings");
    std::string oldName = _ref2Name[t];
    if (oldName != "ipplToyFDTDTimer") {
        boost::timer::cpu_times const elapsed_times(_mainTimer.elapsed());
        boost::timer::nanosecond_type const elapsed(elapsed_times.system
                                                    + elapsed_times.user);
        msg << "ending " << oldName << " at " << std::setprecision(5) << elapsed*1e-9 << endl;

        if (_current == t)
            _current = -1U;
    }
#endif
    return;
}

#endif
