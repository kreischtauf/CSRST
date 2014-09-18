/***************************************************************************
                     ZerothOrderShapeFunction.hh
                         -------------------
    begin                : Thu Feb 23 2012
    copyright            : (C) 2012 by Christof Kraus
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

#ifndef ZEROTHORDERSHAPEFUNCTION_HH
#define ZEROTHORDERSHAPEFUNCTION_HH

#include "defs.hh"

template<class c>
class FieldPatch;

class PartBunch;

class ZerothOrderShapeFunction
{
public:
    ZerothOrderShapeFunction()
    { }

    static
    void getCurrentDensityImplX(VField_Edge_t& JFD,
                                FieldPatch<double> & Jx,
                                const PartBunch & bunch,
                                const double & dt,
                                const Timings::TimerRef & currentTimer,
                                const Timings::TimerRef & commTimer);

    static
    void getCurrentDensityImplY(VField_Edge_t& JFD,
                                FieldPatch<double> & Jy,
                                const PartBunch & bunch,
                                const double & dt,
                                const Timings::TimerRef & currentTimer,
                                const Timings::TimerRef & commTimer);

    static
    FieldPatch<double> getCurrentDensityBaseImplX(VField_Edge_t& JFD,
                                                  const PartBunch & bunch,
                                                  const double & dt,
                                                  const Timings::TimerRef & currentTimer,
                                                  const size_t & from,
                                                  const size_t & to);

    static
    FieldPatch<double> getCurrentDensityBaseImplY(VField_Edge_t& JFD,
                                                  const PartBunch & bunch,
                                                  const double & dt,
                                                  const Timings::TimerRef & currentTimer,
                                                  const size_t & from,
                                                  const size_t & to);

    static
    void getChargeDensityImpl(SField_Vert_t& rho,
                              const PartBunch & bunch,
                              const Timings::TimerRef & chargeTimer,
                              const Timings::TimerRef & commTimer);

    static
    void getChargeDensityDiffImpl(SField_Vert_t& rho,
                                  const PartBunch & bunch,
                                  const Timings::TimerRef & chargeTimer,
                                  const Timings::TimerRef & commTimer);

    static
    void getNumberDensityImpl(SField_Cell_t& dens,
                              const PartBunch & bunch,
                              const Timings::TimerRef & commTimer);

    static
    size_t getNodeBoundaryMargin();

    static
    NDIndex<DIM> getExtraMarginEdge();

    static
    Index getExtraMarginCell();

    static
    Index getExtraMarginVert();

    static
    NDIndex<DIM> getExtraMarginCurrent();
};

inline
size_t ZerothOrderShapeFunction::getNodeBoundaryMargin()
{
    return 0;
}

inline
NDIndex<DIM> ZerothOrderShapeFunction::getExtraMarginEdge()
{
    return NDIndex<DIM>(Index(-1,1),Index(0,1));
}

inline
Index ZerothOrderShapeFunction::getExtraMarginCell()
{
    return Index(-1,1);
}

inline
Index ZerothOrderShapeFunction::getExtraMarginVert()
{
    return Index(0,1);
}

inline
NDIndex<DIM> ZerothOrderShapeFunction::getExtraMarginCurrent()
{
    return NDIndex<DIM>(Index(0,0),Index(0,1));
}

#endif
