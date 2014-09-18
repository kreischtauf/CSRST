/***************************************************************************
                     FirstOrderShapeFunction.hh
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

#ifndef FIRSTORDERSHAPEFUNCTION_HH
#define FIRSTORDERSHAPEFUNCTION_HH

#include "defs.hh"

template<class c>
class FieldPatch;

class PartBunch;

class FirstOrderShapeFunction
{
public:
    FirstOrderShapeFunction()
    { };

    static
    void getCurrentDensityImplX(VField_Edge_t& JFD,
                                const FieldPatch<double> & Jx,
                                const PartBunch & bunch,
                                const double & dt,
                                const Timings::TimerRef & currentTimer,
                                const Timings::TimerRef & commTimer,
                                const size_t & from,
                                const size_t & to);

    static
    void getCurrentDensityImplY(VField_Edge_t& JFD,
                                const FieldPatch<double> & Jy,
                                const PartBunch & bunch,
                                const double & dt,
                                const Timings::TimerRef & currentTimer,
                                const Timings::TimerRef & commTimer,
                                const size_t & from,
                                const size_t & to);

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
    size_t getNodeBoundaryMargin();

    static
    NDIndex<DIM> getExtraMarginEdge();

    static
    Index getExtraMarginCell();

    static
    Index getExtraMarginVert();

    static
    NDIndex<DIM> getExtraMarginCurrent();

private:
    static
    void getCurrentX(double J[],
                     const double & tau,
                     const Vector_t & a,
                     const Vector_t & b,
                     const Vector_t & qoverdxdt,
                     const Vector_t & delta);

    static
    void getCurrentY(double J[],
                     const double & tau,
                     const Vector_t & a,
                     const Vector_t & b,
                     const Vector_t & qoverdxdt,
                     const Vector_t & delta);
};


inline
size_t FirstOrderShapeFunction::getNodeBoundaryMargin()
{
    return 1;
}

inline
NDIndex<DIM> FirstOrderShapeFunction::getExtraMarginEdge()
{
    return NDIndex<DIM>(Index(-1,1),Index(-1,2));
}

inline
Index FirstOrderShapeFunction::getExtraMarginCell()
{
    return Index(-1,1);
}

inline
Index FirstOrderShapeFunction::getExtraMarginVert()
{
    return Index(-1,2);
}

inline
NDIndex<DIM> FirstOrderShapeFunction::getExtraMarginCurrent()
{
    return NDIndex<DIM>(Index(-1,1),Index(-1,2));
}

#endif
