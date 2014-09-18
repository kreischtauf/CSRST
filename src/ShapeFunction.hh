/***************************************************************************
                           ShapeFunction.hh
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

#ifndef SHAPEFUNCTION_HH
#define SHAPEFUNCTION_HH

#include "defs.hh"
#include "ZerothOrderShapeFunction.hh"
#include "FirstOrderShapeFunction.hh"
#include "SecondOrderShapeFunction.hh"
#include "Timings.hh"
#include "FieldPatch.hh"

class PartBunch;

template<class IntOp>
class ShapeFunction
{
public:
    ShapeFunction()
    { }

    void getCurrentDensityX(VField_Edge_t& JFD,
                            const FieldPatch<double> & Jx,
                            const PartBunch & bunch,
                            const double & dt,
                            const size_t & from,
                            const size_t & to) const
    {
        ERRORMSG("Shape function not defined" << endl);
        return;
    }

    void getCurrentDensityY(VField_Edge_t& JFD,
                            const FieldPatch<double> & Jy,
                            const PartBunch & bunch,
                            const double & dt,
                            const size_t & from,
                            const size_t & to) const
    {
        ERRORMSG("Shape function not defined" << endl);
        return;
    }

    FieldPatch<double> getCurrentDensityBaseX(VField_Edge_t& JFD,
                                              const PartBunch & bunch,
                                              const double & dt,
                                              const size_t & from,
                                              const size_t & to) const
    {
        ERRORMSG("Shape function not defined" << endl);
        return;
    }

    FieldPatch<double> getCurrentDensityBaseY(VField_Edge_t& JFD,
                                              const PartBunch & bunch,
                                              const double & dt,
                                              const size_t & from,
                                              const size_t & to) const
    {
        ERRORMSG("Shape function not defined" << endl);
        return;
    }

    void getCurrentDensity(VField_Edge_t& JFD,
                           const PartBunch & bunch,
                           const double & dt)
    {
        ERRORMSG("Shape function not defined" << endl);
        return;
    }

    void getChargeDensity(std::map<coord_t, double> & rho,
                          const std::vector<Particle> & particles,
                          const Vector_t & dx,
                          const Vector_t & origin)
    {
        ERRORMSG("Shape function not defined" << endl);
        return;
    }

    void getChargeDensity(std::map<coord_t, double> & rho,
                          const PartBunch & bunch,
                          const Vector_t & dx,
                          const Vector_t & origin)
    {
        ERRORMSG("Shape function not defined" << endl);
        return;
    }

    static
    size_t getNodeBoundaryMargin()
    {
        ERRORMSG("Shape function not defined" << endl);
        return 0;
    }

    static
    NDIndex<DIM> getExtraMarginEdge()
    {
        ERRORMSG("Shape function not defined" << endl);
        return NDIndex<DIM>(Index(0),Index(0));
    }

    static
    Index getExtraMarginCell()
    {
        ERRORMSG("Shape function not defined" << endl);
        return Index(0,0);
    }

    static
    Index getExtraMarginVert()
    {
        ERRORMSG("Shape function not defined" << endl);
        return Index(0,0);
    }

    static
    NDIndex<DIM> getExtraMarginCurrent()
    {
        ERRORMSG("Shape function not defined" << endl);
        return NDIndex<DIM>(Index(0),Index(0));
    }
};

template<>
class ShapeFunction<IntCIC>
{
public:
    ShapeFunction()
    {
        _current_timer = Timings::getTimer("0st order::current");
        _charge_timer = Timings::getTimer("0st order::charge");
        _communication_timer = Timings::getTimer("0st order::comm");
    }

    void getCurrentDensityX(VField_Edge_t& JFD,
                            FieldPatch<double> & Jx,
                            const PartBunch & bunch,
                            const double & dt) const
    {
        Timings::startTimer(_current_timer);
        ZerothOrderShapeFunction::getCurrentDensityImplX(JFD, Jx, bunch, dt, _current_timer, _communication_timer);
        Timings::stopTimer(_current_timer);
    }


    void getCurrentDensityY(VField_Edge_t& JFD,
                            FieldPatch<double> & Jy,
                            const PartBunch & bunch,
                            const double & dt) const
    {
        Timings::startTimer(_current_timer);
        ZerothOrderShapeFunction::getCurrentDensityImplY(JFD, Jy, bunch, dt, _current_timer, _communication_timer);
        Timings::stopTimer(_current_timer);
    }

    FieldPatch<double> getCurrentDensityBaseX(VField_Edge_t& JFD,
                                              const PartBunch & bunch,
                                              const double & dt,
                                              const size_t & from,
                                              const size_t & to) const
    {
        Timings::startTimer(_current_timer);
        FieldPatch<double> Jx = ZerothOrderShapeFunction::getCurrentDensityBaseImplX(JFD,
                                                                                     bunch,
                                                                                     dt,
                                                                                     _current_timer,
                                                                                     from,
                                                                                     to);
        Timings::stopTimer(_current_timer);
        return Jx;
    }

    FieldPatch<double> getCurrentDensityBaseY(VField_Edge_t& JFD,
                                              const PartBunch & bunch,
                                              const double & dt,
                                              const size_t & from,
                                              const size_t & to) const
    {
        Timings::startTimer(_current_timer);
        FieldPatch<double> Jy = ZerothOrderShapeFunction::getCurrentDensityBaseImplY(JFD,
                                                                                     bunch,
                                                                                     dt,
                                                                                     _current_timer,
                                                                                     from,
                                                                                     to);
        Timings::stopTimer(_current_timer);
        return Jy;
    }

    void getChargeDensity(SField_Vert_t& rho,
                          const PartBunch & bunch) const
    {
        Timings::startTimer(_charge_timer);
        ZerothOrderShapeFunction::getChargeDensityImpl(rho, bunch, _charge_timer, _communication_timer);
        Timings::stopTimer(_charge_timer);
    }

    void getChargeDensityDiff(SField_Vert_t& rho,
                              const PartBunch & bunch) const
    {
        Timings::startTimer(_charge_timer);
        ZerothOrderShapeFunction::getChargeDensityDiffImpl(rho, bunch, _charge_timer, _communication_timer);
        Timings::stopTimer(_charge_timer);
    }

    static
    void getNumberDensity(SField_Cell_t& dens,
                          const PartBunch & bunch)
    {
        Timings::TimerRef comm_timer;
        comm_timer = Timings::getTimer("0st order::comm");
        ZerothOrderShapeFunction::getNumberDensityImpl(dens, bunch, comm_timer);
    }

    static
    size_t getNodeBoundaryMargin()
    {
        return ZerothOrderShapeFunction::getNodeBoundaryMargin();
    }

    static
    NDIndex<DIM> getExtraMarginEdge()
    {
        return ZerothOrderShapeFunction::getExtraMarginEdge();
    }

    static
    Index getExtraMarginCell()
    {
        return ZerothOrderShapeFunction::getExtraMarginCell();
    }

    static
    Index getExtraMarginVert()
    {
        return ZerothOrderShapeFunction::getExtraMarginVert();
    }

    static
    NDIndex<DIM> getExtraMarginCurrent()
    {
        return ZerothOrderShapeFunction::getExtraMarginCurrent();
    }

private:
    Timings::TimerRef _current_timer;
    Timings::TimerRef _charge_timer;
    Timings::TimerRef _communication_timer;
};

template<>
class ShapeFunction<IntS1>
{
public:
    ShapeFunction()
    {
        _current_timer = Timings::getTimer("1st order::current");
        _charge_timer = Timings::getTimer("1st order::charge");
        _communication_timer = Timings::getTimer("1st order::comm");
    }

    void getCurrentDensityX(VField_Edge_t& JFD,
                            const FieldPatch<double> & Jx,
                            const PartBunch & bunch,
                            const double & dt,
                            const size_t & from,
                            const size_t & to) const
    {
        Timings::startTimer(_current_timer);
        FirstOrderShapeFunction::getCurrentDensityImplX(JFD, Jx, bunch, dt, _current_timer, _communication_timer, from, to);
        Timings::stopTimer(_current_timer);
    }

    void getCurrentDensityY(VField_Edge_t& JFD,
                            const FieldPatch<double> & Jy,
                            const PartBunch & bunch,
                            const double & dt,
                            const size_t & from,
                            const size_t & to) const
    {
        Timings::startTimer(_current_timer);
        FirstOrderShapeFunction::getCurrentDensityImplY(JFD, Jy, bunch, dt, _current_timer, _communication_timer, from, to);
        Timings::stopTimer(_current_timer);
    }

    FieldPatch<double> getCurrentDensityBaseX(VField_Edge_t& JFD,
                                              const PartBunch & bunch,
                                              const double & dt,
                                              const size_t & from,
                                              const size_t & to) const
    {
        Timings::startTimer(_current_timer);
        FieldPatch<double> Jx = FirstOrderShapeFunction::getCurrentDensityBaseImplX(JFD,
                                                                                    bunch,
                                                                                    dt,
                                                                                    _current_timer,
                                                                                    from,
                                                                                    to);
        Timings::stopTimer(_current_timer);
        return Jx;
    }

    FieldPatch<double> getCurrentDensityBaseY(VField_Edge_t& JFD,
                                              const PartBunch & bunch,
                                              const double & dt,
                                              const size_t & from,
                                              const size_t & to) const
    {
        Timings::startTimer(_current_timer);
        FieldPatch<double> Jy = FirstOrderShapeFunction::getCurrentDensityBaseImplY(JFD,
                                                                                    bunch,
                                                                                    dt,
                                                                                    _current_timer,
                                                                                    from,
                                                                                    to);
        Timings::stopTimer(_current_timer);
        return Jy;
    }

    void getChargeDensity(SField_Vert_t& rho,
                          const PartBunch & bunch) const
    {
        Timings::startTimer(_charge_timer);
        FirstOrderShapeFunction::getChargeDensityImpl(rho, bunch, _charge_timer, _communication_timer);
        Timings::stopTimer(_charge_timer);
    }

    void getChargeDensityDiff(SField_Vert_t& rho,
                              const PartBunch & bunch) const
    {
        Timings::startTimer(_charge_timer);
        FirstOrderShapeFunction::getChargeDensityDiffImpl(rho, bunch, _charge_timer, _communication_timer);
        Timings::stopTimer(_charge_timer);
    }

    static
    size_t getNodeBoundaryMargin()
    {
        return FirstOrderShapeFunction::getNodeBoundaryMargin();
    }

    static
    NDIndex<DIM> getExtraMarginEdge()
    {
        return FirstOrderShapeFunction::getExtraMarginEdge();
    }

    static
    Index getExtraMarginCell()
    {
        return FirstOrderShapeFunction::getExtraMarginCell();
    }

    static
    Index getExtraMarginVert()
    {
        return FirstOrderShapeFunction::getExtraMarginVert();
    }

    static
    NDIndex<DIM> getExtraMarginCurrent()
    {
        return FirstOrderShapeFunction::getExtraMarginCurrent();
    }

private:
    Timings::TimerRef _current_timer;
    Timings::TimerRef _charge_timer;
    Timings::TimerRef _communication_timer;
};

template<>
class ShapeFunction<IntS2>
{
public:
    ShapeFunction()
    {
        _current_timer = Timings::getTimer("2nd order::current");
        _charge_timer = Timings::getTimer("2nd order::charge");
        _communication_timer = Timings::getTimer("2nd order::comm");
    }

    void getCurrentDensityX(VField_Edge_t& JFD,
                            const FieldPatch<double> & Jx,
                            const PartBunch & bunch,
                            const double & dt,
                            const size_t & from,
                            const size_t & to) const
    {
        Timings::startTimer(_current_timer);
        SecondOrderShapeFunction::getCurrentDensityImplX(JFD, Jx, bunch, dt, _current_timer, _communication_timer, from, to);
        Timings::stopTimer(_current_timer);
    }

    void getCurrentDensityY(VField_Edge_t& JFD,
                            const FieldPatch<double> & Jy,
                            const PartBunch & bunch,
                            const double & dt,
                            const size_t & from,
                            const size_t & to) const
    {
        Timings::startTimer(_current_timer);
        SecondOrderShapeFunction::getCurrentDensityImplY(JFD, Jy, bunch, dt, _current_timer, _communication_timer, from, to);
        Timings::stopTimer(_current_timer);
    }

    FieldPatch<double> getCurrentDensityBaseX(VField_Edge_t& JFD,
                                              const PartBunch & bunch,
                                              const double & dt,
                                              const size_t & from,
                                              const size_t & to) const
    {
        Timings::startTimer(_current_timer);
        FieldPatch<double> Jx = SecondOrderShapeFunction::getCurrentDensityBaseImplX(JFD,
                                                                                     bunch,
                                                                                     dt,
                                                                                     _current_timer,
                                                                                     from,
                                                                                     to);
        Timings::stopTimer(_current_timer);
        return Jx;
    }

    FieldPatch<double> getCurrentDensityBaseY(VField_Edge_t& JFD,
                                              const PartBunch & bunch,
                                              const double & dt,
                                              const size_t & from,
                                              const size_t & to) const
    {
        Timings::startTimer(_current_timer);
        FieldPatch<double> Jy = SecondOrderShapeFunction::getCurrentDensityBaseImplY(JFD,
                                                                                     bunch,
                                                                                     dt,
                                                                                     _current_timer,
                                                                                     from,
                                                                                     to);
        Timings::stopTimer(_current_timer);
        return Jy;
    }

    void getChargeDensity(SField_Vert_t& rho,
                          const PartBunch & bunch) const
    {
        Timings::startTimer(_charge_timer);
        SecondOrderShapeFunction::getChargeDensityImpl(rho, bunch, _charge_timer, _communication_timer);
        Timings::stopTimer(_charge_timer);
    }

    void getChargeDensityDiff(SField_Vert_t& rho,
                              const PartBunch & bunch) const
    {
        Timings::startTimer(_charge_timer);
        SecondOrderShapeFunction::getChargeDensityDiffImpl(rho, bunch, _charge_timer, _communication_timer);
        Timings::stopTimer(_charge_timer);
    }

    static
    size_t getNodeBoundaryMargin()
    {
        return SecondOrderShapeFunction::getNodeBoundaryMargin();
    }

    static
    NDIndex<DIM> getExtraMarginEdge()
    {
        return SecondOrderShapeFunction::getExtraMarginEdge();
    }

    static
    Index getExtraMarginCell()
    {
        return SecondOrderShapeFunction::getExtraMarginCell();
    }

    static
    Index getExtraMarginVert()
    {
        return SecondOrderShapeFunction::getExtraMarginVert();
    }

    static
    NDIndex<DIM> getExtraMarginCurrent()
    {
        return SecondOrderShapeFunction::getExtraMarginCurrent();
    }

private:
    Timings::TimerRef _communication_timer;
    Timings::TimerRef _current_timer;
    Timings::TimerRef _charge_timer;
};

#endif
