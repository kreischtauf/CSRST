/***************************************************************************
                            TETMSolver.cpp
                         -------------------
    begin                : Mon Aug 29 2011
    copyright            : (C) 2011 by Christof Kraus
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

#include "mpi.h"
#include "TETMSolver.hh"
#include "BoundaryCell.hh"
#include "utils.hh"

extern std::ofstream dbg;

TETMSolver::TETMSolver(VField_t & AlphaE,
                       VField_t & AlphaH,
                       PartBunch & bunch,
                       Inform & msg,
                       const VField_Edge_t & EFD,
                       const SolverPrecision_t & dt):
    _pbunch(bunch),
    _alpha_e(AlphaE),
    _alpha_h(AlphaH),
    _dt(dt),
    _me(Ippl::myNode())
{
    std::vector<int> nodes;

    _fieldUpdateTimer = Timings::getTimer("field update");
    _exchangeTimer =  Timings::getTimer("field exchange");
    FieldLayout<DIM> & FL = EFD.getLayout();
    _gDom = FL.getDomain();

    std::vector<NDIndex<DIM> > lDoms;
    Utils::getLocalDomains(FL, lDoms);

    _Nx = _gDom[0].last() - _gDom[0].first() + 1;
    _Ny = _gDom[1].last() - _gDom[1].first() + 1;
    _nx = lDoms[_me][0].last() - lDoms[_me][0].first() + 1;
    _ny = lDoms[_me][1].last() - lDoms[_me][1].first() + 1;


    initMaps(lDoms[_me]);
}

TETMSolver::~TETMSolver()
{ }

void TETMSolver::initMaps(const NDIndex<DIM> & lDom)
{
#ifndef noREFERENCERUN
    _alpha_e = 0.0;
    _alpha_h = 0.0;
#else
    _alpha_e.fillGuardCells();
    _alpha_h.fillGuardCells();
    NDIndex<DIM> elem;

    for (int i = lDom[0].first(); i <= lDom[0].last(); ++ i) {
        elem[0] = Index(i,i);
        for (int j = lDom[1].first(); j <= lDom[1].last(); ++ j) {
            elem[1] = Index(j,j);
            Vector_t ae = _alpha_e.localElement(elem);
            Vector_t ah = _alpha_h.localElement(elem);
            _alpha_e.localElement(elem) = Vector_t(dteps * ae[0], dteps * ae[1]);
            // sigma^{*} / mu_0 \stackrel{!}{=} sigma / epsilon_0
            // therefore here also multiply with dteps
            _alpha_h.localElement(elem) = Vector_t(dtmu * ah[0], dtmu * ah[1]);
        }
    }
#endif

    _alpha_e.fillGuardCells();
    _alpha_h.fillGuardCells();
}

Vector_t TETMSolver::updateHz(const VField_Edge_t& EFD,
                              const NDIndex<DIM>& elem,
                              const BoundaryCell bcell) const {
    Vector_t hr(EFD.get_mesh().get_meshSpacing(0), EFD.get_mesh().get_meshSpacing(1));
    double etax = _dt / (hr(0) * Physics::mu_0);
    double etay = _dt / (hr(1) * Physics::mu_0);
    NDIndex<DIM> neighbourX(elem[0] + 1, elem[1]);
    NDIndex<DIM> neighbourY(elem[0], elem[1] + 1);

    Vector_t HFDUpdate(0.0);
    HFDUpdate(0) = -0.5 * etax * (bcell.lambdaToNeighbours_m(1) * EFD.localElement(neighbourX)(1) -
                                  bcell.lambda_m(1) * EFD.localElement(elem)(1)) / bcell.area_m;
    HFDUpdate(1) = +0.5 * etay * (bcell.lambdaToNeighbours_m(0) * EFD.localElement(neighbourY)(0) -
                                  bcell.lambda_m(0) * EFD.localElement(elem)(0)) / bcell.area_m;

    return HFDUpdate;
}

double TETMSolver::updateEx(const VField_Cell_t& HFD,
                            const VField_Edge_t& JFD,
                            const NDIndex<DIM>& elem) const {
    double hy = HFD.get_mesh().get_meshSpacing(1);
    double dteps = _dt / Physics::epsilon_0;
    NDIndex<DIM> neighbourY(elem[0], elem[1] - 1);

    Vector_t EFDXUpdate = (HFD.localElement(elem) -
                           HFD.localElement(neighbourY)) / hy;
    return dteps * (EFDXUpdate(0) + EFDXUpdate(1) - JFD.localElement(elem)(0));
}

double TETMSolver::updateEy(const VField_Cell_t& HFD,
                            const VField_Edge_t& JFD,
                            const NDIndex<DIM>& elem) const {
    double hx = HFD.get_mesh().get_meshSpacing(0);
    double dteps = _dt / Physics::epsilon_0;
    NDIndex<DIM> neighbourX(elem[0] - 1, elem[1]);

    Vector_t EFDYUpdate = -(HFD.localElement(elem) -
                            HFD.localElement(neighbourX)) / hx;
    return dteps * (EFDYUpdate(0) + EFDYUpdate(1) - JFD.localElement(elem)(1));
}

void outputNodeBoundary(VField_Edge_t & fd)
{
    NDIndex<DIM> ldom = fd.getLayout().getLocalNDIndex();
    NDIndex<DIM> gdom = fd.getLayout().getDomain();

    dbg << ldom << std::endl;

    std::vector<NDIndex<DIM> > boundaries;
    if (ldom[0].first() != gdom[0].first()) {
        NDIndex<DIM> dom;
        dom[0] = Index(ldom[0].first() - 1, ldom[0].first());
        dom[1] = ldom[1];
        boundaries.push_back(dom);
    }
    if (ldom[0].last() != gdom[0].last()) {
        NDIndex<DIM> dom;
        dom[0] = Index(ldom[0].last(), ldom[0].last() + 1);
        dom[1] = ldom[1];
        boundaries.push_back(dom);
    }
    if (ldom[1].first() != gdom[1].first()) {
        NDIndex<DIM> dom;
        dom[0] = ldom[0];
        dom[1] = Index(ldom[1].first() - 1, ldom[1].first());
        boundaries.push_back(dom);
    }
    if (ldom[1].last() != gdom[1].last()) {
        NDIndex<DIM> dom;
        dom[0] = ldom[0];
        dom[1] = Index(ldom[1].last(), ldom[1].last() + 1);
        boundaries.push_back(dom);
    }
    for (unsigned int k = 0; k < boundaries.size(); ++ k) {
        NDIndex<DIM> elem;
        dbg << boundaries[k] << "\n" << std::endl;
        for (int j = boundaries[k][1].first(); j <= boundaries[k][1].last(); ++ j) {
            elem[1] = Index(j,j);
            for (int i = boundaries[k][0].first(); i <= boundaries[k][0].last(); ++ i) {
                elem[0] = Index(i,i);
                Vector_t val = fd.localElement(elem);
                dbg << std::setw(6) << i
                    << std::setw(6) << j
                    << std::setw(14) << val(0)
                    << std::setw(14) << val(1) << "\n";
            }
        }
        dbg << std::endl;
    }
}

void outputNodeBoundary(VField_Cell_t & fd, VField_Edge_t & fd2, const double & etax, const double & etay)
{
    NDIndex<DIM> ldom = fd.getLayout().getLocalNDIndex();
    NDIndex<DIM> gdom = fd.getLayout().getDomain();
    double Z_0 = sqrt(Physics::mu_0 / Physics::epsilon_0);

    std::vector<NDIndex<DIM> > boundaries;
    if (ldom[0].first() != gdom[0].first()) {
        NDIndex<DIM> dom;
        dom[0] = Index(ldom[0].first(), ldom[0].first());
        dom[1] = ldom[1];
        boundaries.push_back(dom);
    }
    if (ldom[0].last() != gdom[0].last()) {
        NDIndex<DIM> dom;
        dom[0] = Index(ldom[0].last(), ldom[0].last());
        dom[1] = ldom[1];
        boundaries.push_back(dom);
    }
    if (ldom[1].first() != gdom[1].first()) {
        NDIndex<DIM> dom;
        dom[0] = ldom[0];
        dom[1] = Index(ldom[1].first(), ldom[1].first());
        boundaries.push_back(dom);
    }
    if (ldom[1].last() != gdom[1].last()) {
        NDIndex<DIM> dom;
        dom[0] = ldom[0];
        dom[1] = Index(ldom[1].last(), ldom[1].last());
        boundaries.push_back(dom);
    }
    for (unsigned int k = 0; k < boundaries.size(); ++ k) {
        NDIndex<DIM> elem, elem21, elem22;
        dbg << boundaries[k] << "\n" << std::endl;
        for (int j = boundaries[k][1].first(); j <= boundaries[k][1].last(); ++ j) {
            elem[1] = Index(j,j);
            elem21[1] = Index(j+1,j+1);
            elem22[1] = elem[1];
            for (int i = boundaries[k][0].first(); i <= boundaries[k][0].last(); ++ i) {
                elem[0] = Index(i,i);
                elem21[0] = elem[0];
                elem22[0] = Index(i+1,i+1);

                Vector_t val = fd.localElement(elem);
                Vector_t val20 = fd2.localElement(elem);
                Vector_t val21 = fd2.localElement(elem21);
                Vector_t val22 = fd2.localElement(elem22);
                dbg << std::setw(6) << i
                    << std::setw(6) << j
                    << std::setw(14) << val(0)
                    << std::setw(14) << val(1)
                    << std::setw(14) << - 0.5 * etax / Z_0 * (val22(1) - val20(1))
                    << std::setw(14) <<   0.5 * etay / Z_0 * (val21(0) - val20(0))
                    << "\n";
            }
        }
        dbg << std::endl;
    }
}
