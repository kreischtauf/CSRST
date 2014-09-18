/***************************************************************************
                     ZerothOrderShapeFunction.cpp
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

#include "ZerothOrderShapeFunction.hh"
#include "PartBunch.hh"
#include "FieldPatch.hh"
#include "Communicator.hh"
#include "utils.hh"

extern std::ofstream dbg;
#define DBGOUT dbg << "ZerothOrderShapeFunction.cpp: " << __LINE__ << "\t"

void ZerothOrderShapeFunction::getCurrentDensityImplX(VField_Edge_t& JFD,
                                                      FieldPatch<double> & myJx,
                                                      const PartBunch & bunch,
                                                      const double & dt,
                                                      const Timings::TimerRef & currentTimer,
                                                      const Timings::TimerRef & commTimer)
{
    NDIndex<DIM> elem;
    const NDIndex<DIM> currentAdd = getExtraMarginCurrent();
    const unsigned int myNode = Ippl::myNode();
    const unsigned int numNodes = Ippl::getNodes();

    const std::vector<size_t> & numLocalParticles = bunch.getNumLocalParticles();

    NDIndex<DIM> lPDom = bunch.getLocalPDomainInclOld(JFD.get_mesh());

    Timings::stopTimer(currentTimer);
    auto localPDomains = Communicator::getAllLocalDomains(lPDom);
    Timings::startTimer(currentTimer);

    NDIndex<DIM> tmp;
    unsigned int kk = 0;
    while (kk < numNodes && numLocalParticles[kk] == 0) ++ kk;
    tmp = localPDomains[kk];

    for (; kk < numNodes; ++ kk) {
        if (numLocalParticles[kk] == 0) continue;

        NDIndex<DIM> & dom = localPDomains[kk];
        for (unsigned int d = 0; d < DIM; ++ d) {
            int lower = std::min(tmp[d].first(), dom[d].first());
            int upper = std::max(tmp[d].last(), dom[d].last());
            tmp[d] = Index(lower, upper);

            dom[d] = Index(dom[d].first() + currentAdd[d].first(),
                           dom[d].last() + currentAdd[d].last());
        }
    }
    const NDIndex<DIM> gPDom = tmp;
    const_cast<PartBunch*>(&bunch)->setGlobalBounds(tmp);

    std::vector<NDIndex<DIM> > localFDomains;
    Utils::getLocalDomains(JFD.getLayout(), localFDomains);
    for (unsigned int k = 0; k < numNodes; ++ k) {
        NDIndex<DIM> & dom = localFDomains[k];
        for (unsigned int d = 0; d < DIM; ++ d) {
            int lower = std::max(dom[d].first(), gPDom[d].first() + currentAdd[d].first());
            int upper = std::min(dom[d].last(),  gPDom[d].last()  + currentAdd[d].last());
            int sign = lower <= upper? 1: upper - lower;
            dom[d] = Index(lower, upper, sign);
        }
    }
    const NDIndex<DIM> lDom = JFD.getLayout().getLocalNDIndex();
    const NDIndex<DIM> & lFDom = localFDomains[myNode];

    Timings::stopTimer(currentTimer);
    Communicator::communicateFields(myJx, localPDomains, localFDomains, commTimer);
    Timings::startTimer(currentTimer);

    JFD[lDom[0]][lDom[1]](0) -= JFD[lDom[0]][lDom[1]](0);

    for (int j = lFDom[1].first(); j <= lFDom[1].last(); ++ j) {
        elem[1] = Index(j, j);
        for (int i = lFDom[0].first(); i <= lFDom[0].last(); ++ i) {
            elem[0] = Index(i, i);
            JFD.localElement(elem)(0) = myJx(i, j);
        }
    }

    JFD.fillGuardCells();
}

void ZerothOrderShapeFunction::getCurrentDensityImplY(VField_Edge_t& JFD,
                                                      FieldPatch<double> & myJy,
                                                      const PartBunch & bunch,
                                                      const double & dt,
                                                      const Timings::TimerRef & currentTimer,
                                                      const Timings::TimerRef & commTimer)
{
    NDIndex<DIM> elem;
    const NDIndex<DIM> currentAdd = getExtraMarginCurrent();
    const unsigned int myNode = Ippl::myNode();
    const unsigned int numNodes = Ippl::getNodes();

    const std::vector<size_t> & numLocalParticles = bunch.getNumLocalParticles();

    NDIndex<DIM> lPDom = bunch.getLocalPDomainInclOld(JFD.get_mesh());

    Timings::stopTimer(currentTimer);
    auto localPDomains = Communicator::getAllLocalDomains(lPDom);
    Timings::startTimer(currentTimer);

    NDIndex<DIM> tmp;
    unsigned int kk = 0;
    while (kk < numNodes && numLocalParticles[kk] == 0) ++ kk;
    tmp = localPDomains[kk];

    for (; kk < numNodes; ++ kk) {
        if (numLocalParticles[kk] == 0) continue;

        NDIndex<DIM> & dom = localPDomains[kk];
        for (unsigned int d = 0; d < DIM; ++ d) {
            int lower = std::min(tmp[d].first(), dom[d].first());
            int upper = std::max(tmp[d].last(), dom[d].last());
            tmp[d] = Index(lower, upper);

            dom[d] = Index(dom[d].first() + currentAdd[1 - d].first(),
                           dom[d].last() + currentAdd[1 - d].last());
        }
    }
    const NDIndex<DIM> gPDom = tmp;
    const_cast<PartBunch*>(&bunch)->setGlobalBounds(tmp);

    std::vector<NDIndex<DIM> > localFDomains;
    Utils::getLocalDomains(JFD.getLayout(), localFDomains);
    for (unsigned int k = 0; k < numNodes; ++ k) {
        NDIndex<DIM> & dom = localFDomains[k];
        for (unsigned int d = 0; d < DIM; ++ d) {
            int lower = std::max(dom[d].first(), gPDom[d].first() + currentAdd[1 - d].first());
            int upper = std::min(dom[d].last(),  gPDom[d].last()  + currentAdd[1 - d].last());
            int sign = lower <= upper? 1: upper - lower;
            dom[d] = Index(lower, upper, sign);
        }
    }
    const NDIndex<DIM> lDom = JFD.getLayout().getLocalNDIndex();
    const NDIndex<DIM> & lFDom = localFDomains[myNode];

    Timings::stopTimer(currentTimer);
    Communicator::communicateFields(myJy, localPDomains, localFDomains, commTimer);
    Timings::startTimer(currentTimer);

    JFD[lDom[0]][lDom[1]](1) -= JFD[lDom[0]][lDom[1]](1);

    for (int j = lFDom[1].first(); j <= lFDom[1].last(); ++ j) {
        elem[1] = Index(j, j);
        for (int i = lFDom[0].first(); i <= lFDom[0].last(); ++ i) {
            elem[0] = Index(i, i);
            JFD.localElement(elem)(1) = myJy(i, j);
        }
    }

    JFD.fillGuardCells();
}

FieldPatch<double> ZerothOrderShapeFunction::getCurrentDensityBaseImplX(VField_Edge_t& JFD,
                                                                        const PartBunch & bunch,
                                                                        const double & dt,
                                                                        const Timings::TimerRef & currentTimer,
                                                                        const size_t & from,
                                                                        const size_t & to)
{
    const NDIndex<DIM> currentAdd = getExtraMarginCurrent();

    const Vector_t dx(JFD.get_mesh().get_meshSpacing(0),
                      JFD.get_mesh().get_meshSpacing(1));
    const Vector_t origin = JFD.get_mesh().get_origin();
    Vector_t Rmin, Rmax;
    bunch.getBoundsInclOld(Rmin, Rmax, from, to);

    NDIndex<DIM> lPDom;
    for (unsigned int d = 0; d < DIM; ++ d) {
        lPDom[d] = Index(std::floor((Rmin(d) - origin(d)) / dx(d)),
                         std::floor((Rmax(d) - origin(d)) / dx(d)));
    }
    NDIndex<DIM> dom;
    for (unsigned int d = 0; d < DIM; ++ d) {
        dom[d] = Index(lPDom[d].first() + currentAdd[d].first(),
                         lPDom[d].last() + currentAdd[d].last());
    }

    Utils::samePos uni(sqrt(dot(dx,dx)));
    const Vector_t oneoverdx = Vector_t(1.0) / dx;
    const Vector_t oneoverdxdt = Vector_t(1.0 / (dx(1) * dt), 1.0 / (dx(0) * dt));
    const double onehalf = 0.5;

    FieldPatch<double> myJx(dom);
    myJx.setOrigin(origin + Vector_t(dom[0].first() * dx(0),
                                     dom[1].first() * dx(1)));
    myJx.setSpacing(dx);

    const ParticleAttrib< Vector_t > & R = bunch.R;
    const ParticleAttrib< Vector_t > & oldR = bunch.oldR;
    const ParticleAttrib< double > & Q = bunch.Q;

    Utils::CompanionSorterItem xy_crossings[10];
    for (size_t i = from; i < to; ++ i) {
        const double & q = Q[i];

        size_t numXs = Utils::calcCellBorderCrossings(xy_crossings, oldR[i], R[i], dx, origin, uni);
        for (size_t j = 0; j < numXs - 1; ++ j) {
            const Vector_t xi = xy_crossings[j].subStep * oneoverdx;
            const Vector_t xf = xy_crossings[j+1].subStep * oneoverdx;
            const Vector_t Dx = q * (xf - xi) * oneoverdxdt;
            const Vector_t center = (xi + xf) * onehalf;
            int I = static_cast<int>(floor(center(0)));
            int J = static_cast<int>(floor(center(1)));

            myJx(I, J) += Dx(0) * (J + 1 - center(1));
            myJx(I, J + 1) += Dx(0) * (center(1) - J);
        }
    }

    return myJx;
}

FieldPatch<double> ZerothOrderShapeFunction::getCurrentDensityBaseImplY(VField_Edge_t& JFD,
                                                                        const PartBunch & bunch,
                                                                        const double & dt,
                                                                        const Timings::TimerRef & currentTimer,
                                                                        const size_t & from,
                                                                        const size_t & to)
{
    const NDIndex<DIM> currentAdd = getExtraMarginCurrent();

    Timings::stopTimer(currentTimer);
    const Vector_t dx(JFD.get_mesh().get_meshSpacing(0),
                      JFD.get_mesh().get_meshSpacing(1));
    const Vector_t origin = JFD.get_mesh().get_origin();
    Vector_t Rmin, Rmax;
    bunch.getBoundsInclOld(Rmin, Rmax, from, to);

    NDIndex<DIM> lPDom;
    for (unsigned int d = 0; d < DIM; ++ d) {
        lPDom[d] = Index(std::floor((Rmin(d) - origin(d)) / dx(d)),
                         std::floor((Rmax(d) - origin(d)) / dx(d)));
    }
    NDIndex<DIM> dom;
    for (unsigned int d = 0; d < DIM; ++ d) {
        dom[d] = Index(lPDom[d].first() + currentAdd[1 - d].first(),
                         lPDom[d].last() + currentAdd[1 - d].last());
    }

    Timings::startTimer(currentTimer);

    Utils::samePos uni(sqrt(dot(dx,dx)));
    const Vector_t oneoverdx = Vector_t(1.0) / dx;
    const Vector_t oneoverdxdt = Vector_t(1.0 / (dx(1) * dt), 1.0 / (dx(0) * dt));
    const double onehalf = 0.5;

    FieldPatch<double> myJy(dom);
    myJy.setOrigin(origin + Vector_t(dom[0].first() * dx(0),
                                     dom[1].first() * dx(1)));
    myJy.setSpacing(dx);

    const ParticleAttrib< Vector_t > & R = bunch.R;
    const ParticleAttrib< Vector_t > & oldR = bunch.oldR;
    const ParticleAttrib< double > & Q = bunch.Q;

    Utils::CompanionSorterItem xy_crossings[10];
    for (size_t i = from; i < to; ++ i) {
        const double & q = Q[i];

        size_t numXs = Utils::calcCellBorderCrossings(xy_crossings, oldR[i], R[i], dx, origin, uni);
        for (size_t j = 0; j < numXs - 1; ++ j) {
            const Vector_t xi = xy_crossings[j].subStep * oneoverdx;
            const Vector_t xf = xy_crossings[j+1].subStep * oneoverdx;
            const Vector_t Dx = q * (xf - xi) * oneoverdxdt;
            const Vector_t center = (xi + xf) * onehalf;
            int I = static_cast<int>(floor(center(0)));
            int J = static_cast<int>(floor(center(1)));

            myJy(I, J) += Dx(1) * (I + 1 - center(0));
            myJy(I + 1, J) += Dx(1) * (center(0) - I);
        }
    }

    return myJy;
}

void ZerothOrderShapeFunction::getChargeDensityImpl(SField_Vert_t& rho,
                                                    const PartBunch & bunch,
                                                    const Timings::TimerRef & chargeTimer,
                                                    const Timings::TimerRef & commTimer)
{
    NDIndex<DIM> elem;
    Mesh_t & mesh = rho.get_mesh();
    const Index vertAdd = getExtraMarginVert();
    const unsigned int myNode = Ippl::myNode();
    const unsigned int numNodes = Ippl::getNodes();
    const std::vector<size_t> & numLocalParticles = bunch.getNumLocalParticles();
    NDIndex<DIM> tmp = bunch.getLocalPDomain(mesh);
    auto localPDomains = Communicator::getAllLocalDomains(tmp);

    unsigned int kk = 0;
    while (kk < numNodes && numLocalParticles[kk] == 0) ++ kk;
    tmp = localPDomains[kk];

    for (; kk < numNodes; ++ kk) {
        if (numLocalParticles[kk] == 0) continue;

        NDIndex<DIM> & dom = localPDomains[kk];
        for (unsigned int d = 0; d < DIM; ++ d) {
            int lower = std::min(tmp[d].first(), dom[d].first());
            int upper = std::max(tmp[d].last(), dom[d].last());
            tmp[d] = Index(lower, upper);

            dom[d] = Index(dom[d].first() + vertAdd.first(),
                           dom[d].last() + vertAdd.last());
        }
    }
    const NDIndex<DIM> gPDom = tmp;
    const NDIndex<DIM> & lPDom = localPDomains[myNode];

    std::vector<NDIndex<DIM> > localFDomains;
    Utils::getLocalDomains(rho.getLayout(), localFDomains);
    for (unsigned int k = 0; k < numNodes; ++ k) {
        NDIndex<DIM> & dom = localFDomains[k];
        for (unsigned int d = 0; d < DIM; ++ d) {
            int lower = std::max(dom[d].first(), gPDom[d].first() + vertAdd.first());
            int upper = std::min(dom[d].last(),  gPDom[d].last()  + vertAdd.last());
            int sign = lower <= upper? 1: upper - lower;
            dom[d] = Index(lower, upper, sign);
        }
    }
    const NDIndex<DIM> & lFDom = localFDomains[myNode];
    const NDIndex<DIM> lDom = rho.getLayout().getLocalNDIndex();

    const Vector_t dx(mesh.get_meshSpacing(0),
                      mesh.get_meshSpacing(1));
    const Vector_t origin = mesh.get_origin();

    FieldPatch<double> myRho(lPDom);
    myRho.setOrigin(origin + Vector_t(lPDom[0].first() * dx(0),
                                      lPDom[1].first() * dx(1)));
    myRho.setSpacing(dx);

    const ParticleAttrib< Vector_t > & R = bunch.R;
    const ParticleAttrib< double > & Q = bunch.Q;

    for (size_t i = 0; i < bunch.getLocalNP(); ++ i) {
        const double & q = Q[i];
        const Vector_t xi((R[i](0) - origin(0)) / dx(0),
                          (R[i](1) - origin(1)) / dx(1));
        const int I = static_cast<int>(floor(xi(0)));
        const int J = static_cast<int>(floor(xi(1)));
        const Vector_t tau(xi(0) - I, xi(1) - J);
        myRho(I, J)         += (1.0 - tau(0)) * (1.0 - tau(1)) * q;
        myRho(I + 1, J)     += tau(0)         * (1.0 - tau(1)) * q;
        myRho(I, J + 1)     += (1.0 - tau(0)) * tau(1)         * q;
        myRho(I + 1, J + 1) += tau(0)         * tau(1)         * q;
    }

    Timings::stopTimer(chargeTimer);
    Communicator::communicateFields(myRho, localPDomains, localFDomains, commTimer, true);
    Timings::startTimer(chargeTimer);

    rho = 0.0;
    for (int j = lFDom[1].first(); j <= lFDom[1].last(); ++ j) {
        elem[1] = Index(j, j);
        for (int i = lFDom[0].first(); i <= lFDom[0].last(); ++ i) {
            elem[0] = Index(i, i);
            rho.localElement(elem) += myRho(i,j) / (dx(0) * dx(1));
        }
    }

    rho.fillGuardCells();
}

void ZerothOrderShapeFunction::getChargeDensityDiffImpl(SField_Vert_t& rho,
                                                        const PartBunch & bunch,
                                                        const Timings::TimerRef & chargeTimer,
                                                        const Timings::TimerRef & commTimer)
{
    NDIndex<DIM> elem;
    const NDIndex<DIM> lDom = rho.getLayout().getLocalNDIndex();
    const NDIndex<DIM> & lPDom = bunch.getLocalPDomainInclOld(rho.get_mesh());
    const NDIndex<DIM> gPDom = bunch.getGlobalPDomainInclOld(rho.get_mesh());
    const Index vertAdd = getExtraMarginVert();
    const Vector_t dx(rho.get_mesh().get_meshSpacing(0),
                      rho.get_mesh().get_meshSpacing(1));
    const Vector_t origin = rho.get_mesh().get_origin();

    FieldPatch<double> myRho(NDIndex<DIM>(Index(lPDom[0].first() + vertAdd.min(), lPDom[0].last() + vertAdd.max()),
                                          Index(lPDom[1].first() + vertAdd.min(), lPDom[1].last() + vertAdd.max())));
    myRho.setOrigin(origin + Vector_t((lPDom[0].first() + vertAdd.min()) * dx(0),
                                      (lPDom[1].first() + vertAdd.min()) * dx(1)));
    myRho.setSpacing(dx);

    const ParticleAttrib< Vector_t > & R = bunch.R;
    const ParticleAttrib< Vector_t > & oldR = bunch.oldR;
    const ParticleAttrib< double > & Q = bunch.Q;

    for (size_t i = 0; i < bunch.getLocalNP(); ++ i) {
        const double & q = Q[i];
        Vector_t xi((R[i](0) - origin(0)) / dx(0),
                    (R[i](1) - origin(1)) / dx(1));
        int I = static_cast<int>(floor(xi(0)));
        int J = static_cast<int>(floor(xi(1)));
        Vector_t tau(xi(0) - I, xi(1) - J);

        myRho(I, J)         += (1.0 - tau(0)) * (1.0 - tau(1)) * q;
        myRho(I + 1, J)     += tau(0)         * (1.0 - tau(1)) * q;
        myRho(I, J + 1)     += (1.0 - tau(0)) * tau(1)         * q;
        myRho(I + 1, J + 1) += tau(0)         * tau(1)         * q;

        xi = Vector_t((oldR[i](0) - origin(0)) / dx(0),
                      (oldR[i](1) - origin(1)) / dx(1));
        I = static_cast<int>(floor(xi(0)));
        J = static_cast<int>(floor(xi(1)));
        tau = Vector_t(xi(0) - I, xi(1) - J);

        myRho(I, J)         -= (1.0 - tau(0)) * (1.0 - tau(1)) * q;
        myRho(I + 1, J)     -= tau(0)         * (1.0 - tau(1)) * q;
        myRho(I, J + 1)     -= (1.0 - tau(0)) * tau(1)         * q;
        myRho(I + 1, J + 1) -= tau(0)         * tau(1)         * q;

    }

    const int lowerI = std::max(lDom[0].first(), gPDom[0].first() + vertAdd.min());
    const int upperI = std::min(lDom[0].last(),  gPDom[0].last() + vertAdd.max());
    const int signI = lowerI <= upperI? 1: upperI - lowerI;

    const int lowerJ = std::max(lDom[1].first(), gPDom[1].first() + vertAdd.min());
    const int upperJ = std::min(lDom[1].last(),  gPDom[1].last() + vertAdd.max());
    const int signJ = lowerJ <= upperJ? 1: upperJ - lowerJ;

    NDIndex<DIM> dom(Index(lowerI, upperI, signI),
                     Index(lowerJ, upperJ, signJ));

    Timings::stopTimer(chargeTimer);
    Communicator::communicateFields(myRho, dom, commTimer);
    Timings::startTimer(chargeTimer);

    rho = 0.0;
    for (int j = lowerJ; j <= upperJ; ++ j) {
        elem[1] = Index(j, j);
        for (int i = lowerI; i <= upperI; ++ i) {
            elem[0] = Index(i, i);
            rho.localElement(elem) += myRho(i,j) / (dx(0) * dx(1));
        }
    }
    rho.fillGuardCells();
}

void ZerothOrderShapeFunction::getNumberDensityImpl(SField_Cell_t& dens,
                                                    const PartBunch & bunch,
                                                    const Timings::TimerRef & commTimer)
{
    NDIndex<DIM> elem;
    const Mesh_t & mesh = dens.get_mesh();
    const Vector_t origin = mesh.get_origin();
    const Vector_t dx(mesh.get_meshSpacing(0),
                      mesh.get_meshSpacing(1));
    const NDIndex<DIM> lDom = dens.getLayout().getLocalNDIndex();
    const NDIndex<DIM> pDom = bunch.getLocalPDomain(mesh);
    const NDIndex<DIM> gPDom = bunch.getGlobalPDomain(mesh);
    const Index cellAdd = getExtraMarginCell();

    FieldPatch<double> myDens(NDIndex<DIM>(Index(pDom[0].first() + cellAdd.first(), pDom[0].last() + cellAdd.last()),
                                           Index(pDom[1].first() + cellAdd.first(), pDom[1].last() + cellAdd.last())));
    myDens.setOrigin(origin);
    myDens.setSpacing(dx);

    const ParticleAttrib< Vector_t > & R = bunch.R;

    for (size_t i = 0; i < bunch.getLocalNP(); ++ i) {
        const Vector_t xi((R[i](0) - origin(0)) / dx(0),
                          (R[i](1) - origin(1)) / dx(1));
        const int I = static_cast<int>(floor(xi(0)));
        const int J = static_cast<int>(floor(xi(1)));
        const Vector_t tau(xi(0) - I, xi(1) - J);
        myDens(I, J)         += (1.0 - tau(0)) * (1.0 - tau(1));
        myDens(I + 1, J)     += tau(0)         * (1.0 - tau(1));
        myDens(I, J + 1)     += (1.0 - tau(0)) * tau(1);
        myDens(I + 1, J + 1) += tau(0)         * tau(1);
    }

    const int lowerI = std::max(lDom[0].first(), gPDom[0].first() + cellAdd.min());
    const int upperI = std::min(lDom[0].last(),  gPDom[0].last() + cellAdd.max());
    const int signI = lowerI <= upperI? 1: upperI - lowerI;

    const int lowerJ = std::max(lDom[1].first(), gPDom[1].first() + cellAdd.min());
    const int upperJ = std::min(lDom[1].last(),  gPDom[1].last() + cellAdd.max());
    const int signJ = lowerJ <= upperJ? 1: upperJ - lowerJ;

    NDIndex<DIM> dom(Index(lowerI, upperI, signI),
                     Index(lowerJ, upperJ, signJ));

    Communicator::communicateFields(myDens, dom, commTimer);

    dens = 0.0;
    for (int j = lowerJ; j <= upperJ; ++ j) {
        elem[1] = Index(j, j);
        for (int i = lowerI; i <= upperI; ++ i) {
            elem[0] = Index(i, i);
            dens.localElement(elem) = myDens(i,j);
        }
    }
    dens.fillGuardCells();
}
