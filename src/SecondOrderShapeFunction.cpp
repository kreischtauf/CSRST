/***************************************************************************
                     SecondOrderShapeFunction.cpp
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

#include "SecondOrderShapeFunction.hh"
#include "PartBunch.hh"
#include "FieldPatch.hh"
#include "Communicator.hh"
#include "utils.hh"

extern std::ofstream dbg;
#define DBGOUT dbg << __FILE__ << ": " << __LINE__ << "\t"

void SecondOrderShapeFunction::getCurrentDensityImplX(VField_Edge_t& JFD,
                                                      const FieldPatch<double> & Jx,
                                                      const PartBunch & bunch,
                                                      const double & dt,
                                                      const Timings::TimerRef & currentTimer,
                                                      const Timings::TimerRef & commTimer,
                                                      const size_t & from,
                                                      const size_t & to)
{
    NDIndex<DIM> elem;
    const NDIndex<DIM> currentAdd = getExtraMarginCurrent();
    const unsigned int myNode = Ippl::myNode();
    const unsigned int numNodes = Ippl::getNodes();

    Timings::stopTimer(currentTimer);
    const std::vector<size_t> & numLocalParticles = bunch.getNumLocalParticles();

    NDIndex<DIM> lPDom = bunch.getLocalPDomainInclOld(JFD.get_mesh());
    FieldPatch<double> myJx = getCurrentDensityBaseImplX(JFD, bunch, dt, currentTimer, from, to);
    if (Jx.size() > 0) {
        myJx.add(Jx);
    }

    auto localPDomains = Communicator::getAllLocalDomains(lPDom);

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

            dom[d] = Index(dom[d].first() + currentAdd[d].min(),
                           dom[d].last() + currentAdd[d].max());
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

void SecondOrderShapeFunction::getCurrentDensityImplY(VField_Edge_t& JFD,
                                                      const FieldPatch<double> & Jy,
                                                      const PartBunch & bunch,
                                                      const double & dt,
                                                      const Timings::TimerRef & currentTimer,
                                                      const Timings::TimerRef & commTimer,
                                                      const size_t & from,
                                                      const size_t & to)
{
    NDIndex<DIM> elem;
    const NDIndex<DIM> currentAdd = getExtraMarginCurrent();
    const unsigned int myNode = Ippl::myNode();
    const unsigned int numNodes = Ippl::getNodes();

    Timings::stopTimer(currentTimer);
    const std::vector<size_t> & numLocalParticles = bunch.getNumLocalParticles();

    NDIndex<DIM> lPDom = bunch.getLocalPDomainInclOld(JFD.get_mesh());
    FieldPatch<double> myJy = getCurrentDensityBaseImplY(JFD, bunch, dt, currentTimer, from, to);
    if (Jy.size() > 0) {
        myJy.add(Jy);
    }

    auto localPDomains = Communicator::getAllLocalDomains(lPDom);

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

            dom[d] = Index(dom[d].first() + currentAdd[d].min(),
                           dom[d].last() + currentAdd[d].max());
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

FieldPatch<double> SecondOrderShapeFunction::getCurrentDensityBaseImplX(VField_Edge_t& JFD,
                                                                        const PartBunch & bunch,
                                                                        const double & dt,
                                                                        const Timings::TimerRef & currentTimer,
                                                                        const size_t & from,
                                                                        const size_t & to)
{
    const NDIndex<DIM> currentAdd = getExtraMarginCurrent();

    Timings::stopTimer(currentTimer);
    NDIndex<DIM> lPDom = bunch.getLocalPDomainInclOld(JFD.get_mesh());
    NDIndex<DIM> dom;
    for (unsigned int d = 0; d < DIM; ++ d) {
        dom[d] = Index(lPDom[d].first() + currentAdd[d].first(),
                         lPDom[d].last() + currentAdd[d].last());
    }

    Timings::startTimer(currentTimer);

    const Vector_t dx(JFD.get_mesh().get_meshSpacing(0),
                      JFD.get_mesh().get_meshSpacing(1));
    const Vector_t origin = JFD.get_mesh().get_origin();
    Utils::samePos uni(sqrt(dot(dx,dx)));
    const Vector_t oneoverdx = Vector_t(1.0) / dx;
    const Vector_t oneoverdxdt = oneoverdx / dt;
    const double onehalf = 0.5;

    FieldPatch<double> myJx(dom);
    myJx.setOrigin(origin + Vector_t(dom[0].first() * dx(0),
                                     dom[1].first() * dx(1)));
    myJx.setSpacing(dx);

    const ParticleAttrib< Vector_t > & R = bunch.R;
    const ParticleAttrib< Vector_t > & oldR = bunch.oldR;
    const ParticleAttrib< double > & Q = bunch.Q;

    double current[12];
    std::fill(current, current + 12, 0.0);

    Utils::CompanionSorterItem xy_crossings[10];
    for (size_t i = from; i < to; ++ i) {
        const double & q = Q[i];

        size_t numXs = Utils::calcCellBorderCrossings(xy_crossings, oldR[i], R[i], dx, origin, uni);

        const Vector_t Dx = R[i] - oldR[i];
        const Vector_t rDx = Dx * oneoverdx;
        const double coveredDistance = sqrt(dot(Dx,Dx));

        for (size_t j = 0; j < numXs - 1; ++ j) {
            const Vector_t & xi = xy_crossings[j].subStep * oneoverdx;
            const Vector_t & xf = xy_crossings[j+1].subStep * oneoverdx;
            const Vector_t center = onehalf * (xi + xf);
            const int I = static_cast<long>(floor(center(0)));
            const int J = static_cast<long>(floor(center(1)));

            const Vector_t segment = xy_crossings[j+1].subStep - xy_crossings[j].subStep;
            const double tau = sqrt(dot(segment, segment)) / coveredDistance;
            const Vector_t a(I + onehalf - xi(0), J + onehalf - xi(1));
            const Vector_t b(I - xi(0), J - xi(1));

            getCurrentX(current, tau, a, b, q * oneoverdxdt, rDx);

            myJx(I - 1, J - 1) += current[0];
            myJx(I,     J - 1) += current[1];
            myJx(I + 1, J - 1) += current[2];
            myJx(I - 1, J    ) += current[3];
            myJx(I,     J    ) += current[4];
            myJx(I + 1, J    ) += current[5];
            myJx(I - 1, J + 1) += current[6];
            myJx(I,     J + 1) += current[7];
            myJx(I + 1, J + 1) += current[8];
            myJx(I - 1, J + 2) += current[9];
            myJx(I,     J + 2) += current[10];
            myJx(I + 1, J + 2) += current[11];

        }
    }

    Timings::stopTimer(currentTimer);

    return myJx;
}

FieldPatch<double> SecondOrderShapeFunction::getCurrentDensityBaseImplY(VField_Edge_t& JFD,
                                                                        const PartBunch & bunch,
                                                                        const double & dt,
                                                                        const Timings::TimerRef & currentTimer,
                                                                        const size_t & from,
                                                                        const size_t & to)
{
    const NDIndex<DIM> currentAdd = getExtraMarginCurrent();

    Timings::stopTimer(currentTimer);
    NDIndex<DIM> lPDom = bunch.getLocalPDomainInclOld(JFD.get_mesh());
    NDIndex<DIM> dom;
    for (unsigned int d = 0; d < DIM; ++ d) {
        dom[d] = Index(lPDom[d].first() + currentAdd[d].first(),
                         lPDom[d].last() + currentAdd[d].last());
    }

    Timings::startTimer(currentTimer);

    const Vector_t dx(JFD.get_mesh().get_meshSpacing(0),
                      JFD.get_mesh().get_meshSpacing(1));
    const Vector_t origin = JFD.get_mesh().get_origin();
    const Vector_t oneoverdx = Vector_t(1.0) / dx;
    const Vector_t oneoverdxdt = oneoverdx / dt;
    Utils::samePos uni(sqrt(dot(dx,dx)));
    const double onehalf = 0.5;

    FieldPatch<double> myJy(dom);
    myJy.setOrigin(origin + Vector_t(dom[0].first() * dx(0),
                                     dom[1].first() * dx(1)));
    myJy.setSpacing(dx);

    const ParticleAttrib< Vector_t > & R = bunch.R;
    const ParticleAttrib< Vector_t > & oldR = bunch.oldR;
    const ParticleAttrib< double > & Q = bunch.Q;

    double current[12];
    std::fill(current, current + 12, 0.0);

    Utils::CompanionSorterItem xy_crossings[10];
    for (size_t i = from; i < to; ++ i) {
        const double & q = Q[i];

        size_t numXs = Utils::calcCellBorderCrossings(xy_crossings, oldR[i], R[i], dx, origin, uni);

        const Vector_t Dx = R[i] - oldR[i];
        const Vector_t rDx = Dx * oneoverdx;
        const double coveredDistance = sqrt(dot(Dx,Dx));

        for (size_t j = 0; j < numXs - 1; ++ j) {
            const Vector_t & xi = xy_crossings[j].subStep * oneoverdx;
            const Vector_t & xf = xy_crossings[j+1].subStep * oneoverdx;
            const Vector_t center = onehalf * (xi + xf);
            const int I = static_cast<long>(floor(center(0)));
            const int J = static_cast<long>(floor(center(1)));

            const Vector_t segment = xy_crossings[j+1].subStep - xy_crossings[j].subStep;
            const double tau = sqrt(dot(segment, segment)) / coveredDistance;
            const Vector_t a(I + onehalf - xi(0), J + onehalf - xi(1));
            const Vector_t b(I - xi(0), J - xi(1));

            getCurrentY(current, tau, a, b, q * oneoverdxdt, rDx);

            myJy(I - 1, J - 1) += current[0];
            myJy(I,     J - 1) += current[1];
            myJy(I + 1, J - 1) += current[2];
            myJy(I + 2, J - 1) += current[3];
            myJy(I - 1, J    ) += current[4];
            myJy(I,     J    ) += current[5];
            myJy(I + 1, J    ) += current[6];
            myJy(I + 2, J    ) += current[7];
            myJy(I - 1, J + 1) += current[8];
            myJy(I,     J + 1) += current[9];
            myJy(I + 1, J + 1) += current[10];
            myJy(I + 2, J + 1) += current[11];

        }
    }

    Timings::stopTimer(currentTimer);

    return myJy;
}

void SecondOrderShapeFunction::getChargeDensityImpl(SField_Vert_t& rho,
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

    const double onesixth = 1.0 / 6.0;
    for (size_t i = 0; i < bunch.getLocalNP(); ++ i) {
        const double & q = Q[i];
        const Vector_t neR = (R[i] - origin) / dx;
        const int I = static_cast<int>(floor(neR(0)));
        const int J = static_cast<int>(floor(neR(1)));

        const double a = neR(0) - I;
        const double b = neR(1) - J;

        const double sxm1 = onesixth * powthree(1.0 - a);
        const double sx = onesixth * (4.0 - 3.0 * a*a * (2.0 - a));
        const double sxp1 = onesixth * (4.0 - 3.0 * (1.0 - a)*(1.0 - a) * (1.0 + a));
        const double sxp2 = onesixth * powthree(a);
        const double sym1 = onesixth * powthree(1.0 - b);
        const double sy = onesixth * (4.0 - 3.0 * b*b * (2.0 - b));
        const double syp1 = onesixth * (4.0 - 3.0 * (1.0 - b)*(1.0 - b) * (1.0 + b));
        const double syp2 = onesixth * powthree(b);

        myRho(I-1, J-1) += q * sxm1 * sym1;
        myRho(I,   J-1) += q * sx   * sym1;
        myRho(I+1, J-1) += q * sxp1 * sym1;
        myRho(I+2, J-1) += q * sxp2 * sym1;
        myRho(I-1, J  ) += q * sxm1 * sy;
        myRho(I,   J  ) += q * sx   * sy;
        myRho(I+1, J  ) += q * sxp1 * sy;
        myRho(I+2, J  ) += q * sxp2 * sy;
        myRho(I-1, J+1) += q * sxm1 * syp1;
        myRho(I,   J+1) += q * sx   * syp1;
        myRho(I+1, J+1) += q * sxp1 * syp1;
        myRho(I+2, J+1) += q * sxp2 * syp1;
        myRho(I-1, J+2) += q * sxm1 * syp2;
        myRho(I,   J+2) += q * sx   * syp2;
        myRho(I+1, J+2) += q * sxp1 * syp2;
        myRho(I+2, J+2) += q * sxp2 * syp2;
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

void SecondOrderShapeFunction::getChargeDensityDiffImpl(SField_Vert_t& rho,
                                                        const PartBunch & bunch,
                                                        const Timings::TimerRef & chargeTimer,
                                                        const Timings::TimerRef & commTimer)
{
    const double onesixth = 1.0 / 6.0;
    NDIndex<DIM> elem;
    const NDIndex<DIM> lDom = rho.getLayout().getLocalNDIndex();
    const NDIndex<DIM> & lPDom = bunch.getLocalPDomainInclOld(rho.get_mesh());
    const NDIndex<DIM> & gPDom = bunch.getGlobalPDomainInclOld(rho.get_mesh());
    const Index vertAdd = getExtraMarginVert();
    const Vector_t dx(rho.get_mesh().get_meshSpacing(0),
                      rho.get_mesh().get_meshSpacing(1));
    const Vector_t origin = rho.get_mesh().get_origin();

    const int minI = lPDom[0].first(), maxI = lPDom[0].last();
    const int minJ = lPDom[1].first(), maxJ = lPDom[1].last();

    FieldPatch<double> myRho(NDIndex<DIM>(Index(minI + vertAdd.min(), maxI + vertAdd.max()),
                                  Index(minJ + vertAdd.min(), maxJ + vertAdd.max())));
    myRho.setOrigin(origin + Vector_t((minI + vertAdd.min()) * dx(0),
                                      (minJ + vertAdd.min()) * dx(1)));
    myRho.setSpacing(dx);

    const ParticleAttrib< Vector_t > & R = bunch.R;
    const ParticleAttrib< Vector_t > & oldR = bunch.oldR;
    const ParticleAttrib< double > & Q = bunch.Q;

    for (size_t i = 0; i < bunch.getLocalNP(); ++ i) {
        const double & q = Q[i];
        const Vector_t neR = (R[i] - origin) / dx;
        const Vector_t olR = (oldR[i] - origin) / dx;
        int I = static_cast<int>(floor(neR(0)));
        int J = static_cast<int>(floor(neR(1)));

        double a = neR(0) - I;
        double b = neR(1) - J;

        double sxm1 = onesixth * powthree(1.0 - a);
        double sx = onesixth * (4.0 - 3.0 * a*a * (2.0 - a));
        double sxp1 = onesixth * (4.0 - 3.0 * (1.0 - a)*(1.0 - a) * (1.0 + a));
        double sxp2 = onesixth * powthree(a);
        double sym1 = onesixth * powthree(1.0 - b);
        double sy = onesixth * (4.0 - 3.0 * b*b * (2.0 - b));
        double syp1 = onesixth * (4.0 - 3.0 * (1.0 - b)*(1.0 - b) * (1.0 + b));
        double syp2 = onesixth * powthree(b);

        myRho(I-1, J-1) += q * sxm1 * sym1;
        myRho(I,   J-1) += q * sx   * sym1;
        myRho(I+1, J-1) += q * sxp1 * sym1;
        myRho(I+2, J-1) += q * sxp2 * sym1;
        myRho(I-1, J  ) += q * sxm1 * sy;
        myRho(I,   J  ) += q * sx   * sy;
        myRho(I+1, J  ) += q * sxp1 * sy;
        myRho(I+2, J  ) += q * sxp2 * sy;
        myRho(I-1, J+1) += q * sxm1 * syp1;
        myRho(I,   J+1) += q * sx   * syp1;
        myRho(I+1, J+1) += q * sxp1 * syp1;
        myRho(I+2, J+1) += q * sxp2 * syp1;
        myRho(I-1, J+2) += q * sxm1 * syp2;
        myRho(I,   J+2) += q * sx   * syp2;
        myRho(I+1, J+2) += q * sxp1 * syp2;
        myRho(I+2, J+2) += q * sxp2 * syp2;

        I = static_cast<int>(floor(olR(0)));
        J = static_cast<int>(floor(olR(1)));

        a = olR(0) - I;
        b = olR(1) - J;

        sxm1 = onesixth * powthree(1.0 - a);
        sx = onesixth * (4.0 - 3.0 * a*a * (2.0 - a));
        sxp1 = onesixth * (4.0 - 3.0 * (1.0 - a)*(1.0 - a) * (1.0 + a));
        sxp2 = onesixth * powthree(a);
        sym1 = onesixth * powthree(1.0 - b);
        sy = onesixth * (4.0 - 3.0 * b*b * (2.0 - b));
        syp1 = onesixth * (4.0 - 3.0 * (1.0 - b)*(1.0 - b) * (1.0 + b));
        syp2 = onesixth * powthree(b);

        myRho(I-1, J-1) -= q * sxm1 * sym1;
        myRho(I,   J-1) -= q * sx   * sym1;
        myRho(I+1, J-1) -= q * sxp1 * sym1;
        myRho(I+2, J-1) -= q * sxp2 * sym1;
        myRho(I-1, J  ) -= q * sxm1 * sy;
        myRho(I,   J  ) -= q * sx   * sy;
        myRho(I+1, J  ) -= q * sxp1 * sy;
        myRho(I+2, J  ) -= q * sxp2 * sy;
        myRho(I-1, J+1) -= q * sxm1 * syp1;
        myRho(I,   J+1) -= q * sx   * syp1;
        myRho(I+1, J+1) -= q * sxp1 * syp1;
        myRho(I+2, J+1) -= q * sxp2 * syp1;
        myRho(I-1, J+2) -= q * sxm1 * syp2;
        myRho(I,   J+2) -= q * sx   * syp2;
        myRho(I+1, J+2) -= q * sxp1 * syp2;
        myRho(I+2, J+2) -= q * sxp2 * syp2;
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

void SecondOrderShapeFunction::getCurrentX(double J[],
                                           const double & tau,
                                           const Vector_t & a,
                                           const Vector_t & b,
                                           const Vector_t & qoverdxdt,
                                           const Vector_t & delta)
{
    const Vector_t current_density(qoverdxdt(1) * delta(0),
                                   qoverdxdt(0) * delta(1));

    const Vector_t deltapow2 = delta*delta;
    const Vector_t deltapow3 = deltapow2*delta;
    const double taupow2 = tau*tau / 2;
    const double taupow3 = taupow2*tau * 2 / 3;
    const double taupow4 = taupow2*taupow2;

    const double half = 0.5;
    const double onethird = 1.0 / 3.0;
    const double onefourth = 0.25;
    const double onefifth = 0.2;
    const double onesixth = half * onethird;

    const double Ax2Im1 = half * (half + a(0))*(half + a(0));
    const double Bx2Im1 = -(half + a(0)) * delta(0) * tau;
    const double Cx2Im1 = half * deltapow2(0) * taupow2;
    const double Ax2I = 0.75 - a(0)*a(0);
    const double Bx2I = 2 * a(0) * delta(0) * tau;
    const double Cx2I = -deltapow2(0) * taupow2;
    const double Ax2Ip1 = half * (0.5 - a(0))*(0.5 - a(0));
    const double Bx2Ip1 = (half - a(0)) * delta(0) * tau;
    const double Cx2Ip1 = half * deltapow2(0) * taupow2;

    const double Ay3Jm1 = (1 + b(1))*(1 + b(1))*(1 + b(1)) * onesixth * tau;
    const double By3Jm1 = -(1 + b(1))*(1 + b(1)) * delta(1) * half * taupow2;
    const double Cy3Jm1 = (1 + b(1)) * deltapow2(1) * half * taupow3;
    const double Dy3Jm1 = -deltapow3(1) * onesixth * taupow4;
    const double Ay3J = (4 - 3 * b(1)*b(1) * (2 + b(1))) * onesixth * tau;
    const double By3J = b(1) * delta(1) * (3 * b(1) + 4) * half * taupow2;
    const double Cy3J = -deltapow2(1) * (2 + 3 * b(1)) * half * taupow3;
    const double Dy3J = deltapow3(1) * half * taupow4;
    const double Ay3Jp1 = (4 - 3 * (b(1) + 1)*(b(1) + 1) * (1 - b(1))) * onesixth * tau;
    const double By3Jp1 = (1 - 2 * b(1) - 3 * b(1)*b(1)) * delta(1) * half * taupow2;
    const double Cy3Jp1 = (1 + 3 * b(1)) * deltapow2(1) * half * taupow3;
    const double Dy3Jp1 = -deltapow3(1) * half * taupow4;
    const double Ay3Jp2 = -b(1)*b(1)*b(1) * onesixth * tau;
    const double By3Jp2 = b(1)*b(1) * delta(1) * half * taupow2;
    const double Cy3Jp2 = -b(1) * deltapow2(1) * half * taupow3;
    const double Dy3Jp2 = deltapow3(1) * onesixth * taupow4;

    // Jx: I-1, J-1
    J[0] = current_density(0) *  ((Ax2Im1 + Bx2Im1 + Cx2Im1) * (Dy3Jm1 + Cy3Jm1 + By3Jm1 + Ay3Jm1)
                                  + (Cx2Im1 * (Dy3Jm1 - Ay3Jm1) - Bx2Im1 * By3Jm1) * onethird
                                  + (Cx2Im1 * Cy3Jm1 - Bx2Im1 * Dy3Jm1) * onefifth
                                  - Bx2Im1 * Cy3Jm1 * onefourth - Bx2Im1 * Ay3Jm1 * half);

    // Jx: I, J-1
    J[1] = current_density(0) *  ((Ax2I + Bx2I + Cx2I) * (Dy3Jm1 + Cy3Jm1 + By3Jm1 + Ay3Jm1)
                                  + (Cx2I * (Dy3Jm1 - Ay3Jm1) - Bx2I * By3Jm1) * onethird
                                  +  (Cx2I * Cy3Jm1 - Bx2I * Dy3Jm1) * onefifth
                                  -  Bx2I * Cy3Jm1 * onefourth - Bx2I * Ay3Jm1 * half);

    // Jx: I+1, J-1
    J[2] = current_density(0) *  ((Ax2Ip1 + Bx2Ip1 + Cx2Ip1) * (Dy3Jm1 + Cy3Jm1 + By3Jm1 + Ay3Jm1)
                                  + (Cx2Ip1 * (Dy3Jm1 - Ay3Jm1) - Bx2Ip1 * By3Jm1) * onethird
                                  + (Cx2Ip1 * Cy3Jm1 - Bx2Ip1 * Dy3Jm1) * onefifth
                                  - Bx2Ip1 * Cy3Jm1 * onefourth - Bx2Ip1 * Ay3Jm1 * half);

    // Jx: I-1, J
    J[3] = current_density(0) *  ((Ax2Im1 + Bx2Im1 + Cx2Im1) * (Dy3J + Cy3J + By3J + Ay3J)
                                  + (Cx2Im1 * (Dy3J - Ay3J) - Bx2Im1 * By3J) * onethird
                                  + (Cx2Im1 * Cy3J - Bx2Im1 * Dy3J) * onefifth
                                  - Bx2Im1 * Cy3J * onefourth - Bx2Im1 * Ay3J * half);

    // Jx: I, J
    J[4] = current_density(0) *  ((Ax2I + Bx2I + Cx2I) * (Dy3J + Cy3J + By3J + Ay3J)
                                  + (Cx2I * (Dy3J - Ay3J) - Bx2I * By3J) * onethird
                                  +  (Cx2I * Cy3J - Bx2I * Dy3J) * onefifth
                                  -  Bx2I * Cy3J * onefourth - Bx2I * Ay3J * half);

    // Jx: I+1, J
    J[5] = current_density(0) *  ((Ax2Ip1 + Bx2Ip1 + Cx2Ip1) * (Dy3J + Cy3J + By3J + Ay3J)
                                  + (Cx2Ip1 * (Dy3J - Ay3J) - Bx2Ip1 * By3J) * onethird
                                  + (Cx2Ip1 * Cy3J - Bx2Ip1 * Dy3J) * onefifth
                                  - Bx2Ip1 * Cy3J * onefourth - Bx2Ip1 * Ay3J * half);

    // Jx: I-1, J+1
    J[6] = current_density(0) *  ((Ax2Im1 + Bx2Im1 + Cx2Im1) * (Dy3Jp1 + Cy3Jp1 + By3Jp1 + Ay3Jp1)
                                  + (Cx2Im1 * (Dy3Jp1 - Ay3Jp1) - Bx2Im1 * By3Jp1) * onethird
                                  + (Cx2Im1 * Cy3Jp1 - Bx2Im1 * Dy3Jp1) * onefifth
                                  - Bx2Im1 * Cy3Jp1 * onefourth - Bx2Im1 * Ay3Jp1 * half);

    // Jx: I, J+1
    J[7] = current_density(0) *  ((Ax2I + Bx2I + Cx2I) * (Dy3Jp1 + Cy3Jp1 + By3Jp1 + Ay3Jp1)
                                  + (Cx2I * (Dy3Jp1 - Ay3Jp1) - Bx2I * By3Jp1) * onethird
                                  +  (Cx2I * Cy3Jp1 - Bx2I * Dy3Jp1) * onefifth
                                  -  Bx2I * Cy3Jp1 * onefourth - Bx2I * Ay3Jp1 * half);

    // Jx: I+1, J+1
    J[8] = current_density(0) *  ((Ax2Ip1 + Bx2Ip1 + Cx2Ip1) * (Dy3Jp1 + Cy3Jp1 + By3Jp1 + Ay3Jp1)
                                  + (Cx2Ip1 * (Dy3Jp1 - Ay3Jp1) - Bx2Ip1 * By3Jp1) * onethird
                                  + (Cx2Ip1 * Cy3Jp1 - Bx2Ip1 * Dy3Jp1) * onefifth
                                  - Bx2Ip1 * Cy3Jp1 * onefourth - Bx2Ip1 * Ay3Jp1 * half);

    // Jx: I-1, J+2
    J[9] = current_density(0) *  ((Ax2Im1 + Bx2Im1 + Cx2Im1) * (Dy3Jp2 + Cy3Jp2 + By3Jp2 + Ay3Jp2)
                                  + (Cx2Im1 * (Dy3Jp2 - Ay3Jp2) - Bx2Im1 * By3Jp2) * onethird
                                  + (Cx2Im1 * Cy3Jp2 - Bx2Im1 * Dy3Jp2) * onefifth
                                  - Bx2Im1 * Cy3Jp2 * onefourth - Bx2Im1 * Ay3Jp2 * half);

    // Jx;  I, J+2
    J[10] = current_density(0) * ((Ax2I + Bx2I + Cx2I) * (Dy3Jp2 + Cy3Jp2 + By3Jp2 + Ay3Jp2)
                                  + (Cx2I * (Dy3Jp2 - Ay3Jp2) - Bx2I * By3Jp2) * onethird
                                  +  (Cx2I * Cy3Jp2 - Bx2I * Dy3Jp2) * onefifth
                                  -  Bx2I * Cy3Jp2 * onefourth - Bx2I * Ay3Jp2 * half);

    // Jx: I+1, J+2
    J[11] = current_density(0) * ((Ax2Ip1 + Bx2Ip1 + Cx2Ip1) * (Dy3Jp2 + Cy3Jp2 + By3Jp2 + Ay3Jp2)
                                  + (Cx2Ip1 * (Dy3Jp2 - Ay3Jp2) - Bx2Ip1 * By3Jp2) * onethird
                                  + (Cx2Ip1 * Cy3Jp2 - Bx2Ip1 * Dy3Jp2) * onefifth
                                  - Bx2Ip1 * Cy3Jp2 * onefourth - Bx2Ip1 * Ay3Jp2 * half);

}

void SecondOrderShapeFunction::getCurrentY(double J[],
                                           const double & tau,
                                           const Vector_t & a,
                                           const Vector_t & b,
                                           const Vector_t & qoverdxdt,
                                           const Vector_t & delta)
{
    const Vector_t current_density(qoverdxdt(1) * delta(0),
                                   qoverdxdt(0) * delta(1));

    const Vector_t deltapow2 = delta*delta;
    const Vector_t deltapow3 = deltapow2*delta;
    const double taupow2 = tau*tau / 2;
    const double taupow3 = taupow2*tau * 2 / 3;
    const double taupow4 = taupow2*taupow2;

    const double half = 0.5;
    const double onethird = 1.0 / 3.0;
    const double onefourth = 0.25;
    const double onefifth = 0.2;
    const double onesixth = half * onethird;

    const double Ax3Im1 = (1 + b(0))*(1 + b(0))*(1 + b(0)) * onesixth * tau;
    const double Bx3Im1 = -(1 + b(0))*(1 + b(0)) * delta(0) * half * taupow2;
    const double Cx3Im1 = (1 + b(0)) * deltapow2(0) * half * taupow3;
    const double Dx3Im1 = -deltapow3(0) * onesixth * taupow4;
    const double Ax3I = (4 - 3 * b(0)*b(0) * (2 + b(0))) * onesixth * tau;
    const double Bx3I = b(0) * delta(0) * (3 * b(0) + 4) * half * taupow2;
    const double Cx3I = -deltapow2(0) * (2 + 3 * b(0)) * half * taupow3;
    const double Dx3I = deltapow3(0) * half * taupow4;
    const double Ax3Ip1 = (4 - 3 * (b(0) + 1)*(b(0) + 1) * (1 - b(0))) * onesixth * tau;
    const double Bx3Ip1 = (1 - 2 * b(0) - 3 * b(0)*b(0)) * delta(0) * half * taupow2;
    const double Cx3Ip1 = (1 + 3 * b(0)) * deltapow2(0) * half * taupow3;
    const double Dx3Ip1 = -deltapow3(0) * half * taupow4;
    const double Ax3Ip2 = -b(0)*b(0)*b(0) * onesixth * tau;
    const double Bx3Ip2 = b(0)*b(0) * delta(0) * half * taupow2;
    const double Cx3Ip2 = -b(0) * deltapow2(0) * half * taupow3;
    const double Dx3Ip2 = deltapow3(0) * onesixth * taupow4;

    const double Ay2Jm1 = (half + a(1))*(half + a(1)) * half;
    const double By2Jm1 = -(half + a(1)) * delta(1) * tau;
    const double Cy2Jm1 = deltapow2(1) * half * taupow2;
    const double Ay2J = 0.75 - a(1)*a(1);
    const double By2J = 2 * a(1) * delta(1) * tau;
    const double Cy2J = -deltapow2(1) * taupow2;
    const double Ay2Jp1 = (half - a(1))*(half - a(1)) * half;
    const double By2Jp1 = (half - a(1)) * delta(1) * tau;
    const double Cy2Jp1 = deltapow2(1) * half * taupow2;

    // Jy: I-1, J-1
    J[ 0] = current_density(1) * ((Ay2Jm1 + By2Jm1 + Cy2Jm1) * (Dx3Im1 + Cx3Im1 + Bx3Im1 + Ax3Im1)
                                  + (Cy2Jm1 * (Dx3Im1 - Ax3Im1) - By2Jm1 * Bx3Im1) * onethird
                                  + (Cy2Jm1 * Cx3Im1 - By2Jm1 * Dx3Im1) * onefifth
                                  - By2Jm1 * Cx3Im1 * onefourth - By2Jm1 * Ax3Im1 * half);

    // Jy: I, J-1
    J[ 1] = current_density(1) * ((Ay2Jm1 + By2Jm1 + Cy2Jm1) * (Dx3I + Cx3I + Bx3I + Ax3I)
                                  + (Cy2Jm1 * (Dx3I - Ax3I) - By2Jm1 * Bx3I) * onethird
                                  + (Cy2Jm1 * Cx3I - By2Jm1 * Dx3I) * onefifth
                                  - By2Jm1 * Cx3I * onefourth - By2Jm1 * Ax3I * half);

    // Jy; I+1, J-1
    J[ 2] = current_density(1) * ((Ay2Jm1 + By2Jm1 + Cy2Jm1) * (Dx3Ip1 + Cx3Ip1 + Bx3Ip1 + Ax3Ip1)
                                  + (Cy2Jm1 * (Dx3Ip1 - Ax3Ip1) - By2Jm1 * Bx3Ip1) * onethird
                                  + (Cy2Jm1 * Cx3Ip1 - By2Jm1 * Dx3Ip1) * onefifth
                                  - By2Jm1 * Cx3Ip1 * onefourth - By2Jm1 * Ax3Ip1 * half);

    // Jy: I+2, J-1
    J[ 3] = current_density(1) * ((Ay2Jm1 + By2Jm1 + Cy2Jm1) * (Dx3Ip2 + Cx3Ip2 + Bx3Ip2 + Ax3Ip2)
                                  + (Cy2Jm1 * (Dx3Ip2 - Ax3Ip2) - By2Jm1 * Bx3Ip2) * onethird
                                  + (Cy2Jm1 * Cx3Ip2 - By2Jm1 * Dx3Ip2) * onefifth
                                  - By2Jm1 * Cx3Ip2 * onefourth - By2Jm1 * Ax3Ip2 * half);

    // Jy: I-1, J
    J[ 4] = current_density(1) * ((Ay2J + By2J + Cy2J) * (Dx3Im1 + Cx3Im1 + Bx3Im1 + Ax3Im1)
                                  + (Cy2J * (Dx3Im1 - Ax3Im1) - By2J * Bx3Im1) * onethird
                                  + (Cy2J * Cx3Im1 - By2J * Dx3Im1) * onefifth
                                  - By2J * Cx3Im1 * onefourth - By2J * Ax3Im1 * half);

    // Jy: I, J
    J[ 5] = current_density(1) * ((Ay2J + By2J + Cy2J) * (Dx3I + Cx3I + Bx3I + Ax3I)
                                  + (Cy2J * (Dx3I - Ax3I) - By2J * Bx3I) * onethird
                                  + (Cy2J * Cx3I - By2J * Dx3I) * onefifth
                                  - By2J * Cx3I * onefourth - By2J * Ax3I * half);

    // Jy: I+1, J
    J[ 6] = current_density(1) * ((Ay2J + By2J + Cy2J) * (Dx3Ip1 + Cx3Ip1 + Bx3Ip1 + Ax3Ip1)
                                  + (Cy2J * (Dx3Ip1 - Ax3Ip1) - By2J * Bx3Ip1) * onethird
                                  + (Cy2J * Cx3Ip1 - By2J * Dx3Ip1) * onefifth
                                  - By2J * Cx3Ip1 * onefourth - By2J * Ax3Ip1 * half);

    // Jy: I+2, J
    J[ 7] = current_density(1) * ((Ay2J + By2J + Cy2J) * (Dx3Ip2 + Cx3Ip2 + Bx3Ip2 + Ax3Ip2)
                                  + (Cy2J * (Dx3Ip2 - Ax3Ip2) - By2J * Bx3Ip2) * onethird
                                  + (Cy2J * Cx3Ip2 - By2J * Dx3Ip2) * onefifth
                                  - By2J * Cx3Ip2 * onefourth - By2J * Ax3Ip2 * half);

    // Jy: I-1, J+1
    J[ 8] = current_density(1) * ((Ay2Jp1 + By2Jp1 + Cy2Jp1) * (Dx3Im1 + Cx3Im1 + Bx3Im1 + Ax3Im1)
                                  + (Cy2Jp1 * (Dx3Im1 - Ax3Im1) - By2Jp1 * Bx3Im1) * onethird
                                  + (Cy2Jp1 * Cx3Im1 - By2Jp1 * Dx3Im1) * onefifth
                                  - By2Jp1 * Cx3Im1 * onefourth - By2Jp1 * Ax3Im1 * half);

    // Jy: I, J+1
    J[ 9] = current_density(1) * ((Ay2Jp1 + By2Jp1 + Cy2Jp1) * (Dx3I + Cx3I + Bx3I + Ax3I)
                                  + (Cy2Jp1 * (Dx3I - Ax3I) - By2Jp1 * Bx3I) * onethird
                                  + (Cy2Jp1 * Cx3I - By2Jp1 * Dx3I) * onefifth
                                  - By2Jp1 * Cx3I * onefourth - By2Jp1 * Ax3I * half);

    // Jy: I+1, J+1
    J[10] = current_density(1) * ((Ay2Jp1 + By2Jp1 + Cy2Jp1) * (Dx3Ip1 + Cx3Ip1 + Bx3Ip1 + Ax3Ip1)
                                  + (Cy2Jp1 * (Dx3Ip1 - Ax3Ip1) - By2Jp1 * Bx3Ip1) * onethird
                                  + (Cy2Jp1 * Cx3Ip1 - By2Jp1 * Dx3Ip1) * onefifth
                                  - By2Jp1 * Cx3Ip1 * onefourth - By2Jp1 * Ax3Ip1 * half);

    // Jy: I+2, J+1
    J[11] = current_density(1) * ((Ay2Jp1 + By2Jp1 + Cy2Jp1) * (Dx3Ip2 + Cx3Ip2 + Bx3Ip2 + Ax3Ip2)
                                  + (Cy2Jp1 * (Dx3Ip2 - Ax3Ip2) - By2Jp1 * Bx3Ip2) * onethird
                                  + (Cy2Jp1 * Cx3Ip2 - By2Jp1 * Dx3Ip2) * onefifth
                                  - By2Jp1 * Cx3Ip2 * onefourth - By2Jp1 * Ax3Ip2 * half);
}
