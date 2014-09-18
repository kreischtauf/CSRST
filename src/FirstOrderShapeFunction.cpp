/***************************************************************************
                     FirstOrderShapeFunction.cpp
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

#include "FirstOrderShapeFunction.hh"
#include "PartBunch.hh"
#include "FieldPatch.hh"
#include "Communicator.hh"
#include "utils.hh"

extern std::ofstream dbg;
#define DBGOUT dbg << __FILE__ << ": " << __LINE__ << "\t"

void FirstOrderShapeFunction::getCurrentDensityImplX(VField_Edge_t& JFD,
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

void FirstOrderShapeFunction::getCurrentDensityImplY(VField_Edge_t& JFD,
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

FieldPatch<double> FirstOrderShapeFunction::getCurrentDensityBaseImplX(VField_Edge_t& JFD,
                                                                       const PartBunch & bunch,
                                                                       const double & dt,
                                                                       const Timings::TimerRef & currentTimer,
                                                                       const size_t & from,
                                                                       const size_t & to)
{
    const NDIndex<DIM> currentAdd = getExtraMarginCurrent();

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

    FieldPatch<double> myJx(dom);
    myJx.setOrigin(origin + Vector_t(dom[0].first() * dx(0),
                                     dom[1].first() * dx(1)));
    myJx.setSpacing(dx);

    const ParticleAttrib< Vector_t > & R = bunch.R;
    const ParticleAttrib< Vector_t > & oldR = bunch.oldR;
    const ParticleAttrib< double > & Q = bunch.Q;

    double current[6];
    std::fill(current, current + 6, 0.0);

    Utils::CompanionSorterItem xy_crossings[10];

    for (size_t i = from; i < to; ++ i) {
        const double & q = Q[i];

        const Vector_t Dx = R[i] - oldR[i];
        const Vector_t delta = Dx * oneoverdx;
        const double coveredDistance = sqrt(dot(Dx, Dx));

        size_t numXs = Utils::calcDualCellBorderCrossings(xy_crossings, oldR[i], R[i], dx, origin, uni);
        for (size_t j = 0; j < numXs - 1; ++ j) {
            const Vector_t xi = xy_crossings[j].subStep * oneoverdx;
            const Vector_t xf = xy_crossings[j+1].subStep * oneoverdx;
            const Vector_t segment = xy_crossings[j+1].subStep - xy_crossings[j].subStep;
            const double tau = sqrt(dot(segment, segment)) / coveredDistance;
            const Vector_t center = (xi + xf) * onehalf;
            const int I = static_cast<int>(floor(center(0) + onehalf));
            const int J = static_cast<int>(floor(center(1) + onehalf));
            const Vector_t a(I + onehalf - xi(0), J + onehalf - xi(1));
            const Vector_t b(I - xi(0), J - xi(1));

            getCurrentX(current, tau, a, b, q * oneoverdxdt, delta);

            myJx(I-1,J-1) += current[0];
            myJx(I,  J-1) += current[1];
            myJx(I-1,J)   += current[2];
            myJx(I,  J)   += current[3];
            myJx(I-1,J+1) += current[4];
            myJx(I,  J+1) += current[5];

        }
    }

    Timings::stopTimer(currentTimer);

    return myJx;
}

FieldPatch<double> FirstOrderShapeFunction::getCurrentDensityBaseImplY(VField_Edge_t& JFD,
                                                                       const PartBunch & bunch,
                                                                       const double & dt,
                                                                       const Timings::TimerRef & currentTimer,
                                                                       const size_t & from,
                                                                       const size_t & to)
{
    const NDIndex<DIM> currentAdd = getExtraMarginCurrent();

    NDIndex<DIM> lPDom = bunch.getLocalPDomainInclOld(JFD.get_mesh());
    NDIndex<DIM> dom;
    for (unsigned int d = 0; d < DIM; ++ d) {
        dom[d] = Index(lPDom[d].first() + currentAdd[1 - d].first(),
                       lPDom[d].last() + currentAdd[1 - d].last());
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

    double current[6];
    std::fill(current, current + 6, 0.0);

    Utils::CompanionSorterItem xy_crossings[10];
    for (size_t i = from; i < to; ++ i) {
        const double & q = Q[i];

        const Vector_t Dx = R[i] - oldR[i];
        const Vector_t delta = Dx * oneoverdx;
        const double coveredDistance = sqrt(dot(Dx, Dx));

        size_t numXs = Utils::calcDualCellBorderCrossings(xy_crossings, oldR[i], R[i], dx, origin, uni);
        for (size_t j = 0; j < numXs - 1; ++ j) {
            const Vector_t xi = xy_crossings[j].subStep * oneoverdx;
            const Vector_t xf = xy_crossings[j+1].subStep * oneoverdx;
            const Vector_t segment = xy_crossings[j+1].subStep - xy_crossings[j].subStep;
            const double tau = sqrt(dot(segment, segment)) / coveredDistance;
            const Vector_t center = (xi + xf) * onehalf;
            const int I = static_cast<int>(floor(center(0) + onehalf));
            const int J = static_cast<int>(floor(center(1) + onehalf));
            const Vector_t a(I + onehalf - xi(0), J + onehalf - xi(1));
            const Vector_t b(I - xi(0), J - xi(1));

            getCurrentY(current, tau, a, b, q * oneoverdxdt, delta);

            myJy(I-1,J-1) += current[0];
            myJy(I,  J-1) += current[1];
            myJy(I+1,J-1) += current[2];
            myJy(I-1,J)   += current[3];
            myJy(I  ,J)   += current[4];
            myJy(I+1,J)   += current[5];

        }
    }


    Timings::stopTimer(currentTimer);

    return myJy;
}

void FirstOrderShapeFunction::getChargeDensityImpl(SField_Vert_t& rho,
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
        const Vector_t neR = (R[i] - origin) / dx;
        const int I = static_cast<int>(floor(neR(0) + 0.5));
        const int J = static_cast<int>(floor(neR(1) + 0.5));

        const double a = neR(0) - I;
        const double b = neR(1) - J;

        const double sxm = q * (0.5 - a)*(0.5 - a) * 0.5;
        const double sxp = q * (0.5 + a)*(0.5 + a) * 0.5;
        const double sx = q - sxm - sxp;
        const double sym = (0.5 - b)*(0.5 - b) * 0.5;
        const double syp = (0.5 + b)*(0.5 + b) * 0.5;
        const double sy = 1.0 - sym - syp;

        myRho(I-1, J-1) += sxm * sym;
        myRho(I,   J-1) += sx  * sym;
        myRho(I+1, J-1) += sxp * sym;
        myRho(I-1, J  ) += sxm * sy;
        myRho(I,   J  ) += sx  * sy;
        myRho(I+1, J  ) += sxp * sy;
        myRho(I-1, J+1) += sxm * syp;
        myRho(I,   J+1) += sx  * syp;
        myRho(I+1, J+1) += sxp * syp;
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

void FirstOrderShapeFunction::getChargeDensityDiffImpl(SField_Vert_t& rho,
                                                       const PartBunch & bunch,
                                                       const Timings::TimerRef & chargeTimer,
                                                       const Timings::TimerRef & commTimer)
{
    NDIndex<DIM> elem;
    Mesh_t & mesh = rho.get_mesh();
    const NDIndex<DIM> lDom = rho.getLayout().getLocalNDIndex();
    const NDIndex<DIM> & lPDom = bunch.getLocalPDomainInclOld(mesh);
    const NDIndex<DIM> & gPDom = bunch.getGlobalPDomainInclOld(mesh);
    const Index vertAdd = getExtraMarginVert();
    const Vector_t dx(mesh.get_meshSpacing(0),
                      mesh.get_meshSpacing(1));
    const Vector_t origin = mesh.get_origin();

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
        int I = static_cast<int>(floor(neR(0) + 0.5));
        int J = static_cast<int>(floor(neR(1) + 0.5));

        double a = neR(0) - I;
        double b = neR(1) - J;

        double sxm = q * (0.5 - a)*(0.5 - a) * 0.5;
        double sxp = q * (0.5 + a)*(0.5 + a) * 0.5;
        double sx = q - sxm - sxp;
        double sym = (0.5 - b)*(0.5 - b) * 0.5;
        double syp = (0.5 + b)*(0.5 + b) * 0.5;
        double sy = 1.0 - sym - syp;

        myRho(I-1, J-1) += sxm * sym;
        myRho(I,   J-1) += sx  * sym;
        myRho(I+1, J-1) += sxp * sym;
        myRho(I-1, J  ) += sxm * sy;
        myRho(I,   J  ) += sx  * sy;
        myRho(I+1, J  ) += sxp * sy;
        myRho(I-1, J+1) += sxm * syp;
        myRho(I,   J+1) += sx  * syp;
        myRho(I+1, J+1) += sxp * syp;

        I = static_cast<int>(floor(olR(0) + 0.5));
        J = static_cast<int>(floor(olR(1) + 0.5));

        a = olR(0) - I;
        b = olR(1) - J;

        sxm = q * (0.5 - a)*(0.5 - a) * 0.5;
        sxp = q * (0.5 + a)*(0.5 + a) * 0.5;
        sx = q - sxm - sxp;
        sym = (0.5 - b)*(0.5 - b) * 0.5;
        syp = (0.5 + b)*(0.5 + b) * 0.5;
        sy = 1 - sym - syp;

        myRho(I-1, J-1) -= sxm * sym;
        myRho(I,   J-1) -= sx  * sym;
        myRho(I+1, J-1) -= sxp * sym;
        myRho(I-1, J  ) -= sxm * sy;
        myRho(I,   J  ) -= sx  * sy;
        myRho(I+1, J  ) -= sxp * sy;
        myRho(I-1, J+1) -= sxm * syp;
        myRho(I,   J+1) -= sx  * syp;
        myRho(I+1, J+1) -= sxp * syp;
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

void FirstOrderShapeFunction::getCurrentX(double J[],
                                          const double & tau,
                                          const Vector_t & a,
                                          const Vector_t & b,
                                          const Vector_t & qoverdxdt,
                                          const Vector_t & delta)
{
    const Vector_t current_density(qoverdxdt(1) * delta(0),
                                   qoverdxdt(0) * delta(1));

    const static double half = 0.5;
    const static double onethird = 1.0 / 3.0;
    const static double onefourth = half * half;
    const double taupow2 = tau*tau * half;
    const double taupow3 = taupow2*tau * 2 * onethird;

    // The following integral has to be solved (see also thesis of Barthelme, page 54)
    //
    // \int_{0}^{\tau} S1(a_x-delta_x * t) * S2(b_y - delta_y * t) dt
    //
    // which will yield something of the form
    // Ax1 * Ay2 * tau + (Ax1 * By2 + Bx1 * Ay2) * tau^2/2 + (Ax1 * Cy2 + Bx1 * By2) * tau^3/3 + Bx1 * Cy2 * tau^4/4
    //
    //  where                                                                          Ax1
    //       A stands for zeroth order (in tau) factor, B for second, C for third etc--|||
    //  the direction from the neares grid point thus either x or y (or z when in 3D)---||
    //                                                  1 stands for S1, 2 for S2 etc----|
    // These factors are computed first since they can be reused several times
    // The above names are completed by a code that shows for which node in the vicinity of the particle
    // the factor is computed. Thus the complete name could be                                     Ax1Im1 where
    // here I, the nearest grid point in x direction, is meant (which is a duplicate information)-----||
    // the number of cells away from the nearest grid point, where m stands for minus and p for plus---|
    // The current density is the computed by a scheme as below:
    // JxIm1 = current_density(0) * (Bx1Im1 * Cy2Jm1                     * pow(tau,4) / 4 +
    //                               (Bx1Im1 * By2Jm1 + Ax1Im1 * Cy2Jm1) * pow(tau,3) / 3 +
    //                               (Bx1Im1 * Ay2Jm1 + Ax1Im1 * By2Jm1) * pow(tau,2) / 2 +
    //                               (                  Ax1Im1 * Ay2Jm1) * tau);
    // This results in quite a substantial number of multiplications. This number can be reduced by
    // multiplying the factors by tau.
    // The factors of S1 are then Ax1I --> Ax1I,
    //                            Bx1I --> Bx1I * tau
    // and the factors of S2: Ay2J --> Ay2J * tau,
    //                        By2J --> By2J * tau^2 / 2
    //                        Cy2J --> Cy2J * tau^3 / 3
    //
    // The result is then
    // JxIm1 = current_density(0) * ((Ax1Im1 + Bx1Im1) * (Ay2Jm1 + By2Jm1 + Cy2Jm1)
    //                               - Bx1Im1 * (Ay2Jm1 / 2 + By2Jm1 / 3 + Cy2Jm1 / 4));

    const double Ax1Im1 = a(0);
    const double Bx1Im1 = -delta(0) * tau;
    const double Ax1I = 1.0 - a(0);
    const double Bx1I = delta(0) * tau;

    const double Ay2Jm1 = half* (half+ b(1))*(half+ b(1)) * tau;
    const double By2Jm1 = -(half+ b(1)) * delta(1) * taupow2;
    const double Cy2Jm1 = half* delta(1)*delta(1) * taupow3;
    const double Ay2J = (0.75 - b(1)*b(1)) * tau;
    const double By2J = 2 * b(1) * delta(1) * taupow2;
    const double Cy2J = -delta(1)*delta(1) * taupow3;
    const double Ay2Jp1 = half* (half- b(1))*(half- b(1)) * tau;
    const double By2Jp1 = (half- b(1)) * delta(1) * taupow2;
    const double Cy2Jp1 = half* delta(1)*delta(1) * taupow3;

    // Jx: I-1, J-1
    J[0] = current_density(0) * ((Ax1Im1 + Bx1Im1) * (Ay2Jm1 + By2Jm1 + Cy2Jm1)
                                 - Bx1Im1 * (half * Ay2Jm1 + onethird * By2Jm1 + onefourth * Cy2Jm1));

    // Jx: I, J-1
    J[1] = current_density(0) * ((Ax1I + Bx1I) * (Ay2Jm1 + By2Jm1 + Cy2Jm1)
                                 - Bx1I * (half * Ay2Jm1 + onethird * By2Jm1 + onefourth * Cy2Jm1));

    // Jx: I-1, J
    J[2] = current_density(0) * ((Ax1Im1 + Bx1Im1) * (Ay2J + By2J + Cy2J)
                                 - Bx1Im1 * (half * Ay2J + onethird * By2J + onefourth * Cy2J));

    // Jx: I, J
    J[3] = current_density(0) * ((Ax1I + Bx1I) * (Ay2J + By2J + Cy2J)
                                 - Bx1I * (half * Ay2J + onethird * By2J + onefourth * Cy2J));

    // Jx: I-1, J+1
    J[4] = current_density(0) * ((Ax1Im1 + Bx1Im1) * (Ay2Jp1 + By2Jp1 + Cy2Jp1)
                                 - Bx1Im1 * (half * Ay2Jp1 + onethird * By2Jp1 + onefourth * Cy2Jp1));

    // Jx: I, J+1
    J[5] = current_density(0) * ((Ax1I + Bx1I) * (Ay2Jp1 + By2Jp1 + Cy2Jp1)
                                 - Bx1I * (half * Ay2Jp1 + onethird * By2Jp1 + onefourth * Cy2Jp1));
}

void FirstOrderShapeFunction::getCurrentY(double J[],
                                          const double & tau,
                                          const Vector_t & a,
                                          const Vector_t & b,
                                          const Vector_t & qoverdxdt,
                                          const Vector_t & delta)
{
    const Vector_t current_density(qoverdxdt(1) * delta(0),
                                   qoverdxdt(0) * delta(1));

    const static double half = 0.5;
    const static double onethird = 1.0 / 3.0;
    const static double onefourth = half * half;
    const double taupow2 = tau*tau * half;
    const double taupow3 = taupow2*tau * 2 * onethird;

    const double Ax2Im1 = half* (half+ b(0))*(half+ b(0)) * tau;
    const double Bx2Im1 = -(half+ b(0)) * delta(0) * taupow2;
    const double Cx2Im1 = half* delta(0)*delta(0) * taupow3;
    const double Ax2I = (0.75 - b(0)*b(0)) * tau;
    const double Bx2I = 2 * b(0) * delta(0) * taupow2;
    const double Cx2I = -delta(0)*delta(0) * taupow3;
    const double Ax2Ip1 = half* (half- b(0))*(half- b(0)) * tau;
    const double Bx2Ip1 = (half- b(0)) * delta(0) * taupow2;
    const double Cx2Ip1 = half* delta(0)*delta(0) * taupow3;

    const double Ay1Jm1 = a(1);
    const double By1Jm1 = -delta(1) * tau;
    const double Ay1J = 1.0 - a(1);
    const double By1J = delta(1) * tau;

    // Jy: I-1, J-1
    J[0] = current_density(1) * ((Ay1Jm1 + By1Jm1) * (Ax2Im1 + Bx2Im1 + Cx2Im1)
                                 - By1Jm1 * (half * Ax2Im1 + onethird * Bx2Im1 + onefourth * Cx2Im1));

    // Jy: I, J-1
    J[1] = current_density(1) * ((Ay1Jm1 + By1Jm1) * (Ax2I + Bx2I + Cx2I)
                                 - By1Jm1 * (half * Ax2I + onethird * Bx2I + onefourth * Cx2I));

    // Jy: I+1, J-1
    J[2] = current_density(1) * ((Ay1Jm1 + By1Jm1) * (Ax2Ip1 + Bx2Ip1 + Cx2Ip1)
                                 - By1Jm1 * (half * Ax2Ip1 + onethird * Bx2Ip1 + onefourth * Cx2Ip1));

    // Jy: I-1, J
    J[3] = current_density(1) * ((Ay1J + By1J) * (Ax2Im1 + Bx2Im1 + Cx2Im1)
                                 - By1J * (half * Ax2Im1 + onethird * Bx2Im1 + onefourth * Cx2Im1));

    // Jy: I, J
    J[4] = current_density(1) * ((Ay1J + By1J) * (Ax2I + Bx2I + Cx2I)
                                 - By1J * (half * Ax2I + onethird * Bx2I + onefourth * Cx2I));

    // Jy: I+1, J
    J[5] = current_density(1) * ((Ay1J + By1J) * (Ax2Ip1 + Bx2Ip1 + Cx2Ip1)
                                 - By1J * (half * Ax2Ip1 + onethird * Bx2Ip1 + onefourth * Cx2Ip1));
}
