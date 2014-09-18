/***************************************************************************
                            PartBunch.cpp
                         -------------------
    begin                : Mon Oct 5 2009
    copyright            : (C) 2009 by Christof Kraus
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

#include "ShapeFunction.hh"
#include "PartBunch.hh"
#include "Commands/BinaryVTKFile.hh"
#include "Commands/ASCIIVTKFile.hh"
#include "Distribution.hh"
#include "FieldPatch.hh"
#include "Communicator.hh"
#include "Bend.hh"
#include "Timings.hh"
#include "utils.hh"

#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>

#include <algorithm>
#include <numeric>
#include "mpi.h"

extern std::ofstream dbg;
#define DBGOUT dbg << "PartBunch.cpp: " << __LINE__ << "\t"
#define noMASTERSLAVEP2P

bool ascSecond(const std::pair<unsigned long, double> & a, const std::pair<unsigned long, double> & b) {
    return a.second < b.second;
}

bool descSecond(const std::pair<unsigned long, double> & a, const std::pair<unsigned long, double> & b) {
    return a.second > b.second;
}

PartBunch::PartBunch(FieldLayout<DIM> & edgeFL,
                     FieldLayout<DIM> & cellFL,
                     Mesh_t & mesh,
                     Distribution & distro,
                     const Bend & bend,
                     const double & charge_mass_ratio):
    _pl(NULL),
    _dt(0.0),
    _charge(distro.getTotalCharge() / distro.size()),
    _charge_mass_ratio(charge_mass_ratio),
    _momentum(0.0),
    _localNums(Ippl::getNodes(), 0),
    _localNumsSet(false),
    _bend(bend),
    _mesh(mesh),
    _iAmSlave(false),
    _iAmMaster(false)
{
    _others_timer = Timings::getTimer("other stuff");
    Timings::startTimer(_others_timer);

#ifdef GETUNIFORMLAYOUT
    _pl = new playout_t();
#else
    _pl = new playout_t(edgeFL, mesh);
#endif

    initialize(_pl);

    addAttribute(Ef);
    addAttribute(Hf);
    addAttribute(M);

    _push_timer = Timings::getTimer("PBunch::push");
    _fieldInterp_timer = Timings::getTimer("field interp");
    _bounds_timer = Timings::getTimer("bounds");
    _comm_timer = Timings::getTimer("field interp comm");
    _rebalance_timer = Timings::getTimer("PBunch::rebalance");
    _statParamTimer = IpplTimings::getTimer("Statistics");
    // _conservation_timer = Timings::getTimer("PBunch::calcCC");

    Utils::getLocalDomains(edgeFL, _localFDomains);
    _localPDomains.resize(Ippl::getNodes());

    Vector_t dx(_mesh.get_meshSpacing(0),
                _mesh.get_meshSpacing(1));

    getDistribution(distro, dx);
    Timings::stopTimer(_others_timer);
#ifdef GETUNIFORMLAYOUT
    rebalance();
#else

    Timings::startTimer(_rebalance_timer);
    this->update();

    if(false) {
        e_dim_tag decomp[] = {PARALLEL, PARALLEL, SERIAL};
        Vector_t densOrigin = mesh.get_origin();
        FieldLayout_Cell_t densFL(mesh, decomp);
        SField_Cell_t density(mesh, edgeFL, GuardCellSizes<DIM>(1));

        Timings::stopTimer(_rebalance_timer);
        ShapeFunction<IntCIC>::getNumberDensity(density, *this);
        Timings::startTimer(_rebalance_timer);

        _lPDom = calcRepartition(densFL, density);

        cellFL.Repartition(&_lPDom, &_lPDom + 1);
        // NDIndex<DIM> densGDom = edgeFL.getDomain();
        // for (unsigned int d = 0; d < DIM; ++ d) {
        //     if (densGDom[d].last() == _lPDom[d].last() + 1) {
        //         _lPDom[d] = Index(_lPDom[d].first(), _lPDom[d].last() + 1);
        //     }
        // }
        edgeFL.Repartition(&_lPDom, &_lPDom + 1);
    } else {
        this->update();
        NDIndex<DIM> eldom = edgeFL.getLocalNDIndex();
        NDIndex<DIM> egdom = edgeFL.getDomain();
        // for (unsigned int d = 0; d < DIM; ++ d) {
        //     if (eldom[d].last() == egdom[d].last()) {
        //         eldom[d] = Index(eldom[d].first(), eldom[d].last() - 1);
        //     }
        // }
        cellFL.Repartition(&eldom, &eldom + 1);
    }
    Timings::stopTimer(_rebalance_timer);

#endif
    Timings::startTimer(_others_timer);
    this->update();

    _masterStart = getLocalNum();
    _masterEnd = _masterStart;

    Timings::stopTimer(_others_timer);
    sortParticles();
    Timings::startTimer(_others_timer);

    _momentum = distro.get_pmean();
    Timings::stopTimer(_others_timer);
}

void
PartBunch::commNumLocalParticles() const
{
    boost::mpi::communicator world;
    const size_t myLocalNum = getLocalNP();

    std::vector<size_t> * lNums = const_cast<std::vector<size_t> *>(&_localNums);
    bool * lNumsSet = const_cast<bool * >(&_localNumsSet);
    boost::mpi::all_gather(world, myLocalNum, *lNums);
    *lNumsSet = true;
}

void
PartBunch::makeReady(const double & dt)
{
    Timings::startTimer(_others_timer);

    const unsigned int myNode = Ippl::myNode();

    commNumLocalParticles();
    for (size_t lnp = 0; lnp < _localNums[myNode]; ++ lnp) {
        Vector_t dX(dt * Physics::c * P[lnp] * Utils::Q_rsqrt(1 + dot(P[lnp], P[lnp])));
        oldR[lnp] = R[lnp] - dX;
    }
    Timings::stopTimer(_others_timer);
}

void PartBunch::drift_particles(const double & dt)
{
    Timings::startTimer(_push_timer);
    double fullDistance = dt * Physics::c;
    for (size_t lnp = 0; lnp < getLocalNP(); ++ lnp) {
        R[lnp] += fullDistance * P[lnp] * Utils::Q_rsqrt(1 + dot(P[lnp], P[lnp]));
    }
    Timings::stopTimer(_push_timer);
}

void PartBunch::move_by(const Vector_t & dX)
{
    Timings::startTimer(_others_timer);
    for (size_t lnp = 0; lnp < getLocalNP(); ++ lnp) {
        R[lnp] = R[lnp] + dX;
    }
    Timings::stopTimer(_others_timer);
}

void PartBunch::rotate(const Vector_t & X,
                       const double& phi)
{
    Timings::startTimer(_others_timer);
    Tenzor<double, DIM> M(cos(phi), -sin(phi),
                          sin(phi), cos(phi));

    for (size_t lnp = 0; lnp < getLocalNP(); ++ lnp) {
        R[lnp] = dot(M, R[lnp] - X) + X;
        P[lnp] = dot(M, P[lnp]);
    }
    this->update();

    // _momentum = dot(M, _momentum);
    Timings::stopTimer(_others_timer);
}

void PartBunch::push_particles(VField_Edge_t& EFD,
                               VField_Cell_t& HFD,
                               VField_Edge_t& JFD,
                               const double & dt,
                               const unsigned int & Jcomp)
{
    Timings::startTimer(_push_timer);
    size_t lN = getLocalNum();
    double halfDistance = 0.5 * dt * Physics::c;
    // const unsigned int myNode = Ippl::myNode();
    // std::vector<NDIndex<DIM> > localFDomains;
    // Utils::getLocalDomains(EFD.getLayout(), localFDomains);

    for (size_t lnp = 0; lnp < lN; ++ lnp) {
        //push particles half a time step
        oldR[lnp] = R[lnp];
        R[lnp] += halfDistance * P[lnp] * Utils::Q_rsqrt(1.0 + dot(P[lnp],P[lnp]));
    }

    Timings::stopTimer(_push_timer);

    std::vector<FieldPatch<double> > emFields(3), emFieldsMaster(3);

    getFields(emFields[0], emFields[1], emFields[2], EFD, HFD);

    Ef = 0.0;
    Hf = 0.0;

#ifdef MASTERSLAVEP2P
    if (_iAmSlave) {
        int tag = Ippl::Comm->next_tag(IPPL_APP_TAG0, IPPL_APP_CYCLE);
        MPI_Status status;
        MPI_Request req;
        std::vector<NDIndex<DIM> > localFDomains;
        Utils::getLocalDomains(EFD.getLayout(), localFDomains);
        NDIndex<DIM> & rFDom = localFDomains[_masterWorldRank];

        NDIndex<DIM> marginEdge = _shapeFunction.getExtraMarginEdge();
        Index marginCell = _shapeFunction.getExtraMarginCell();
        unsigned int sizeExPatch = (rFDom[0].length() + marginEdge[0].length());
        unsigned int sizeEyPatch = (rFDom[0].length() + marginEdge[1].length());
        unsigned int sizeHzPatch = (rFDom[0].length() + marginCell.length());
        for (unsigned int d = 1; d < DIM; ++ d) {
            unsigned int d2 = (d + 1) % DIM;
            sizeExPatch *= (rFDom[d].length() + marginEdge[d].length());
            sizeEyPatch *= (rFDom[d].length() + marginEdge[d2].length());
            sizeHzPatch *= (rFDom[d].length() + marginCell.length());
        }

        Timings::startTimer(_comm_timer);

        unsigned long sizeFieldPatches = sizeof(double) * (sizeExPatch + sizeEyPatch + sizeHzPatch);
        sizeFieldPatches += 6 * DIM * (sizeof(int) + sizeof(double)) + 3;

        std::vector<char> container(sizeFieldPatches,0);
        MPI_Irecv(&(container[0]), sizeFieldPatches, MPI_CHAR, _commRoot, tag, _commMasterSlave, &req);//&status);

        Timings::stopTimer(_comm_timer);

        emFields[0].gather(Ef, this->R, 0, IntOp_t(), 0, _masterStart);
        emFields[1].gather(Ef, this->R, 1, IntOp_t(), 0, _masterStart);
        emFields[2].gather(Hf, this->R, IntOp_t(), 0, _masterStart);

        Timings::startTimer(_comm_timer);
        MPI_Wait(&req, &status);
        size_t startIndex = 0;
        emFieldsMaster[0].deserialize(container, startIndex);
        emFieldsMaster[1].deserialize(container, startIndex);
        emFieldsMaster[2].deserialize(container, startIndex);
        Timings::stopTimer(_comm_timer);

        emFieldsMaster[0].gather(Ef, this->R, 0, IntOp_t(), _masterStart, _masterEnd);
        emFieldsMaster[1].gather(Ef, this->R, 1, IntOp_t(), _masterStart, _masterEnd);
        emFieldsMaster[2].gather(Hf, this->R, IntOp_t(), _masterStart, _masterEnd);
    } else if (_iAmMaster) {
        int tag = Ippl::Comm->next_tag(IPPL_APP_TAG0, IPPL_APP_CYCLE);
        MPI_Request req;
        std::vector<MPI_Request> requests;
        std::vector<char> container;
        int sizeComm;
        MPI_Comm_size(_commMasterSlave, &sizeComm);

        Timings::startTimer(_comm_timer);

        emFields[0].serialize(container);
        emFields[1].serialize(container);
        emFields[2].serialize(container);
        unsigned long sizeFieldPatches = container.size();

        for (int k = 0; k < sizeComm - 1; ++ k) {
            MPI_Isend(&(container[0]), sizeFieldPatches, MPI_CHAR, k, tag, _commMasterSlave, &req);
            requests.push_back(req);
        }
        Timings::stopTimer(_comm_timer);

        emFields[0].gather(Ef, this->R, 0, IntOp_t(), 0, _masterStart);
        emFields[1].gather(Ef, this->R, 1, IntOp_t(), 0, _masterStart);
        emFields[2].gather(Hf, this->R, IntOp_t(), 0, _masterStart);

        std::vector<MPI_Status> status(requests.size());
        MPI_Waitall(requests.size(), &requests[0], &status[0]);

    } else {
        emFields[0].gather(Ef, this->R, 0, IntOp_t(), 0, _masterStart);
        emFields[1].gather(Ef, this->R, 1, IntOp_t(), 0, _masterStart);
        emFields[2].gather(Hf, this->R, IntOp_t(), 0, _masterStart);
    }
#else
    if (_iAmSlave) {
        Timings::startTimer(_comm_timer);

        unsigned long sizeFieldPatches = 0;
        MPI_Bcast(&sizeFieldPatches, 1, MPI_UNSIGNED_LONG, _commRoot, _commMasterSlave);

        std::vector<char> container(sizeFieldPatches,0);
        MPI_Bcast(&(container[0]), sizeFieldPatches, MPI_CHAR, _commRoot, _commMasterSlave);

        size_t startIndex = 0;
        emFieldsMaster[0].deserialize(container, startIndex);
        emFieldsMaster[1].deserialize(container, startIndex);
        emFieldsMaster[2].deserialize(container, startIndex);
        Timings::stopTimer(_comm_timer);

        emFieldsMaster[0].gather(Ef, this->R, 0, IntOp_t(), _masterStart, _masterEnd);
        emFieldsMaster[1].gather(Ef, this->R, 1, IntOp_t(), _masterStart, _masterEnd);
        emFieldsMaster[2].gather(Hf, this->R, IntOp_t(), _masterStart, _masterEnd);
    } else if (_iAmMaster) {
        std::vector<char> container;

        Timings::startTimer(_comm_timer);

        emFields[0].serialize(container);
        emFields[1].serialize(container);
        emFields[2].serialize(container);
        unsigned long sizeFieldPatches = container.size();

        MPI_Bcast(&sizeFieldPatches, 1, MPI_UNSIGNED_LONG, _commRoot, _commMasterSlave);
        MPI_Bcast(&(container[0]), sizeFieldPatches, MPI_CHAR, _commRoot, _commMasterSlave);
        Timings::stopTimer(_comm_timer);
    }

    emFields[0].gather(Ef, this->R, 0, IntOp_t(), 0, _masterStart);
    emFields[1].gather(Ef, this->R, 1, IntOp_t(), 0, _masterStart);
    emFields[2].gather(Hf, this->R, IntOp_t(), 0, _masterStart);
#endif

    Timings::startTimer(_push_timer);
    auto externalBf = _bend.getExternalField(this->R);
    double fact1 = halfDistance * _charge_mass_ratio;
    double fact2 = fact1 * Physics::c;
    // _momentum = Vector_t(0.0);
    for (size_t lnp = 0; lnp < lN; ++ lnp) {
        // factor 0.5 is due to TETM that calls this method twice per time step
        double B = Physics::mu_0 * Hf[lnp] + 0.5 * externalBf[lnp];

        // apply kick
        Vector_t kick = fact1 * Ef[lnp];
        Vector_t um = P[lnp] + kick;
        double a = fact2 * B * Utils::Q_rsqrt(1.0 + dot(um, um));
        Vector_t s = um + a * Vector_t(um(1), -um(0));
        double tmp = 1.0 + a * a;
        um(0) = ((1.0 + a*a) * s(0) + a * s(1)) / tmp;
        um(1) = ((1.0 + a*a) * s(1) - a * s(0)) / tmp;
        P[lnp] = um + kick;

        // _momentum += P[lnp];
        // push particle
        R[lnp] += halfDistance * P[lnp] * Utils::Q_rsqrt(1.0 + dot(P[lnp],P[lnp]));
    }

    // if (getTotalNum() > 0) {
    //     DBGOUT << "using collective communication" << std::endl;
    //     reduce(&_momentum[0], &_momentum[0] + DIM, &_momentum[0], OpAddAssign());
    //     _momentum /= getTotalNum();
    // }

    Timings::stopTimer(_push_timer);

    if (Jcomp == 0) {
        FieldPatch<double> Jx;

        if (_iAmSlave) {
            FieldPatch<double> JxMaster = _shapeFunction.getCurrentDensityBaseX(JFD, *this, dt, _masterStart, _masterEnd);
            Timings::startTimer(_comm_timer);

            std::vector<char> container;
            JxMaster.serialize(container);
            unsigned long sizeFieldPatch = container.size();

            int tag = Ippl::Comm->next_tag(IPPL_APP_TAG0, IPPL_APP_CYCLE);
            MPI_Request req;
            std::vector<MPI_Request> requests;

            MPI_Isend(&(container[0]), sizeFieldPatch, MPI_CHAR, _commRoot, tag, _commMasterSlave, &req);
            requests.push_back(req);

            Timings::stopTimer(_comm_timer);

            Jx = _shapeFunction.getCurrentDensityBaseX(JFD, *this, dt, 0, _masterStart);

            std::vector<MPI_Status> status(requests.size());
            MPI_Waitall(requests.size(), &requests[0], &status[0]);
        } else if (_iAmMaster) {
            Jx = _shapeFunction.getCurrentDensityBaseX(JFD, *this, dt, 0, _masterStart);

            Timings::startTimer(_comm_timer);
            int sizeComm;
            MPI_Comm_size(_commMasterSlave, &sizeComm);

            int tag = Ippl::Comm->next_tag(IPPL_APP_TAG0, IPPL_APP_CYCLE);

            NDIndex<DIM> lFDom = JFD.getLayout().getLocalNDIndex();
            unsigned long sizeFieldPatch = sizeof(double) * (lFDom[0].length() + 2);
            for (unsigned int d = 1; d < DIM; ++ d)
                sizeFieldPatch *= (lFDom[d].length() + 2);
            sizeFieldPatch += 2 * DIM * (sizeof(int) + sizeof(double)) + 1;

            std::vector<char> container(sizeFieldPatch);
            for (int k = 0; k < sizeComm - 1; ++ k) {
                size_t startIndex = 0;
                MPI_Status status;

                MPI_Recv(&(container[0]), sizeFieldPatch, MPI_CHAR, MPI_ANY_SOURCE, tag, _commMasterSlave, &status);

                FieldPatch<double> contributedField;
                contributedField.deserialize(container, startIndex);
                Jx.add(contributedField);
            }
            Timings::stopTimer(_comm_timer);
        } else {
            Jx = _shapeFunction.getCurrentDensityBaseX(JFD, *this, dt, 0, _masterStart);
        }

        _shapeFunction.getCurrentDensityX(JFD, Jx, *this, dt);
    } else if (Jcomp == 1) {
        FieldPatch<double> Jy;

        if (_iAmSlave) {
            FieldPatch<double> JyMaster = _shapeFunction.getCurrentDensityBaseY(JFD, *this, dt, _masterStart, _masterEnd);
            Timings::startTimer(_comm_timer);

            std::vector<char> container;
            JyMaster.serialize(container);
            unsigned long sizeFieldPatch = container.size();

            int tag = Ippl::Comm->next_tag(IPPL_APP_TAG0, IPPL_APP_CYCLE);
            MPI_Request req;
            std::vector<MPI_Request> requests;

            MPI_Isend(&(container[0]), sizeFieldPatch, MPI_CHAR, _commRoot, tag, _commMasterSlave, &req);
            requests.push_back(req);

            Timings::stopTimer(_comm_timer);

            Jy = _shapeFunction.getCurrentDensityBaseY(JFD, *this, dt, 0, _masterStart);

            std::vector<MPI_Status> status(requests.size());
            MPI_Waitall(requests.size(), &requests[0], &status[0]);
        } else if (_iAmMaster) {
            Jy = _shapeFunction.getCurrentDensityBaseY(JFD, *this, dt, 0, _masterStart);

            Timings::startTimer(_comm_timer);
            int sizeComm;
            MPI_Comm_size(_commMasterSlave, &sizeComm);

            int tag0 = Ippl::Comm->next_tag(IPPL_APP_TAG0, IPPL_APP_CYCLE);

            NDIndex<DIM> lFDom = JFD.getLayout().getLocalNDIndex();
            unsigned long sizeFieldPatch = sizeof(double) * (lFDom[0].length() + 2);
            for (unsigned int d = 1; d < DIM; ++ d)
                sizeFieldPatch *= (lFDom[d].length() + 2);
            sizeFieldPatch += 2 * DIM * (sizeof(int) + sizeof(double)) + 1;

            std::vector<char> container(sizeFieldPatch);
            for (int k = 0; k < sizeComm - 1; ++ k) {
                size_t startIndex = 0;
                MPI_Status status;

                MPI_Recv(&(container[0]), sizeFieldPatch, MPI_CHAR, MPI_ANY_SOURCE, tag0, _commMasterSlave, &status);

                FieldPatch<double> contributedField;
                contributedField.deserialize(container, startIndex);
                Jy.add(contributedField);
            }
            Timings::stopTimer(_comm_timer);
        } else {
            Jy = _shapeFunction.getCurrentDensityBaseY(JFD, *this, dt, 0, _masterStart);
        }

        _shapeFunction.getCurrentDensityY(JFD, Jy, *this, dt);
    } else {
        DBGOUT << "wrong component of current requested" << std::endl;
        PAssert(Jcomp < 2);
    }

    // static int iteration = 0;
    // boost::shared_ptr<SField_Vert_t> rho = Utils::getScalarVertField(JFD);
    // _shapeFunction.getChargeDensityDiff(*rho, *this);

    // calcChargeConservation(*rho,
    //                        JFD,
    //                        dt,
    //                        iteration);
    // ++ iteration;
}

void PartBunch::getDistribution(Distribution & distro, const Vector_t & dx)
{
    if (getLocalNum())
        destroy(getLocalNum(), 0, true);
    performDestroy();

    for (size_t i = 0; i < distro.size(); ++ i) {
        create(1);
        R[i] = distro.getR(i);
        oldR[i] = distro.getOldR(i);
        P[i] = distro.getP(i);
        Q[i] = distro.getQ(i);
        M[i] = std::abs(Q[i] / _charge_mass_ratio);
    }
    distro.purge();
    _masterStart = getLocalNum();
    _masterEnd = getLocalNum();
}

void PartBunch::calcChargeConservation(const SField_Vert_t & rho,
                                       const VField_Edge_t & JFD,
                                       const double & dt,
                                       const int & iteration)
{
    JFD.fillGuardCells();
    NDIndex<DIM> elem, elemmx, elemmy;
    NDIndex<DIM> lDom = JFD.getLayout().getLocalNDIndex();
    Vector_t dx(_mesh.get_meshSpacing(0),
                _mesh.get_meshSpacing(1));

    Timings::startTimer(_conservation_timer);
    for (int j = lDom[1].first(); j <= lDom[1].last(); ++ j) {
        elem[1] = elemmx[1] = Index(j,j);
        elemmy[1] = Index(j-1,j-1);
        for (int i = lDom[0].first(); i <= lDom[0].last(); ++ i) {
            elem[0] = elemmy[0] = Index(i,i);
            elemmx[0] = Index(i-1,i-1);

            double localRho = rho.localElement(elem);
            double balance = localRho / dt;
            balance += (JFD.localElement(elem)(0) - JFD.localElement(elemmx)(0)) / dx(0);
            balance += (JFD.localElement(elem)(1) - JFD.localElement(elemmy)(1)) / dx(1);

            if (std::abs(balance) > 1e-5 * std::abs(localRho) / dt) {
                dbg << "####################################################################\n";

                dbg << "charge is not conserved on node <" << i << "," << j << ">\n"
                    << "remainder: " << balance << ",\t ratio: " << std::abs(balance/localRho) * dt << "\n"
                    << "iteration: " << iteration << "\n\n";

                dbg << "local rho: " << localRho / dt << "\n";
                dbg << std::setw(30) << std::setprecision(8) << JFD.localElement(elem)(1) / dx(1) << "\n";
                dbg << std::setw(14) << std::setprecision(8) << JFD.localElement(elemmx)(0) / dx(0);
                dbg << std::setw(30) << std::setprecision(8) << JFD.localElement(elem)(0) / dx(0) << "\n";
                dbg << std::setw(30) << std::setprecision(8) << JFD.localElement(elemmy)(1) / dx(1) << "\n\n";

                dbg << "####################################################################" << std::endl;
            }
        }
    }
    Timings::stopTimer(_conservation_timer);
}

void PartBunch::calcBeamParameters(PartBunchState & state) const {
    using Physics::c;

    Vector_t eps2, fac, rsqsum, psqsum, rpsum;
    const double m0 = Physics::m_e * 1000;

    IpplTimings::startTimer(_statParamTimer);

    const size_t locNp = this->getLocalNum();
    const double N =  static_cast<double>(this->getTotalNum());

    const double zero = 0.0;
    if(N == 0) {
        for(unsigned int i = 0 ; i < DIM; i++) {
            state._rmean(i) = 0.0;
            state._pmean(i) = 0.0;
            state._rrms(i) = 0.0;
            state._prms(i) = 0.0;
            state._epsNorm(i) = 0.0;
        }
        state._rprms = 0.0;
        state._Ekin = 0.0;

        IpplTimings::stopTimer(_statParamTimer);
        return;
    }

    calcMoments(state);

    for(unsigned int i = 0 ; i < DIM; i++) {
        state._rmean(i) = state._centroid[2 * i] / N;
        state._pmean(i) = state._centroid[(2 * i) + 1] / N;
    }

    state._centroid[0] = 0.0;
    state._centroid[2] = 0.0;
    state._centroid[1] = dot(state._pmean, state._pmean);
    state._centroid[3] = 0.0;

    for(unsigned int i = 0 ; i < DIM; i++) {
        rsqsum(i) = state._moments(2 * i, 2 * i);
        psqsum(i) = state._moments((2 * i) + 1, (2 * i) + 1) - N * state._centroid[(2 * i) + 1];
        if(psqsum(i) < 0) {
            psqsum(i) = 0;
        }
        rpsum(i) = state._moments((2 * i), (2 * i) + 1);
    }
    eps2 = (rsqsum * psqsum - rpsum * rpsum) / (N * N);
    rpsum /= N;

    for(unsigned int i = 0 ; i < DIM; i++) {
        state._rrms(i) = sqrt(rsqsum(i) / N);
        state._prms(i) = sqrt(psqsum(i) / N);
        state._epsNorm(i)  = sqrt(std::max(eps2(i), zero));
        double tmp = state._rrms(i) * state._prms(i);
        fac(i) = (tmp == 0) ? zero : 1.0 / tmp;
    }
    state._rprms = rpsum * fac;

    state._Dy = state._moments(2, 1) / N;
    state._DDy = state._moments(3, 1) / N;

    double gamma = 0.0;
    for(size_t i = 0; i < locNp; i++)
        gamma += sqrt(1.0 + dot(P[i], P[i]));

    reduce(gamma, gamma, OpAddAssign());
    gamma /= N;
    state._Ekin = (gamma - 1.0) * m0;
    // calculate energy spread
    state._dE = state._prms(0) * sqrt(state._Ekin * (state._Ekin + 2. * m0) / (1. + state._Ekin * (state._Ekin + 2. * m0) / (m0*m0)));

    IpplTimings::stopTimer(_statParamTimer);
}

void PartBunch::calcMoments(PartBunchState & state) const {

    const unsigned int locNum = this->getLocalNum();
    const unsigned int totNum = this->getTotalNum();
    double tmp;
    long double part[2 * DIM];
    Utils::KahanAccumulation loc_sums[2 * DIM];
    Utils::KahanAccumulation loc_sqr_sums[4 * DIM*DIM];
    long double loc_centroid[2 * DIM];
    long double loc_moment[4 * DIM*DIM];
    long double moments[4 * DIM*DIM];
    double limits[2 * DIM];
    double abslimits[2 * DIM];
    double cosphi;
    double sinphi;
    double avgmomentum;
    Vector_t avgP;
    Vector_t avgR;

    for (unsigned int d = 0; d < DIM; ++ d) {
        limits[d] = R[0](d);
        limits[d + DIM] = R[0](d);
    }

    for (unsigned int i = 0; i < 2 * DIM; ++i) {
        loc_sums[i].sum = 0.0;
        loc_sums[i].correction = 0.0;
        for (unsigned int j = 0; j < 2 * DIM; ++j) {
            loc_sqr_sums[i * 2 * DIM + j].sum = 0.0;
            loc_sqr_sums[i * 2 * DIM + j].correction = 0.0;
        }
    }

    std::fill(loc_centroid, loc_centroid + 2 * DIM, 0.0);
    std::fill(loc_moment, loc_moment + 4 * DIM*DIM, 0.0);

    for (unsigned long k = 0; k < locNum; ++k) {
        for (unsigned int d = 0; d < DIM; ++ d) {
            part[2 * d    ] = this->R[k](d);
            part[2 * d + 1] = this->P[k](d);

            if (this->R[k](d) < limits[d]) limits[d] = this->R[k](d);
            if (this->R[k](d) > limits[d + DIM]) limits[d + DIM] = this->R[k](d);
        }

        for (unsigned int i = 0; i < 2 * DIM; ++i) {
            loc_sums[i] = Utils::KahanSum(loc_sums[i], part[i]);
        }
    }

    for (unsigned int d = 0; d < DIM; ++ d) {
        limits[d + DIM] = -limits[d + DIM];
    }

    for (unsigned int i = 0; i < 2 * DIM; ++i) {
        loc_centroid[i] = loc_sums[i].sum;
    }

    DBGOUT << "using collective communication" << std::endl;
    MPI_Allreduce(loc_centroid, state._centroid, 2 * DIM, MPI_LONG_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(limits, abslimits, 2 * DIM, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    for (unsigned int d = 0; d < DIM; ++ d) {
        state._rmin(d) = abslimits[d];
        state._rmax(d) = -abslimits[d + DIM];
    }

    avgR = Vector_t(state._centroid[0] / totNum,
                    state._centroid[2] / totNum);
    avgP = Vector_t(state._centroid[1] / totNum,
                    state._centroid[3] / totNum);
    avgmomentum = std::sqrt(dot(avgP, avgP));
    cosphi = avgP(0) / avgmomentum;
    sinphi = avgP(1) / avgmomentum;

    for(unsigned long k = 0; k < locNum; ++k) {
        part[1] =  cosphi * this->P[k](0) + sinphi * this->P[k](1);
        part[3] = -sinphi * this->P[k](0) + cosphi * this->P[k](1);
        part[0] = this->R[k](0) - avgR(0);
        part[2] = this->R[k](1) - avgR(1);
        tmp     =  cosphi * part[0] + sinphi * part[2];
        part[2] = -sinphi * part[0] + cosphi * part[2];
        part[0] = tmp;

        for(unsigned int i = 0; i < 2 * DIM; ++i) {
            for(unsigned int j = 0; j <= i; ++j) {
                loc_sqr_sums[i * 2 * DIM + j] = Utils::KahanSum(loc_sqr_sums[i * 2 * DIM + j], part[i] * part[j]);
            }
        }
    }

    for(unsigned int i = 0; i < 2 * DIM; ++i) {
        for(unsigned int j = 0; j <= i; ++j) {
            loc_moment[i * 2 * DIM + j] = loc_sqr_sums[i * 2 * DIM + j].sum;
            loc_moment[j * 2 * DIM + i] = loc_moment[i * 2 * DIM + j];
        }
    }

    DBGOUT << "using collective communication" << std::endl;
    MPI_Allreduce(loc_moment, moments, 4 * DIM*DIM, MPI_LONG_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    // reduce(&(loc_moment[0]), &(loc_moment[0]) + 4 * DIM * DIM,
    //        &(moments[0]), OpAddAssign());

    for(unsigned int i = 0; i < 2 * DIM; ++i) {
        for(unsigned int j = 0; j <= i; ++j) {
            state._moments(i, j) = moments[i * 2 * DIM + j];
        }
    }
}

void PartBunch::get_rmean(Vector_t & spos) const
{
    Vector_t lspos(0.0);
    for (size_t lnp = 0; lnp < getLocalNum(); ++ lnp) {
        lspos += R[lnp];
    }
    DBGOUT << "using collective communication" << std::endl;
    reduce(&lspos[0], &lspos[0] + DIM, &lspos[0], OpAddAssign());
    spos = Vector_t(lspos[0] / getTotalNum(),
                    lspos[1] / getTotalNum());
}

double PartBunch::get_qtotal() const
{
    double qtotal = 0.0;
    for (size_t lnp = 0; lnp < getLocalNum(); ++ lnp) {
        qtotal += Q[lnp];
    }
    DBGOUT << "using collective communication" << std::endl;
    reduce(qtotal, qtotal, OpAddAssign());
    return qtotal;
}

void PartBunch::save(const VField_Edge_t & field, const std::string & fname, const int comp)
{
    NDIndex<DIM> lDom = field.getLayout().getLocalNDIndex();
    NDIndex<DIM> elem;

    std::ofstream out(fname.c_str());

    for (int j = lDom[1].first(); j <= lDom[1].last(); ++ j) {
        elem[1] = Index(j, j);
        for (int i = lDom[0].first(); i <= lDom[0].last(); ++ i) {
            elem[0] = Index(i, i);

            out << std::setw(18) << std::setprecision(10) << field.localElement(elem)(comp);
        }
        out << "\n";
    }
    out.close();
}

NDIndex<DIM> PartBunch::getLocalPDomain(const Mesh_t & mesh) const
{
    Timings::startTimer(_bounds_timer);
    const Vector_t origin = mesh.get_origin();
    const Vector_t dx(mesh.get_meshSpacing(0),
                      mesh.get_meshSpacing(1));

    Vector_t lRmin(0.0), lRmax(0.0);
    if (_masterStart > 0)
        getBounds(lRmin, lRmax, 0, _masterStart);

    if (_iAmSlave) {
        Vector_t lRminMaster, lRmaxMaster;
        getBounds(lRminMaster, lRmaxMaster, _masterStart, _masterEnd);
        std::vector<double> sendbuf(2 * DIM, 0.0);
        for (unsigned int d = 0; d < DIM; ++ d) {
            sendbuf[d] = lRminMaster(d);
            sendbuf[d + DIM] = -lRmaxMaster(d);
        }
        MPI_Reduce(&(sendbuf[0]), &(sendbuf[0]), 2 * DIM, MPI_DOUBLE, MPI_MIN, _commRoot, _commMasterSlave);
    } else if (_iAmMaster) {
        std::vector<double> receivebuf(2 * DIM, 0.0);
        for (unsigned int d = 0; d < DIM; ++ d) {
            receivebuf[d] = lRmin(d);
            receivebuf[d + DIM] = -lRmax(d);
        }
        MPI_Reduce(MPI_IN_PLACE, &(receivebuf[0]), 2 * DIM, MPI_DOUBLE, MPI_MIN, _commRoot, _commMasterSlave);
        for (unsigned int d = 0; d < DIM; ++ d) {
            lRmin(d) = receivebuf[d];
            lRmax(d) = -receivebuf[d + DIM];
        }
    }

    NDIndex<DIM> ret;
    for (unsigned int d = 0; d < DIM; ++ d) {
        ret[d] = Index(std::floor((lRmin(d) - origin(d)) / dx(d)),
                       std::floor((lRmax(d) - origin(d)) / dx(d)));
    }

    Timings::stopTimer(_bounds_timer);
    return ret;
}

NDIndex<DIM> PartBunch::getLocalPDomainInclOld(const Mesh_t & mesh) const
{
    if (getLocalNum() == 0) return NDIndex<DIM>(Index(0), Index(0));

    Timings::startTimer(_bounds_timer);
    const Vector_t origin = mesh.get_origin();
    const Vector_t dx(mesh.get_meshSpacing(0),
                      mesh.get_meshSpacing(1));

    Vector_t lRmin, lRmax;
    if (_masterStart > 0)
        getBoundsInclOld(lRmin, lRmax, 0, _masterStart);

    if (_iAmSlave) {
        Vector_t lRminMaster, lRmaxMaster;
        getBoundsInclOld(lRminMaster, lRmaxMaster, _masterStart, _masterEnd);
        std::vector<double> sendbuf(2 * DIM, 0.0);
        for (unsigned int d = 0; d < DIM; ++ d) {
            sendbuf[d] = lRminMaster(d);
            sendbuf[d + DIM] = -lRmaxMaster(d);
        }
        MPI_Reduce(&(sendbuf[0]), &(sendbuf[0]), 2 * DIM, MPI_DOUBLE, MPI_MIN, _commRoot, _commMasterSlave);
    } else if (_iAmMaster) {
        std::vector<double> receivebuf(2 * DIM, 0.0);
        for (unsigned int d = 0; d < DIM; ++ d) {
            receivebuf[d] = lRmin(d);
            receivebuf[d + DIM] = -lRmax(d);
        }
        MPI_Reduce(MPI_IN_PLACE, &(receivebuf[0]), 2 * DIM, MPI_DOUBLE, MPI_MIN, _commRoot, _commMasterSlave);
        for (unsigned int d = 0; d < DIM; ++ d) {
            lRmin(d) = receivebuf[d];
            lRmax(d) = -receivebuf[d + DIM];
        }
    }

    NDIndex<DIM> ret;
    for (unsigned int d = 0; d < DIM; ++ d) {
        lRmin(d) = (lRmin(d) - origin(d)) / dx(d);
        lRmax(d) = (lRmax(d) - origin(d)) / dx(d);
        ret[d] = Index(std::floor(lRmin(d)), std::floor(lRmax(d)));
    }
    Timings::stopTimer(_bounds_timer);
    return ret;
}

NDIndex<DIM> PartBunch::getGlobalPDomain(const Mesh_t & mesh) const
{
    Timings::startTimer(_bounds_timer);

    const unsigned int numNodes = Ippl::getNodes();
    NDIndex<DIM> ret = getLocalPDomain(mesh);
    const auto localPDomains = Communicator::getAllLocalDomains(ret);

    if (!_localNumsSet) commNumLocalParticles();

    unsigned int kk = 0;
    while (kk < numNodes && _localNums[kk] == 0) ++ kk;
    ret = localPDomains[kk];

    for (; kk < numNodes; ++ kk) {
        if (_localNums[kk] == 0) continue;

        const NDIndex<DIM> & dom = localPDomains[kk];
        for (unsigned int d = 0; d < DIM; ++ d) {
            int lower = std::min(ret[d].first(), dom[d].first());
            int upper = std::max(ret[d].last(), dom[d].last());
            ret[d] = Index(lower, upper);
        }
    }

    Timings::stopTimer(_bounds_timer);
    return ret;
}

NDIndex<DIM> PartBunch::getGlobalPDomainInclOld(const Mesh_t & mesh) const
{
    Timings::startTimer(_bounds_timer);

    const unsigned int numNodes = Ippl::getNodes();
    NDIndex<DIM> ret = getLocalPDomainInclOld(mesh);
    const auto localPDomains = Communicator::getAllLocalDomains(ret);

    if (!_localNumsSet) commNumLocalParticles();

    unsigned int kk = 0;
    while (kk < numNodes && _localNums[kk] == 0) ++ kk;
    ret = localPDomains[kk];

    for (; kk < numNodes; ++ kk) {
        if (_localNums[kk] == 0) continue;

        const NDIndex<DIM> & dom = localPDomains[kk];
        for (unsigned int d = 0; d < DIM; ++ d) {
            int lower = std::min(ret[d].first(), dom[d].first());
            int upper = std::max(ret[d].last(), dom[d].last());
            ret[d] = Index(lower, upper);
        }
    }

    Timings::stopTimer(_bounds_timer);
    return ret;
}

void PartBunch::getFields(FieldPatch<double> & Ex,
                          FieldPatch<double> & Ey,
                          FieldPatch<double> & Hz,
                          const VField_Edge_t& EFD,
                          const VField_Cell_t& HFD)
{
    Timings::startTimer(_comm_timer);
    const unsigned int numNodes = Ippl::getNodes();
    const Index cellAdd = ShapeFunction<IntOp_t>::getExtraMarginCell();
    const NDIndex<DIM> edgeAdd = ShapeFunction<IntOp_t>::getExtraMarginEdge();

    _lPDom = getLocalPDomain(_mesh);
    const std::vector<NDIndex<DIM> > localPDomains = Communicator::getAllLocalDomains(_lPDom);
    std::vector<NDIndex<DIM> > localPDomains_hz(localPDomains);
    std::vector<NDIndex<DIM> > localPDomains_ex(localPDomains);
    std::vector<NDIndex<DIM> > localPDomains_ey(localPDomains);

    _gPDom = _lPDom;
    for (unsigned int k = 0; k < localPDomains.size(); ++ k) {
        for (unsigned int d = 0; d < DIM; ++ d) {
            int lower = std::min(_gPDom[d].first(), localPDomains[k][d].first());
            int upper = std::max(_gPDom[d].last(), localPDomains[k][d].last());
            _gPDom[d] = Index(lower, upper);

            NDIndex<DIM> & lhz = localPDomains_hz[k];
            lhz[d] = Index(lhz[d].first() + cellAdd.first(), lhz[d].last() + cellAdd.last());
            NDIndex<DIM> & lex = localPDomains_ex[k];
            lex[d] = Index(lex[d].first() + edgeAdd[d].first(), lex[d].last() + edgeAdd[d].last());
            NDIndex<DIM> & ley = localPDomains_ey[k];
            ley[d] = Index(ley[d].first() + edgeAdd[(1 + d) % DIM].first(), ley[d].last() + edgeAdd[(1 + d) % DIM].last());
        }
    }
    Timings::stopTimer(_comm_timer);

    Timings::startTimer(_fieldInterp_timer);
    const NDIndex<DIM> & lFDom_e = EFD.getLayout().getLocalNDIndex();
    const NDIndex<DIM> & lFDom_h = HFD.getLayout().getLocalNDIndex();

    const int gcs = GUARDCELLSIZE;

    int lowerI = std::max(lFDom_h[0].first() - gcs, _gPDom[0].first() + cellAdd.min());
    int upperI = std::min(lFDom_h[0].last() + gcs,  _gPDom[0].last() + cellAdd.max());
    int signI = lowerI <= upperI ? 1 : upperI - lowerI;

    int lowerJ = std::max(lFDom_h[1].first() - gcs, _gPDom[1].first() + cellAdd.min());
    int upperJ = std::min(lFDom_h[1].last() + gcs,  _gPDom[1].last() + cellAdd.max());
    int signJ = lowerJ <= upperJ ? 1 : upperJ - lowerJ;

    Hz.resize(NDIndex<DIM>(Index(lowerI, upperI, signI),
                           Index(lowerJ, upperJ, signJ)));

    NDIndex<DIM> elem;
    for (int j = lowerJ; j <= upperJ; ++ j) {
        elem[1] = Index(j, j);
        for (int i = lowerI; i <= upperI; ++ i) {
            elem[0] = Index(i, i);

            Hz(i, j) = HFD.localElement(elem)(0) + HFD.localElement(elem)(1);
        }
    }

    lowerI = std::max(lFDom_e[0].first() - gcs, _gPDom[0].first() + edgeAdd[0].min());
    upperI = std::min(lFDom_e[0].last() + gcs,  _gPDom[0].last() + edgeAdd[0].max());
    signI = lowerI <= upperI ? 1 : upperI - lowerI;

    lowerJ = std::max(lFDom_e[1].first() - gcs, _gPDom[1].first() + edgeAdd[1].min());
    upperJ = std::min(lFDom_e[1].last() + gcs,  _gPDom[1].last() + edgeAdd[1].max());
    signJ = lowerJ <= upperJ ? 1 : upperJ - lowerJ;

    Ex.resize(NDIndex<DIM>(Index(lowerI, upperI, signI),
                           Index(lowerJ, upperJ, signJ)));

    for (int j = lowerJ; j <= upperJ; ++ j) {
        elem[1] = Index(j, j);
        for (int i = lowerI; i <= upperI; ++ i) {
            elem[0] = Index(i, i);

            Ex(i, j) = EFD.localElement(elem)(0);
        }
    }

    lowerI = std::max(lFDom_e[0].first() - gcs, _gPDom[0].first() + edgeAdd[1].min());
    upperI = std::min(lFDom_e[0].last() + gcs,  _gPDom[0].last() + edgeAdd[1].max());
    signI = lowerI <= upperI ? 1 : upperI - lowerI;

    lowerJ = std::max(lFDom_e[1].first() - gcs, _gPDom[1].first() + edgeAdd[0].min());
    upperJ = std::min(lFDom_e[1].last() + gcs,  _gPDom[1].last() + edgeAdd[0].max());
    signJ = lowerJ <= upperJ ? 1 : upperJ - lowerJ;

    Ey.resize(NDIndex<DIM>(Index(lowerI, upperI, signI),
                           Index(lowerJ, upperJ, signJ)));

    for (int j = lowerJ; j <= upperJ; ++ j) {
        elem[1] = Index(j, j);
        for (int i = lowerI; i <= upperI; ++ i) {
            elem[0] = Index(i, i);

            Ey(i, j) = EFD.localElement(elem)(1);
        }
    }

    const Vector_t origin = _mesh.get_origin();
    Vector_t dx(_mesh.get_meshSpacing(0),
                _mesh.get_meshSpacing(1));
    Vector_t patchOrigin = origin + Vector_t(_lPDom[0].first() * dx(0),
                                             _lPDom[1].first() * dx(1));

    Hz.setOrigin(patchOrigin + cellAdd.min() * dx);
    Hz.setSpacing(dx);

    Ex.setOrigin(patchOrigin + Vector_t(edgeAdd[0].first() * dx(0),
                                        edgeAdd[1].first() * dx(1)));
    Ex.setSpacing(dx);

    Ey.setOrigin(patchOrigin + Vector_t(edgeAdd[1].first() * dx(0), +
                                        edgeAdd[0].first() * dx(1)));
    Ey.setSpacing(dx);

    Timings::stopTimer(_fieldInterp_timer);

    Timings::startTimer(_comm_timer);
    std::vector<NDIndex<DIM> > localFDomains_ex;
    Utils::getLocalDomains(EFD.getLayout(), localFDomains_ex);

    std::vector<NDIndex<DIM> > localFDomains_ey(localFDomains_ex);

    std::vector<NDIndex<DIM> > localFDomains_hz;
    Utils::getLocalDomains(HFD.getLayout(), localFDomains_hz);

    for (unsigned int k = 0; k < numNodes; ++ k) {
        for (unsigned int d = 0; d < DIM; ++ d) {
            int lower = std::max(localFDomains_ex[k][d].first() - gcs, _gPDom[d].first() + edgeAdd[d].min());
            int upper = std::min(localFDomains_ex[k][d].last() + gcs, _gPDom[d].last() + edgeAdd[d].max());
            int sign = lower <= upper? 1 : upper - lower;
            localFDomains_ex[k][d] = Index(lower, upper, sign);

            lower = std::max(localFDomains_ey[k][d].first() - gcs, _gPDom[d].first() + edgeAdd[1-d].min());
            upper = std::min(localFDomains_ey[k][d].last() + gcs, _gPDom[d].last() + edgeAdd[1-d].max());
            sign = lower <= upper? 1 : upper - lower;
            localFDomains_ey[k][d] = Index(lower, upper, sign);

            lower = std::max(localFDomains_hz[k][d].first() - gcs, _gPDom[d].first() + cellAdd.min());
            upper = std::min(localFDomains_hz[k][d].last() + gcs, _gPDom[d].last() + cellAdd.max());
            sign = lower <= upper? 1 : upper - lower;
            localFDomains_hz[k][d] = Index(lower, upper, sign);
        }
    }

    Timings::stopTimer(_comm_timer);

    std::vector<FieldPatch<double> *> localFields = {&Hz, &Ex, &Ey};
    std::vector<std::vector<NDIndex<DIM> > > allLocalFDomains = {localFDomains_hz, localFDomains_ex, localFDomains_ey};
    std::vector<std::vector<NDIndex<DIM> > > allLocalPDomains = {localPDomains_hz, localPDomains_ex, localPDomains_ey};
    communicateFields(localFields, allLocalFDomains, allLocalPDomains);
}

void PartBunch::getFieldsNEW(FieldPatch<double> & Ex,
                             FieldPatch<double> & Ey,
                             FieldPatch<double> & Hz,
                             const VField_Edge_t& EFD,
                             const VField_Cell_t& HFD)
{
    Timings::startTimer(_fieldInterp_timer);
    std::vector<NDIndex<DIM> > localFDomains;
    const FieldLayout<DIM> & FL = EFD.getLayout();
    const NDIndex<DIM> gDom = FL.getDomain();
    const int gcs = GUARDCELLSIZE;
    const Index cellAdd = ShapeFunction<IntOp_t>::getExtraMarginCell();
    const NDIndex<DIM> edgeAdd = ShapeFunction<IntOp_t>::getExtraMarginEdge();
    const unsigned int myNode = Ippl::myNode();

    Utils::getLocalDomains(FL, localFDomains);
    std::vector<NDIndex<DIM> > origFDomains(localFDomains);
    std::vector<NDIndex<DIM> > localPDomains(localFDomains);

    int minLower = std::min(cellAdd.first(), edgeAdd[0].first());
    minLower = std::min(minLower, edgeAdd[1].first());
    int maxUpper = std::max(cellAdd.last(), edgeAdd[0].last());
    maxUpper = std::max(maxUpper, edgeAdd[1].last());
    Index maxAdd(minLower - 1, maxUpper + 1);

    int sendLayerThickness = std::max(maxAdd.last(), -maxAdd.first()) - gcs;

    for (unsigned int k = 0; k < localFDomains.size(); ++ k) {
        for (unsigned int d = 0; d < DIM; ++ d) {
            int lower = localFDomains[k][d].first();
            int upper = localFDomains[k][d].last();
            if (localFDomains[k][d].first() != gDom[d].first()) {
                lower = localFDomains[k][d].first() - gcs;
            }
            if (localFDomains[k][d].last() != gDom[d].last()) {
                upper = localFDomains[k][d].last() + gcs;
            }
            localFDomains[k][d] = Index(lower, upper);

            lower = localPDomains[k][d].first();
            upper = localPDomains[k][d].last();
            if (localPDomains[k][d].first() != gDom[d].first()) {
                lower = localPDomains[k][d].first() + maxAdd.first();
            }
            if (localPDomains[k][d].last() != gDom[d].last()) {
                upper = localPDomains[k][d].last() + maxAdd.last();
            }
            localPDomains[k][d] = Index(lower, upper);
        }
    }

    const NDIndex<DIM> & lFDom = localFDomains[myNode];
    Hz.resize(lFDom);
    Ex.resize(lFDom);
    Ey.resize(lFDom);

    NDIndex<DIM> elem;
    for (int j = lFDom[1].first(); j <= lFDom[1].last(); ++ j) {
        elem[1] = Index(j,j);
        for (int i = lFDom[0].first(); i <= lFDom[0].last(); ++ i) {
            elem[0] = Index(i,i);
            Hz(i,j) = HFD.localElement(elem)(0) + HFD.localElement(elem)(1);
            Ex(i,j) = EFD.localElement(elem)(0);
            Ey(i,j) = EFD.localElement(elem)(1);
        }
    }

    const Vector_t origin = _mesh.get_origin();
    Vector_t dx(_mesh.get_meshSpacing(0),
                _mesh.get_meshSpacing(1));
    Vector_t patchOrigin = origin + Vector_t(localFDomains[myNode][0].first() * dx(0),
                                             localFDomains[myNode][1].first() * dx(1));

    Hz.setOrigin(patchOrigin);
    Hz.setSpacing(dx);

    Ex.setOrigin(patchOrigin);
    Ex.setSpacing(dx);

    Ey.setOrigin(patchOrigin);
    Ey.setSpacing(dx);

    for (unsigned int k = 0; k < localFDomains.size(); ++ k) {
        for (unsigned int d = 0; d < DIM; ++ d) {
            int lower = origFDomains[k][d].first();
            int upper = origFDomains[k][d].last();
            if (origFDomains[k][d].first() != gDom[d].first()) {
                lower = origFDomains[k][d].first() + gcs;
            }
            if (origFDomains[k][d].last() != gDom[d].last()) {
                upper = origFDomains[k][d].last() - gcs;
            }
            localFDomains[k][d] = Index(lower, upper);
        }
    }

    Timings::stopTimer(_fieldInterp_timer);

    if (sendLayerThickness > 0) {
        std::vector<FieldPatch<double> *> localFields = {&Hz, &Ex, &Ey};
        std::vector<std::vector<NDIndex<DIM> > > allLocalFDomains = {localFDomains, localFDomains, localFDomains};
    std::vector<std::vector<NDIndex<DIM> > > allLocalPDomains = {localPDomains, localPDomains, localPDomains};
        communicateFields(localFields, allLocalFDomains, allLocalPDomains);
    }
}

void PartBunch::communicateFields(std::vector<FieldPatch<double> * > localFields,
                                  const std::vector<std::vector<NDIndex<DIM> > > & allLocalFDomains,
                                  const std::vector<std::vector<NDIndex<DIM> > > & allLocalPDomains) const
{
    Timings::startTimer(_comm_timer);

    const unsigned int numNodes = Ippl::getNodes();
    const unsigned int myNode = Ippl::myNode();

    std::vector<std::vector<char> > sendData(numNodes);
    std::vector<std::vector<char> > receiveData(numNodes);
    std::vector<MPI_Request> requests;
    int tag = Ippl::Comm->next_tag(F_GUARD_CELLS_TAG, F_TAG_CYCLE);

    for (unsigned int k = 0; k < numNodes; ++ k) {
        if (k == myNode) continue;
        std::vector<char> & data = receiveData[k];

        for (unsigned int fidx = 0; fidx < 3; ++ fidx) {
            const NDIndex<DIM> & oLFDom = allLocalFDomains[fidx][k];
            const NDIndex<DIM> & mLPDom = allLocalPDomains[fidx][myNode];

            bool inside = true;
            NDIndex<DIM> dom;
            for (unsigned int d = 0; d < DIM; ++ d) {
                int lower = std::max(mLPDom[d].first(), oLFDom[d].first());
                int upper = std::min(mLPDom[d].last(),  oLFDom[d].last());
                if (lower > upper) {
                    inside = false;
                    break;
                }
                dom[d] = Index(lower, upper);
            }
            if (!inside) continue;

            data.resize(data.size() + 1 + 2 * DIM * sizeof(int) + dom.size() * sizeof(double));
        }

        if (receiveData[k].size() == 0) continue;

        MPI_Request req;
        MPI_Irecv(&(receiveData[k][0]), receiveData[k].size(), MPI_BYTE, k, tag, MPI_COMM_WORLD, &req);
        requests.push_back(req);
    }

    for (unsigned int k = 0; k < numNodes; ++ k) {
        if (k == myNode) continue;
        unsigned int start = sendData[k].size();
        std::vector<char> & data = sendData[k];

        for (unsigned int fidx = 0; fidx < 3; ++ fidx) {
            const NDIndex<DIM> & oLPDom = allLocalPDomains[fidx][k];
            const NDIndex<DIM> & mLFDom = allLocalFDomains[fidx][myNode];

            bool inside = true;
            NDIndex<DIM> dom;
            for (unsigned int d = 0; d < DIM; ++ d) {
                int lower = std::max(oLPDom[d].first(), mLFDom[d].first());
                int upper = std::min(oLPDom[d].last(),  mLFDom[d].last());
                if (lower > upper) {
                    inside = false;
                    break;
                }
                dom[d] = Index(lower, upper);
            }
            if (!inside) continue;

            const char* buffer;

            char FType = fidx;
            data.insert(data.end(), &FType, &FType + 1);

            for (unsigned int d = 0; d < DIM; ++ d) {
                int lower = dom[d].first();
                buffer = reinterpret_cast<const char*>(&lower);
                data.insert(data.end(), buffer, buffer + sizeof(int));

                int upper = dom[d].last();
                buffer = reinterpret_cast<const char*>(&upper);
                data.insert(data.end(), buffer, buffer + sizeof(int));
            }

            double * rawData = new double[dom.size()];
            const FieldPatch<double> & fld = *(localFields[fidx]);
            unsigned int ii = 0;
            for (int j = dom[1].first(); j <= dom[1].last(); ++ j) {
                for (int i = dom[0].first(); i <= dom[0].last(); ++ i) {
                    rawData[ii++] = fld(i,j);
                }
            }
            buffer = reinterpret_cast<const char*>(rawData);
            data.insert(data.end(), buffer, buffer + dom.size() * sizeof(double));

            delete[] rawData;
        }

        if (sendData[k].size() - start == 0) continue;

        MPI_Send(&(sendData[k][0]), sendData[k].size(), MPI_BYTE, k, tag, MPI_COMM_WORLD);
    }

    std::vector<MPI_Status> status(requests.size());
    MPI_Waitall(requests.size(), &requests[0], &status[0]);

    for (unsigned int fidx = 0; fidx < 3; ++ fidx) {
        localFields[fidx]->resize(allLocalPDomains[fidx][myNode]);
    }

    for (unsigned int k = 0; k < numNodes; ++ k) {
        if (k == myNode) continue;
        std::vector<char> & data = receiveData[k];
        unsigned int size = data.size();
        if (size == 0) continue;

        NDIndex<DIM> dom;
        unsigned it = 0;
        while (it < size) {
            unsigned int FType = data[it++];
            FieldPatch<double> & fp = *localFields[FType];

            for (unsigned int d = 0; d < DIM; ++ d) {
                int lower = *reinterpret_cast<const int*>(&data[it]);
                it += sizeof(int);
                int upper = *reinterpret_cast<const int*>(&data[it]);
                it += sizeof(int);
                dom[d] = Index(lower, upper);
            }

            unsigned int sizeData = dom.size();

            double * rawData = new double[sizeData];
            for (unsigned int l = 0; l < sizeData; ++ l) {
                rawData[l] = *reinterpret_cast<const double*>(&data[it]);
                it += sizeof(double);
            }

            int ii = 0;
            for (int j = dom[1].first(); j <= dom[1].last(); ++ j) {
                for (int i = dom[0].first(); i <= dom[0].last(); ++ i) {
                    fp(i,j) = rawData[ii++];
                }
            }

            delete[] rawData;
        }
    }

    Timings::stopTimer(_comm_timer);
}

void PartBunch::rebalance(bool doRedistribute)
{
    _gPDom = getGlobalPDomain(_mesh);

    Timings::startTimer(_rebalance_timer);

    e_dim_tag decomp[] = {PARALLEL, PARALLEL, SERIAL};
    Vector_t spacing(_mesh.get_meshSpacing(0), _mesh.get_meshSpacing(1));
    Vector_t origin = _mesh.get_origin();

    const Index CellAdd = ZerothOrderShapeFunction::getExtraMarginCell();
    NDIndex<DIM> densDom(Index(_gPDom[0].first() + CellAdd.first(), _gPDom[0].last() + CellAdd.last()),
                         Index(_gPDom[1].first() + CellAdd.first(), _gPDom[1].last() + CellAdd.last()));

    for (unsigned int d = 0; d < DIM; ++ d) {
        spacing[d] = _mesh.get_meshSpacing(d);
        origin[d] += densDom[d].first() * spacing[d];
        densDom[d] = Index(0, densDom[d].length());
    }

    Mesh_t mesh(densDom[0], densDom[1], &spacing[0], origin);
    FieldLayout_Cell_t FL(mesh, decomp);
    SField_Cell_t density(mesh, FL, GuardCellSizes<DIM>(1));

    Timings::stopTimer(_rebalance_timer);
    ShapeFunction<IntCIC>::getNumberDensity(density, *this);
    // BinaryVtkFile vtkFile;
    // vtkFile.addScalarField(density, "density");
    // vtkFile.writeFile("Data/density");
    Timings::startTimer(_rebalance_timer);

    _lPDom = calcRepartition(FL, density);

    if (doRedistribute)
        redistribute(mesh);

    Timings::stopTimer(_rebalance_timer);
}

void PartBunch::redistribute(Mesh_t & mesh)
{
    // This function is more or less copied from IPPLs
    // void ParticleSpatialLayout<...>::update(...)
    // and adapted to the problem at hand.

    const unsigned int numNodes = Ippl::getNodes();
    size_t localNum = getLocalNum();
    size_t destroyNum = getDestroyNum();
    size_t totalNum;

    performDestroy();
    localNum -= destroyNum;

    if (numNodes == 1) {
        setTotalNum(localNum);
        setLocalNum(localNum);
        return;
    }

    localNum = swap_particles(localNum, mesh);
    totalNum = localNum;

    int tag1 = Ippl::Comm->next_tag(P_SPATIAL_LAYOUT_TAG, P_LAYOUT_CYCLE);
    int tag2 = Ippl::Comm->next_tag(P_SPATIAL_RETURN_TAG, P_LAYOUT_CYCLE);
    if (Ippl::myNode() != 0) {
        Message *msg = new Message;

        msg->put(localNum);
        Ippl::Comm->send(msg, 0, tag1);

        int node = 0;
        Message* recmsg = Ippl::Comm->receive_block(node, tag2);
        recmsg->get(totalNum);
        delete recmsg;
    } else {
        int notrecvd = numNodes - 1;
        while (notrecvd > 0) {
            int node = Communicate::COMM_ANY_NODE;
            Message *recmsg = Ippl::Comm->receive_block(node, tag1);
            size_t remNodeCount = 0;
            recmsg->get(remNodeCount);
            delete recmsg;
            notrecvd--;

            totalNum += remNodeCount;
        }

        Message *msg = new Message;
        msg->put(totalNum);
        Ippl::Comm->broadcast_others(msg, tag2);
    }
    setTotalNum(totalNum);
    setLocalNum(localNum);
    balanceParticleNumbers();
}

void PartBunch::getAllLocalPDomains()
{
    Timings::startTimer(_others_timer);
    const unsigned int numNodes = Ippl::getNodes();
    if (numNodes > 1) {
        boost::mpi::communicator world;
        std::vector<int> data(2*DIM,0);
        for (unsigned int i = 0; i < DIM; ++ i) {
            data[2 * i    ] = _lPDom[i].first();
            data[2 * i + 1] = _lPDom[i].last();
        }

        std::vector<std::vector<int> > globalData;
        boost::mpi::all_gather(world, data, globalData);

        for (size_t i = 0; i < numNodes; ++ i) {
            NDIndex<DIM> & dom = _localPDomains[i];
            for (size_t d = 0; d < DIM; ++ d) {
                dom[d] = Index(globalData[i][2 * d],
                               globalData[i][2 * d + 1]);
            }
        }
    } else {
        _localPDomains[0] = _lPDom;
    }
    Timings::stopTimer(_others_timer);
}

size_t PartBunch::swap_particles(size_t localNum,
                                 Mesh_t & mesh)
{
    // This method is more or less copied from IPPLs
    // size_t ParticleSpatialLayout<...>::new_swap_particles(...)
    // and adapted to the problem at hand.

    Ippl::Comm->barrier();

    unsigned numNodes = Ippl::getNodes();
    unsigned myNode = Ippl::myNode();

    int *msgsend = new int[numNodes];
    std::fill(msgsend, msgsend+numNodes, 0);
    int *msgrecv = new int[numNodes];
    std::fill(msgrecv, msgrecv+numNodes, 0);

    NDRegion<double, DIM> pLoc;

    std::multimap<unsigned, unsigned> p2n; //<node ID, particle ID>

    RLayout_t::iterator rDom, remoteEnd = _localPDomains.end();

    for (unsigned int ip=0; ip<localNum; ++ip) {

        NDIndex<DIM> vert = mesh.getCellContaining(R[ip]);

        unsigned destination = myNode;
        if (_lPDom.contains(vert))
            continue;

        for (size_t i = 0; i < _localPDomains.size(); ++ i)
            if (_localPDomains[i].contains(vert))
                destination = i;

        msgsend[destination] = 1;

        p2n.insert(std::pair<unsigned, unsigned>(destination, ip));
    }

    //reduce message count so every node knows how many messages to receive
    DBGOUT << "using collective communication" << std::endl;
    MPI_Allreduce(msgsend, msgrecv, numNodes, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    int tag = Ippl::Comm->next_tag(P_SPATIAL_TRANSFER_TAG,P_LAYOUT_CYCLE);

    typename std::multimap<unsigned, unsigned>::iterator i = p2n.begin();

    Format *format = getFormat();


    std::vector<MPI_Request> requests;
    std::vector<MsgBuffer*> buffers;

    while (i!=p2n.end()) {
        unsigned cur_destination = i->first;

        MsgBuffer *msgbuf = new MsgBuffer(format, p2n.count(i->first));

        for (; i!=p2n.end() && i->first == cur_destination; ++i) {
            Message msg;
            putMessage(msg, i->second);
            destroy(1, i->second);
            msgbuf->add(&msg);
        }

        MPI_Request request = Ippl::Comm->raw_isend( msgbuf->getBuffer(), msgbuf->getSize(), cur_destination, tag);

        //remember request and buffer so we can delete them later
        requests.push_back(request);
        buffers.push_back(msgbuf);
    }

    localNum -= getDestroyNum();
    performDestroy();

    //receive new particles
    for (int k = 0; k<msgrecv[myNode]; ++k) {
        int node = Communicate::COMM_ANY_NODE;
        char *buffer = 0;
        int bufsize = Ippl::Comm->raw_probe_receive(buffer, node, tag);
        MsgBuffer recvbuf(format, buffer, bufsize);

        Message *msg = recvbuf.get();
        while (msg != 0) {
            localNum += getSingleMessage(*msg);
            delete msg;
            msg = recvbuf.get();
        }


    }

    //wait for communication to finish and clean up buffers
    MPI_Waitall(requests.size(), &(requests[0]), MPI_STATUSES_IGNORE);
    for (unsigned int j = 0; j<buffers.size(); ++j) {
        delete buffers[j];
    }

    delete[] msgsend;
    delete[] msgrecv;
    delete format;

    return localNum;
}

void PartBunch::balanceParticleNumbers()
{
    boost::mpi::communicator world;

    unsigned int totalNumParticles = getTotalNum();
    unsigned long localNum = getLocalNum();
    unsigned int  myCol, myRow, ind;
    unsigned int mynode = Ippl::myNode();
    unsigned int numNodes = Ippl::getNodes();
    const auto & node2index = Communicator::getNode2Index();
    const auto & index2node = Communicator::getIndex2Node();
    const auto procMeshSize = Communicator::getNodeMeshSize();
    ind = node2index[mynode];
    myCol = ind % procMeshSize.first;
    myRow = ind / procMeshSize.first;

    size_t averageNumParticles = std::floor(static_cast<double>(totalNumParticles) / Ippl::getNodes());
    size_t missing = totalNumParticles - averageNumParticles * Ippl::getNodes();
    if (static_cast<size_t>(Ippl::myNode()) < missing) {
        averageNumParticles ++;
    }
    std::vector<unsigned long> numLocalParticles(numNodes);
    boost::mpi::all_gather(world, localNum, numLocalParticles);
    std::vector<long> unbalanced(numNodes,0.0);
    for (unsigned int i = 0; i < numNodes; ++ i)
        unbalanced[i] = static_cast<long>(numLocalParticles[i]) - static_cast<long>(averageNumParticles);
    std::vector<std::pair<long,long> > particleMoves(numNodes);
    unsigned int nodeA, nodeB, nodeC;
    for (unsigned int i = 0; i < procMeshSize.first - 1; ++ i) {
        for (unsigned int j = 0; j < procMeshSize.second - 1; ++ j) {
            ind = j * procMeshSize.first + i;
            nodeA = index2node[ind];
            if (unbalanced[nodeA] == 0) continue;
            nodeB = index2node[ind + procMeshSize.first];
            nodeC = index2node[ind + 1];

            long absA = std::abs(unbalanced[nodeA]);
            int signA = unbalanced[nodeA] / absA;
            long B = unbalanced[nodeB];
            long C = unbalanced[nodeC];

            particleMoves[nodeA].first = signA * (std::min(absA, std::max(0L, signA * (C - B)))
                                                  + std::ceil((absA - std::min(absA, std::abs(B - C)))/ 2));
            particleMoves[nodeA].second = signA * (std::min(absA, std::max(0L, signA * (B - C)))
                                                   + std::floor((absA - std::min(absA, std::abs(B - C)))/ 2));
            unbalanced[nodeB] += particleMoves[nodeA].first;
            unbalanced[nodeC] += particleMoves[nodeA].second;
        }
        ind = (procMeshSize.second - 1) * procMeshSize.first + i;
        nodeA = index2node[ind];
        if (unbalanced[nodeA] == 0) continue;
        nodeC = index2node[ind + 1];
        particleMoves[nodeA].first = 0;
        particleMoves[nodeA].second = unbalanced[nodeA];
        unbalanced[nodeC] += particleMoves[nodeA].second;

    }
    for (unsigned int j = 0; j < procMeshSize.second - 1; ++ j) {
        ind = (j + 1) * procMeshSize.first - 1;
        nodeA = index2node[ind];
        nodeB = index2node[ind + procMeshSize.first];
        particleMoves[nodeA].first = unbalanced[nodeA];
        particleMoves[nodeA].second = 0;
        unbalanced[nodeB] += particleMoves[nodeA].first;
    }

    int tag = Ippl::Comm->next_tag(P_SPATIAL_TRANSFER_TAG,P_LAYOUT_CYCLE);
    Format *format = getFormat();
    std::vector<MPI_Request> requests;
    std::vector<MsgBuffer*> buffers;

    unsigned int destination;
    if (particleMoves[mynode].first > 0) {
        ind = (myRow + 1) * procMeshSize.first + myCol;
        destination = index2node[ind];
        std::vector<std::pair<unsigned long, double> > ycomps(localNum); // <particle ID, y component of R>
        for (unsigned long ip = 0; ip < localNum; ++ ip)
            ycomps[ip] = std::make_pair(ip, R[ip](1));
        std::sort(ycomps.begin(), ycomps.end(), descSecond);

        MsgBuffer *msgbuf = new MsgBuffer(format, particleMoves[mynode].first);
        for (unsigned int i = 0; i < particleMoves[mynode].first; ++ i) {
            Message msg;
            unsigned long ip = ycomps[i].first;
            putMessage(msg, ip);
            destroy(1, ip);
            msgbuf->add(&msg);
        }

        MPI_Request request = Ippl::Comm->raw_isend( msgbuf->getBuffer(), msgbuf->getSize(), destination, tag);
        requests.push_back(request);
        buffers.push_back(msgbuf);
    }

    localNum -= getDestroyNum();
    performDestroy();

    if (myRow > 0) {
        ind = (myRow - 1) * procMeshSize.first + myCol;
        destination = index2node[ind];
        if (particleMoves[destination].first < 0) {
            std::vector<std::pair<unsigned long, double> > ycomps(localNum); // <particle ID, y component of R>
            for (unsigned long ip = 0; ip < localNum; ++ ip)
                ycomps[ip] = std::make_pair(ip, R[ip](1));
            std::sort(ycomps.begin(), ycomps.end(), ascSecond);

            MsgBuffer *msgbuf = new MsgBuffer(format, -particleMoves[destination].first);
            for (unsigned int i = 0; i < -particleMoves[destination].first; ++ i) {
                Message msg;
                unsigned long ip = ycomps[i].first;
                putMessage(msg, ip);
                destroy(1, ip);
                msgbuf->add(&msg);
            }

            MPI_Request request = Ippl::Comm->raw_isend( msgbuf->getBuffer(), msgbuf->getSize(), destination, tag);
            requests.push_back(request);
            buffers.push_back(msgbuf);
        }
    }

    localNum -= getDestroyNum();
    performDestroy();

    if (particleMoves[mynode].second > 0) {
        ind = myRow * procMeshSize.first + myCol + 1;
        destination = index2node[ind];
        std::vector<std::pair<unsigned long, double> > xcomps(localNum); // <particle ID, y component of R>
        for (unsigned long ip = 0; ip < localNum; ++ ip)
            xcomps[ip] = std::make_pair(ip, R[ip](0));
        std::sort(xcomps.begin(), xcomps.end(), descSecond);

        MsgBuffer *msgbuf = new MsgBuffer(format, particleMoves[mynode].second);
        for (unsigned int i = 0; i < particleMoves[mynode].second; ++ i) {
            Message msg;
            unsigned long ip = xcomps[i].first;
            putMessage(msg, ip);
            destroy(1, ip);
            msgbuf->add(&msg);
        }

        MPI_Request request = Ippl::Comm->raw_isend( msgbuf->getBuffer(), msgbuf->getSize(), destination, tag);
        requests.push_back(request);
        buffers.push_back(msgbuf);
    }

    localNum -= getDestroyNum();
    performDestroy();

    if (myCol > 0) {
        ind = myRow * procMeshSize.first + myCol - 1;
        destination = index2node[ind];
        if (particleMoves[destination].second < 0) {
            std::vector<std::pair<unsigned long, double> > xcomps(localNum); // <particle ID, y component of R>
            for (unsigned long ip = 0; ip < localNum; ++ ip)
                xcomps[ip] = std::make_pair(ip, R[ip](0));
            std::sort(xcomps.begin(), xcomps.end(), ascSecond);


            MsgBuffer *msgbuf = new MsgBuffer(format, -particleMoves[destination].second);
            for (unsigned int i = 0; i < -particleMoves[destination].second; ++ i) {
                Message msg;
                unsigned long ip = xcomps[i].first;
                putMessage(msg, ip);
                destroy(1, ip);
                msgbuf->add(&msg);
            }

            MPI_Request request = Ippl::Comm->raw_isend( msgbuf->getBuffer(), msgbuf->getSize(), destination, tag);
            requests.push_back(request);
            buffers.push_back(msgbuf);
        }
    }

    localNum -= getDestroyNum();
    performDestroy();

    unsigned numRecv = 0;
    if (particleMoves[mynode].first < 0) ++ numRecv;
    if (particleMoves[mynode].second < 0) ++ numRecv;
    if (myCol > 0) {
        ind = myRow * procMeshSize.first + myCol - 1;
        destination = index2node[ind];
        if (particleMoves[destination].second > 0) ++ numRecv;
    }
    if (myRow > 0) {
        ind = (myRow - 1) * procMeshSize.first + myCol;
        destination = index2node[ind];
        if (particleMoves[destination].first > 0) ++ numRecv;
    }

    //receive new particles
    for (unsigned int k = 0; k < numRecv; ++ k) {
        int node = Communicate::COMM_ANY_NODE;
        char *buffer = 0;

        int bufsize = Ippl::Comm->raw_probe_receive(buffer, node, tag);
        MsgBuffer recvbuf(format, buffer, bufsize);

        Message *msg = recvbuf.get();
        while (msg != 0) {
            localNum += getSingleMessage(*msg);
            delete msg;
            msg = recvbuf.get();
        }
    }

    //wait for communication to finish and clean up buffers
    MPI_Waitall(requests.size(), &(requests[0]), MPI_STATUSES_IGNORE);
    for (unsigned int j = 0; j<buffers.size(); ++j) {
        delete buffers[j];
    }

    delete format;

    setLocalNum(localNum);
}

void PartBunch::sortParticles()
{
    double dy = _mesh.get_meshSpacing(1);
    Vector_t origin = _mesh.get_origin();
    NDIndex<DIM> lPDom = getLocalPDomain(_mesh);
    _others_timer = Timings::getTimer("other stuff");
    Timings::startTimer(_others_timer);
    unsigned long localNum = getLocalNum();
    std::vector<std::vector<Particle> > yslices(lPDom[1].length());

    for (unsigned long ip = 0; ip < localNum; ++ ip) {
        unsigned int idy = static_cast<int>(std::floor((R[ip](1) - origin(1)) / dy)) - lPDom[1].first();
        yslices[idy].push_back(Particle(oldR[ip], R[ip], P[ip], Q[ip]));
    }

    unsigned long ip = 0;
    for (unsigned int idy = 0; idy < yslices.size(); ++ idy) {
        std::sort(yslices[idy].begin(), yslices[idy].end(),
                  [](const Particle & a, const Particle & b) { return a.getR()(0) < b.getR()(0); });
        for (auto it = yslices[idy].begin(); it != yslices[idy].end(); ++ it, ++ ip) {
            R[ip] = it->getR();
            oldR[ip] = it->getOldR();
            P[ip] = it->getP();
            Q[ip] = it->getQ();
        }
    }
    Timings::stopTimer(_others_timer);
}

NDIndex<DIM> PartBunch::calcRepartition(const FieldLayout_Cell_t & FL, const SField_Cell_t & rho)
{
    if (Ippl::getNodes() == 1) return FL.getDomain();

    std::vector<NDIndex<DIM> > vectorLocalDomains;
    Utils::getLocalDomains(FL, vectorLocalDomains);
    const auto & index2node = Communicator::getIndex2Node();
    const auto procMeshSize = Communicator::getNodeMeshSize();

    boost::mpi::communicator world;
    NDIndex<DIM> elem;
    NDIndex<DIM> gDom = FL.getDomain();
    NDIndex<DIM> & lDom = vectorLocalDomains[Ippl::myNode()];
    std::vector<double> localLineDensityX(gDom[0].length(), 0.0);
    std::vector<double> localLineDensityY(gDom[1].length(), 0.0);
    std::vector<double> lineDensityX(gDom[0].length(), 0.0);
    std::vector<double> lineDensityY(gDom[1].length(), 0.0);

    for (int j = lDom[1].first(); j <= lDom[1].last(); ++ j) {
        elem[1] = Index(j,j);
        for (int i = lDom[0].first(); i <= lDom[0].last(); ++ i) {
            elem[0] = Index(i,i);
            localLineDensityX[i] += rho.localElement(elem);
            localLineDensityY[j] += rho.localElement(elem);
        }
    }
    DBGOUT << "using collective communication" << std::endl;
    boost::mpi::all_reduce(world, &localLineDensityX[0], gDom[0].length(), &lineDensityX[0], std::plus<double>());
    DBGOUT << "using collective communication" << std::endl;
    boost::mpi::all_reduce(world, &localLineDensityY[0], gDom[1].length(), &lineDensityY[0], std::plus<double>());

    std::partial_sum(lineDensityX.begin(), lineDensityX.end(), localLineDensityX.begin());
    std::partial_sum(lineDensityY.begin(), lineDensityY.end(), localLineDensityY.begin());

    std::vector<unsigned int> cutsX(procMeshSize.first + 1,0);
    std::vector<unsigned int> cutsY(procMeshSize.second + 1,0);

    if (false) {
        unsigned int startX = 0;
        unsigned int lastX = procMeshSize.first;
        unsigned int headX = 0;
        unsigned int maxLengthX = gDom[0].length() * 1.4142 / procMeshSize.first;

        cutsX[startX] = gDom[0].first() - 1;
        cutsX[lastX] = gDom[0].last();

        while (2 * (startX + 1) < procMeshSize.first || procMeshSize.second <= 2) {
            unsigned int firstIdx = cutsX[startX] + 1;
            unsigned int lastIdx = cutsX[lastX];
            unsigned int cutnr = 1 + startX;
            double stepSize = 1.0 / (lastX - startX) * (localLineDensityX[lastIdx - 1] - headX);
            double cutat = stepSize + headX;

            for (unsigned int i = firstIdx; i < lastIdx; ++ i) {
                if (localLineDensityX[i] > cutat) {
                    if (localLineDensityX[i] - cutat < cutat - localLineDensityX[i-1]) {
                        cutsX[cutnr] = i;
                    } else {
                        cutsX[cutnr] = i-1;
                    }
                    ++ cutnr;
                    cutat += stepSize;
                }
            }

            bool cont = false;
            if (cutsX[startX + 1] - cutsX[startX] > maxLengthX) {
                cutsX[startX + 1] = cutsX[startX] + maxLengthX;
                cont = true;
            }

            if (cutsX[lastX] - cutsX[lastX - 1] > maxLengthX) {
                cutsX[lastX - 1] = cutsX[lastX] - maxLengthX;
                cont = true;
            }

            headX = localLineDensityX[cutsX[startX + 1]];
            ++ startX;
            --lastX;

            if (cont) continue;

            break;
        }

        unsigned int startY = 0;
        unsigned int lastY = procMeshSize.second;
        unsigned int headY = 0;
        unsigned int maxLengthY = gDom[1].length() * 1.4142 / procMeshSize.second;

        cutsY[startY] = gDom[1].first() - 1;
        cutsY[lastY] = gDom[1].last();

        while (2 * (startY + 1) < procMeshSize.second || procMeshSize.second <= 2) {
            unsigned int firstIdx = cutsY[startY] + 1;
            unsigned int lastIdx = cutsY[lastY];
            unsigned int cutnr = 1 + startY;
            double stepSize = 1.0 / (lastY - startY) * (localLineDensityY[lastIdx - 1] - headY);
            double cutat = stepSize + headY;

            for (unsigned int i = firstIdx; i < lastIdx; ++ i) {
                if (localLineDensityY[i] > cutat) {
                    if (localLineDensityY[i] - cutat < cutat - localLineDensityY[i-1]) {
                        cutsY[cutnr] = i;
                    } else {
                        cutsY[cutnr] = i-1;
                    }
                    ++ cutnr;
                    cutat += stepSize;
                }
            }

            bool cont = false;
            if (cutsY[startY + 1] - cutsY[startY] > maxLengthY) {
                cutsY[startY + 1] = cutsY[startY] + maxLengthY;
                cont = true;
            }

            if (cutsY[lastY] - cutsY[lastY - 1] > maxLengthY) {
                cutsY[lastY - 1] = cutsY[lastY] - maxLengthY;
                cont = true;
            }

            headY = localLineDensityY[cutsY[startY + 1]];
            ++ startY;
            --lastY;

            if (cont) continue;

            break;
        }
    } else {
        unsigned int cutnr = 1;
        double cutat = 1.0 / procMeshSize.first * localLineDensityX.back();
        for (unsigned int i = 0; i < lineDensityX.size(); ++ i) {
            if (localLineDensityX[i] > cutat) {
                if (localLineDensityX[i] - cutat < cutat - localLineDensityX[i-1]) {
                    cutsX[cutnr] = i;
                } else {
                    cutsX[cutnr] = i-1;
                }
                ++ cutnr;
                cutat = cutnr * 1.0 / procMeshSize.first * localLineDensityX.back();
            }
        }
        cutsX[0] = gDom[0].first() - 1;
        cutsX[procMeshSize.first] = gDom[0].last();

        cutnr = 1;
        cutat = 1.0 / procMeshSize.second * localLineDensityY.back();
        for (unsigned int i = 0; i < lineDensityY.size(); ++ i) {
            if (localLineDensityY[i] > cutat) {
                if (localLineDensityY[i] - cutat < cutat - localLineDensityY[i-1]) {
                    cutsY[cutnr] = i;
                } else {
                    cutsY[cutnr] = i-1;
                }
                ++ cutnr;
                cutat = cutnr * 1.0 / procMeshSize.second * localLineDensityY.back();
            }
        }
        cutsY[0] = gDom[1].first() - 1;
        cutsY[procMeshSize.second] = gDom[1].last();
    }

    for (unsigned int i = 0; i < _localPDomains.size(); ++ i) {
        auto it = std::find(index2node.begin(), index2node.end(), i);
        unsigned int ind = it - index2node.begin();
        unsigned int col = ind % procMeshSize.first;
        unsigned int row = ind / procMeshSize.first;
        _localPDomains[i][0] = Index(cutsX[col] + 1, cutsX[col+1]);
        _localPDomains[i][1] = Index(cutsY[row] + 1, cutsY[row+1]);
    }

    return _localPDomains[Ippl::myNode()];
}

struct MasterSlaveRelation {
    unsigned int _master;
    std::vector<int> _slaves;
    double _quorum;

    MasterSlaveRelation(unsigned int m, double n):
        _master(m),
        _quorum(n)
    { }
};

void PartBunch::balanceParticles()
{
    boost::mpi::communicator world;
    unsigned int myNode = Ippl::myNode();
    unsigned int numNodes = Ippl::getNodes();
    double numParticlesPerCore = 1.0 * getTotalNum() / Ippl::getNodes();
    std::vector<MasterSlaveRelation> rel;
    std::vector<std::pair<unsigned int, size_t> > slaves;

    for (unsigned int k = 0; k < numNodes; ++ k) {
        if (_localNums[k] > numParticlesPerCore) {
            rel.push_back(MasterSlaveRelation(k, 100 * _localNums[k] / numParticlesPerCore));
        } else if (_localNums[k] < numParticlesPerCore) {
            if (k == myNode) _iAmSlave = true;
            slaves.push_back(std::make_pair(k, _localNums[k]));
        }
    }

    std::sort(slaves.begin(), slaves.end(),
              [](const std::pair<unsigned int, size_t> & a, const std::pair<unsigned int, size_t> & b) -> bool {
                  return a.second > b.second;
              });

    while (slaves.size() > 0) {
        std::sort(rel.begin(), rel.end(),
                  [](const MasterSlaveRelation& a, const MasterSlaveRelation& b) ->bool {
                      return a._quorum > b._quorum;
                  });

        if (rel[0]._master == myNode) {
            _iAmMaster = true;
            _commRoot = rel[0]._master;
        }
        if (slaves.back().first == myNode) {
            _commRoot = rel[0]._master;
        }
        rel[0]._slaves.push_back(slaves.back().first);
        rel[0]._quorum = (rel[0]._quorum * rel[0]._slaves.size() +
                          100 * slaves.back().second / numParticlesPerCore) / (rel[0]._slaves.size() + 1);
        slaves.pop_back();
    }

    unsigned int myGroup = -1, test = -1;
    for (unsigned int i = 0; i < rel.size(); ++ i) {
        if (myNode == rel[i]._master) {
            myGroup = i;
            break;
        }
        auto it = find(rel[i]._slaves.begin(), rel[i]._slaves.end(), myNode);
        if (it != rel[i]._slaves.end()) {
            myGroup = i;
            break;
        }
    }
    if (myGroup == test) {
        if (Ippl::getNodes() > 1)
            DBGOUT << "something went wrong!" << std::endl;
    } else {
        rel[myGroup]._slaves.push_back(rel[myGroup]._master);
        MPI_Group worldgroup, newgroup;
        MPI_Comm_group(MPI_COMM_WORLD, &worldgroup);
        MPI_Group_incl(worldgroup, rel[myGroup]._slaves.size(), &(rel[myGroup]._slaves[0]), &newgroup);
        MPI_Comm_create(MPI_COMM_WORLD, newgroup, &_commMasterSlave);
        _commRoot = rel[myGroup]._slaves.size() - 1;
        // int myRank = 0;
        // if (_iAmMaster) {
        //     MPI_Comm_rank(_commMasterSlave, &myRank);
        // }
        // MPI_Allreduce(&myRank, &_commRoot, 1, MPI_INT, MPI_SUM, _commMasterSlave);
        // MPI_Comm_rank(_commMasterSlave, &myRank);
        // unsigned long localNum = _localNums[myNode];
        // std::vector<unsigned long> localNumsGroup(rel[myGroup]._slaves.size(),0);
        // MPI_Gather(&localNum, 1, MPI_UNSIGNED_LONG, &(localNumsGroup[0]), 1, MPI_UNSIGNED_LONG, _commRoot, _commMasterSlave);

        size_t sumParticles = 0;
        for (const unsigned int k: rel[myGroup]._slaves) {
            sumParticles += _localNums[k];
        }
        int tag = Ippl::Comm->next_tag(IPPL_APP_TAG0, IPPL_APP_CYCLE);
        numParticlesPerCore = 1.0 * sumParticles / rel[myGroup]._slaves.size();
        unsigned int missing = sumParticles - std::floor(numParticlesPerCore) * rel[myGroup]._slaves.size();

        for (unsigned k = 0; k < rel[myGroup]._slaves.size() - 1; ++ k) {
            unsigned int worldrank = rel[myGroup]._slaves[k];
            size_t numParticles = std::floor(numParticlesPerCore) + (missing > k? 1: 0) - _localNums[worldrank];
            if (worldrank == myNode) {
                std::vector<double> receiveParticles(numParticles * (3 * DIM + 1));
                world.recv(rel[myGroup]._master, tag, &(receiveParticles[0]), numParticles * (3 * DIM + 1));
                _masterStart = _localNums[worldrank];
                _masterEnd = _localNums[worldrank] + numParticles;
                _masterWorldRank = rel[myGroup]._master;

                create(numParticles);
                size_t l = 0;
                for (size_t lnp = _masterStart; lnp < _masterEnd; ++ lnp) {
                    for (unsigned int d = 0; d < DIM; ++ d)
                        R[lnp](d) = receiveParticles[l++];
                    for (unsigned int d = 0; d < DIM; ++ d)
                        oldR[lnp](d) = receiveParticles[l++];
                    for (unsigned int d = 0; d < DIM; ++ d)
                        P[lnp](d) = receiveParticles[l++];
                    Q[lnp] = receiveParticles[l++];
                }

            } else if (_iAmMaster) {
                std::vector<double> sendParticles(numParticles * (3 * DIM + 1));
                size_t l = 0;
                for (size_t lnp = getLocalNum() - 1; lnp >= getLocalNum() - numParticles; -- lnp) {
                    for (unsigned int d = 0; d < DIM; ++ d)
                        sendParticles[l++] = R[lnp](d);
                    for (unsigned int d = 0; d < DIM; ++ d)
                        sendParticles[l++] = oldR[lnp](d);
                    for (unsigned int d = 0; d < DIM; ++ d)
                        sendParticles[l++] = P[lnp](d);
                    sendParticles[l++] = Q[lnp];
                }
                world.send(worldrank, tag, &(sendParticles[0]), numParticles * (3 * DIM + 1));

                destroy(numParticles, getLocalNum() - numParticles, true);
            }
        }

        if (!_iAmSlave) {
            _masterStart = getLocalNum();
            _masterEnd = getLocalNum();
        }
    }
    for (const auto & r: rel) {
        DBGOUT << r._master << "\t" << "\t" << r._quorum << std::endl;
        for (unsigned int i: r._slaves)
            DBGOUT << i << std::endl;
        dbg << std::endl;
    }
}
