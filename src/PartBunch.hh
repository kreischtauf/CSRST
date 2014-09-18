/***************************************************************************
                             PartBunch.hh
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

#ifndef OPAL_PartBunch_HH
#define OPAL_PartBunch_HH

#include <fstream>
#include <string>
#include <vector>
#include "Physics.hh"
#include "Bunch.hh"
#include "ShapeFunction.hh"

class Bend;
class Distribution;
template<class T>
class FieldPatch;

class PartBunchState;

class PartBunch: public Bunch {

public:
    PartBunch(FieldLayout<DIM> & edgeFL,
              FieldLayout<DIM> & cellFL,
              Mesh_t & mesh,
              Distribution & distro,
              const Bend & bend,
              const double& charge_mass_ratio);

    /*virtual*/
    void makeReady(const double & dt);

    void commNumLocalParticles() const;

    /*virtual*/
    void push_particles(VField_Edge_t& EFD,
                        VField_Cell_t& BFD,
                        VField_Edge_t& JFD,
                        const double & dt,
                        const unsigned int & Jcomp);

    /*virtual*/
    void create_particle(const Vector_t & X,
                         const Vector_t & PX,
                         const double & q);

    void calcBeamParameters(PartBunchState & state) const;

    /*virtual*/
    void get_rmean(Vector_t & spos) const;

    /*virtual*/
    void get_pmean(Vector_t & pmean) const {
        pmean = _momentum;
    }

    /*virtual*/
    double get_qtotal() const;

    /*virtual*/
    size_t getLocalNP() const {
        return getLocalNum();
    }

    void getBounds(Vector_t & Rmin, Vector_t & Rmax, size_t from, size_t to) const;
    void getBoundsInclOld(Vector_t & Rmin, Vector_t & Rmax, size_t from, size_t to) const;

    const std::vector<size_t> & getNumLocalParticles() const;
    NDIndex<DIM> getLocalPDomain(const Mesh_t & mesh) const;
    NDIndex<DIM> getLocalPDomainInclOld(const Mesh_t & mesh) const;
    NDIndex<DIM> getGlobalPDomain(const Mesh_t & mesh) const;
    NDIndex<DIM> getGlobalPDomainInclOld(const Mesh_t & mesh) const;

    /*virtual*/
    void scatterQ(SField_t & rho) const;

    /*virtual*/
    void drift_particles(const double & dt);

    /*virtual*/
    void move_by(const Vector_t & dX);

    /*virtual*/
    void rotate(const Vector_t & X,
                const double& phi);

    void rebalance(bool doRedistribute = true);

    void setGlobalBounds(const NDIndex<DIM> & dom) {
        _globalBounds = dom;
    }

    NDIndex<DIM> getGlobalBounds() const {
        return _globalBounds;
    }

    void balanceParticles();

    void saveR();
    void restoreR();
private:
    void getDistribution(Distribution & distro, const Vector_t & dx);

    void calcChargeConservation(const SField_Vert_t & rho,
                                const VField_Edge_t & JFD,
                                const double & dt,
                                const int & iteration);

    void save(const VField_Edge_t & field, const std::string & fname, const int comp);

    void getFields(FieldPatch<double> & Ex,
                   FieldPatch<double> & Ey,
                   FieldPatch<double> & Hz,
                   const VField_Edge_t& EFD,
                   const VField_Cell_t& HFD);

    void getFieldsNEW(FieldPatch<double> & Ex,
                      FieldPatch<double> & Ey,
                      FieldPatch<double> & Hz,
                      const VField_Edge_t& EFD,
                      const VField_Cell_t& HFD);

    void communicateFields(std::vector<FieldPatch<double> * > localFields,
                           const std::vector<std::vector<NDIndex<DIM> > > & allLocalFDomains,
                           const std::vector<std::vector<NDIndex<DIM> > > & allLocalPDomains) const;

    typedef std::vector<NDIndex<DIM> > RLayout_t;
    void getAllLocalPDomains();

    void redistribute(Mesh_t & mesh);
    NDIndex<DIM> calcRepartition(const FieldLayout_Cell_t & FL, const SField_Cell_t & rho);
    size_t swap_particles(size_t localNum, Mesh_t & mesh);
    void balanceParticleNumbers();
    void sortParticles();

    void calcMoments(PartBunchState & state) const;

    ParticleAttrib< Vector_t > Ef;
    ParticleAttrib< double > Hf;
    ParticleAttrib< double > M;

    NDIndex<DIM> _lPDom;
    NDIndex<DIM> _gPDom;
    RLayout_t _localPDomains;
    RLayout_t _localFDomains;

    playout_t * _pl;

    ShapeFunction<IntOp_t> _shapeFunction;

    /// holds the timestep in seconds
    double _dt;

    double _charge;
    double _charge_mass_ratio;
    Vector_t _momentum;
    NDIndex<DIM> _globalBounds;

    std::vector<size_t> _localNums;
    bool _localNumsSet;

    const Bend & _bend;
    const Mesh_t & _mesh;

    size_t _masterStart;
    size_t _masterEnd;
    int _commRoot;
    int _masterWorldRank;
    MPI_Comm _commMasterSlave;
    bool _iAmSlave;
    bool _iAmMaster;

    Timings::TimerRef _push_timer;
    Timings::TimerRef _conservation_timer;
    Timings::TimerRef _rebalance_timer;
    Timings::TimerRef _fieldInterp_timer;
    Timings::TimerRef _comm_timer;
    Timings::TimerRef _bounds_timer;
    Timings::TimerRef _others_timer;
    Timings::TimerRef _statParamTimer;
public:
    static bool Vector_0_sort_ASC(const Vector_t& v1, const Vector_t& v2)
    {
        return (v1(0) < v2(0));
    }

    static bool Vector_0_sort_DESC(const Vector_t& v1, const Vector_t& v2)
    {
        return (v1(0) > v2(0));
    }
};

inline
const std::vector<size_t> & PartBunch::getNumLocalParticles() const
{
    return _localNums;
}

inline
void PartBunch::create_particle(const Vector_t & X,
                                const Vector_t & PX,
                                const double & q)
{
    size_t nPart = this->getLocalNum();
    create(1);
    this->R[nPart] = X;
    this->P[nPart] = PX;
    this->Q[nPart] = q;
    this->Ef[nPart] = Vector_t(0.0);
    this->Hf[nPart] = 0.0;
}

inline
void PartBunch::scatterQ(SField_t & rho) const
{
    _shapeFunction.getChargeDensity(rho, *this);
}

inline
void PartBunch::saveR()
{
    this->Ef = this->R;
}

inline
void PartBunch::restoreR()
{
    this->R = this->Ef;
    this->Ef = 0.0;
}

inline
void PartBunch::getBounds(Vector_t & Rmin, Vector_t & Rmax, size_t from, size_t to) const
{
    to = std::min(to, getLocalNum());
    from = std::min(from, to);
    if (from < to) {
        Rmin = R[from];
        Rmax = R[from];
        for (size_t i = from; i < to; ++ i) {
            for (unsigned int d = 0; d < DIM; ++ d) {
                if (R[i](d) < Rmin(d)) Rmin(d) = R[i](d);
                if (R[i](d) > Rmax(d)) Rmax(d) = R[i](d);
            }
        }
    } else {
        Rmin = Vector_t(0.0);
        Rmax = Vector_t(0.0);
    }
}

inline
void PartBunch::getBoundsInclOld(Vector_t & Rmin, Vector_t & Rmax, size_t from, size_t to) const
{
    to = std::min(to, getLocalNum());
    from = std::min(from, to);
    if (from < to) {
        Rmin = R[from];
        Rmax = R[from];
        for (size_t i = from; i < to; ++ i) {
            for (unsigned int d = 0; d < DIM; ++ d) {
                if (R[i](d) < Rmin(d)) Rmin(d) = R[i](d);
                if (oldR[i](d) < Rmin(d)) Rmin(d) = oldR[i](d);
                if (R[i](d) > Rmax(d)) Rmax(d) = R[i](d);
                if (oldR[i](d) > Rmax(d)) Rmax(d) = oldR[i](d);
            }
        }
    } else {
        Rmin = Vector_t(0.0);
        Rmax = Vector_t(0.0);
    }
}


class PartBunchState {
public:

    const double & get_meanEnergy() const;
    const Vector_t & get_rmin() const;
    const Vector_t & get_rmax() const;
    const Vector_t & get_rrms() const;
    const Vector_t & get_prms() const;
    const Vector_t & get_rmean() const;
    const Vector_t & get_pmean() const;
    const Vector_t & get_norm_emit() const;
    const double & get_Dy() const;
    const double & get_DDy() const;
    const double & get_dE() const;

private:
    friend class PartBunch;

    Vector_t _rmin;
    Vector_t _rmax;
    Vector_t _rmean;
    Vector_t _pmean;
    Vector_t _rrms;
    Vector_t _prms;
    Vector_t _epsNorm;
    SymTenzor<long double, 2 * DIM> _moments;
    long double _centroid[2 * DIM];
    double _Ekin;
    Vector_t _rprms;

    double _Dy;
    double _DDy;
    double _dE;
};

inline
const double & PartBunchState::get_meanEnergy() const
{
    return _Ekin;
}

inline
const Vector_t & PartBunchState::get_rmin() const
{
    return _rmin;
}

inline
const Vector_t & PartBunchState::get_rmax() const
{
    return _rmax;
}

inline
const Vector_t & PartBunchState::get_rrms() const
{
    return _rrms;
}

inline
const Vector_t & PartBunchState::get_prms() const
{
    return _prms;
}

inline
const Vector_t & PartBunchState::get_rmean() const
{
    return _rmean;
}

inline
const Vector_t & PartBunchState::get_pmean() const
{
    return _pmean;
}

inline
const Vector_t & PartBunchState::get_norm_emit() const
{
    return _epsNorm;
}

inline
const double & PartBunchState::get_Dy() const
{
    return _Dy;
}

inline
const double & PartBunchState::get_DDy() const
{
    return _DDy;
}

inline
const double & PartBunchState::get_dE() const
{
    return _dE;
}
#endif // PARTBUNCH_HH
