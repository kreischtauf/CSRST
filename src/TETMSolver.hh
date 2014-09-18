/***************************************************************************
                            TETMSolver.hh
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

#ifndef TETMSOLVER_HH
#define TETMSOLVER_HH

#include <vector>

#include "defs.hh"
#include "PartBunch.hh"
#include "BoundaryCell.hh"
#include "Physics.hh"
#include "Communicator.hh"

/** \brief implements the TETM solver
 *   <a href="http://dx.doi.org/10.1016/j.jcp.2005.01.003">TE/TM scheme for computation of electromagnetic fields in accelerators</a>
 *   by I. Zagorodnov and T. Weiland and adapts it to
 *   time dependent problems.
 */

void outputNodeBoundary(VField_Edge_t & fd);
void outputNodeBoundary(VField_Cell_t & fd, VField_Edge_t & fd2, const double & etax, const double & etay);

class TETMSolver {
public:
    TETMSolver(VField_t & AlphaE,
               VField_t & AlphaH,
               PartBunch & bunch,
               Inform & msg,
               const VField_Edge_t & EFD,
               const SolverPrecision_t & dt);
    ~TETMSolver();

    double getDT() const;

    template<class G, class BC>
    void updateField_firstStep(VField_Edge_t & EFD,
                               VField_Cell_t & Hz,
                               VField_Edge_t & JFD,
                               VField_Edge_t & mmE,
                               VField_Cell_t & mmH,
                               VField_Edge_t & mmJ,
                               NDIndex<DIM> & ddim,
                               const BC& bc,
                               G& geo);

    template<class G, class BC>
    void updateField_secondStep(VField_Edge_t& EFD,
                                VField_Cell_t & Hz,
                                VField_Edge_t & JFD,
                                VField_Edge_t & mmE,
                                VField_Cell_t & mmH,
                                VField_Edge_t & mmJ,
                                NDIndex<DIM> & ddim,
                                const BC& bc,
                                G& geo);

    void initMaps(const NDIndex<DIM> & lDom);

    Vector_t updateHz(const VField_Edge_t& EFD,
                      const NDIndex<DIM>& elem,
                      const BoundaryCell bcell) const;

    double updateEx(const VField_Cell_t& HFD,
                    const VField_Edge_t& JFD,
                    const NDIndex<DIM>& elem) const;

    double updateEy(const VField_Cell_t& HFD,
                    const VField_Edge_t& JFD,
                    const NDIndex<DIM>& elem) const;
private:
    PartBunch & _pbunch;
    VField_t & _alpha_e;
    VField_t & _alpha_h;

    NDIndex<DIM> _gDom;

    Idx_t _Nx;
    Idx_t _Ny;

    Idx_t _nx;
    Idx_t _ny;

    double _dt;

    Timings::TimerRef _fieldUpdateTimer;
    Timings::TimerRef _exchangeTimer;
    const int _me;
};

inline
double TETMSolver::getDT() const {
    return _dt;
}

template<class G, class BC>
void TETMSolver::updateField_firstStep(VField_Edge_t & EFD,
                                       VField_Cell_t & HFD,
                                       VField_Edge_t & JFD,
                                       VField_Edge_t & mmE,
                                       VField_Cell_t & mmH,
                                       VField_Edge_t & mmJ,
                                       NDIndex<DIM> & ddim,
                                       const BC& bc,
                                       G& geo)
{
    Vector_t hr(EFD.get_mesh().get_meshSpacing(0), EFD.get_mesh().get_meshSpacing(1));
    double etax = Physics::c * _dt / hr(0);
    double etay = Physics::c * _dt / hr(1);
    double Z_0 = sqrt(Physics::mu_0 / Physics::epsilon_0);
    double dteps = Physics::c * _dt * Z_0;
    Field<bool, DIM> & insideMask = geo.getInsideMask();

    const NDIndex<DIM> & lDom = EFD.getLayout().getLocalNDIndex();
    const NDIndex<DIM> & mmlDomE = mmE.getLayout().getLocalNDIndex();
    const NDIndex<DIM> & mmlDomH = mmH.getLayout().getLocalNDIndex();

    Index II = lDom[0];
    Index JJ = lDom[1];
    Index mmIIe = mmlDomE[0];
    Index mmJJe = mmlDomE[1];
    Index mmIIh = mmlDomH[0];
    Index mmJJh = mmlDomH[1];

    _pbunch.push_particles(mmE, mmH, mmJ, _dt, 0);

    Timings::startTimer(_fieldUpdateTimer);
    assign(HFD[II][JJ](0), where(insideMask[II][JJ],
                                 HFD[II][JJ](0) - 0.5 * etax / Z_0 * (EFD[II+1][JJ](1) - EFD[II][JJ](1)),
                                 HFD[II][JJ](0)));
    assign(HFD[II][JJ](1), where(insideMask[II][JJ],
                                 HFD[II][JJ](1) + 0.5 * etay / Z_0 * (EFD[II][JJ+1](0) - EFD[II][JJ](0)),
                                 HFD[II][JJ](1)));
    assign(mmH[mmIIh][mmJJh](0), mmH[mmIIh][mmJJh](0) - 0.5 * etax / Z_0 * (mmE[mmIIh+1][mmJJh](1) - mmE[mmIIh][mmJJh](1)));
    assign(mmH[mmIIh][mmJJh](1), mmH[mmIIh][mmJJh](1) + 0.5 * etay / Z_0 * (mmE[mmIIh][mmJJh+1](0) - mmE[mmIIh][mmJJh](0)));

    Timings::stopTimer(_fieldUpdateTimer);
    Timings::startTimer(_exchangeTimer);

    Communicator::exchangeBoundaries(HFD,
                                     mmH,
                                     Vektor<int,DIM>(ddim[0].min(), ddim[1].min()),
                                     HFieldToDouble());

    Timings::stopTimer(_exchangeTimer);
    Timings::startTimer(_fieldUpdateTimer);

    bc.updateHz(HFD, EFD, *this, geo);
    HFD.fillGuardCells();

    assign(EFD[II][JJ](0), where(insideMask[II][JJ],
                                 EFD[II][JJ](0) + (etay * Z_0 * (HFD[II][JJ](0) - HFD[II][JJ-1](0) +
                                                                 HFD[II][JJ](1) - HFD[II][JJ-1](1)) -
                                                   dteps * JFD[II][JJ](0)),
                                 EFD[II][JJ](0)));
    assign(mmE[mmIIe][mmJJe](0), mmE[mmIIe][mmJJe](0) + (etay * Z_0 * (mmH[mmIIe][mmJJe](0) - mmH[mmIIe][mmJJe-1](0) +
                                                                       mmH[mmIIe][mmJJe](1) - mmH[mmIIe][mmJJe-1](1)) -
                                                         dteps * mmJ[mmIIe][mmJJe](0)));
    Timings::stopTimer(_fieldUpdateTimer);
    Timings::startTimer(_exchangeTimer);

    Communicator::exchangeBoundaries(EFD,
                                     mmE,
                                     Vektor<int,DIM>(ddim[0].min(), ddim[1].min()),
                                     ExFieldToDouble());
    Timings::stopTimer(_exchangeTimer);
    Timings::startTimer(_fieldUpdateTimer);

    bc.updateEx(EFD, HFD, JFD, *this, geo);
    EFD.fillGuardCells();

    assign(HFD[II][JJ](0), where(insideMask[II][JJ],
                                 HFD[II][JJ](0) - 0.5 * etax / Z_0 * (EFD[II+1][JJ](1) - EFD[II][JJ](1)),
                                 HFD[II][JJ](0)));
    assign(HFD[II][JJ](1), where(insideMask[II][JJ],
                                 HFD[II][JJ](1) + 0.5 * etay / Z_0 * (EFD[II][JJ+1](0) - EFD[II][JJ](0)),
                                 HFD[II][JJ](1)));
    assign(mmH[mmIIh][mmJJh](0), mmH[mmIIh][mmJJh](0) - 0.5 * etax / Z_0 * (mmE[mmIIh+1][mmJJh](1) - mmE[mmIIh][mmJJh](1)));
    assign(mmH[mmIIh][mmJJh](1), mmH[mmIIh][mmJJh](1) + 0.5 * etay / Z_0 * (mmE[mmIIh][mmJJh+1](0) - mmE[mmIIh][mmJJh](0)));

    Timings::stopTimer(_fieldUpdateTimer);
    Timings::startTimer(_exchangeTimer);

    Communicator::exchangeBoundaries(HFD,
                                     mmH,
                                     Vektor<int,DIM>(ddim[0].min(), ddim[1].min()),
                                     HFieldToDouble());
    Timings::stopTimer(_exchangeTimer);
    Timings::startTimer(_fieldUpdateTimer);

    bc.updateHz(HFD, EFD, *this, geo);
    HFD.fillGuardCells();

    Timings::stopTimer(_fieldUpdateTimer);

    _pbunch.drift_particles(-0.5 * _dt);
}

template<class G, class BC>
void TETMSolver::updateField_secondStep(VField_Edge_t & EFD,
                                        VField_Cell_t & HFD,
                                        VField_Edge_t & JFD,
                                        VField_Edge_t & mmE,
                                        VField_Cell_t & mmH,
                                        VField_Edge_t & mmJ,
                                        NDIndex<DIM> & ddim,
                                        const BC& bc,
                                        G& geo)
{
    Vector_t hr(EFD.get_mesh().get_meshSpacing(0), EFD.get_mesh().get_meshSpacing(1));
    double etax = Physics::c * _dt / hr(0);
    double Z_0 = sqrt(Physics::mu_0 / Physics::epsilon_0);
    double dteps = Physics::c * _dt * Z_0;
    Field<bool, DIM> & insideMask = geo.getInsideMask();

    const NDIndex<DIM> & lDom = EFD.getLayout().getLocalNDIndex();
    const NDIndex<DIM> & mmlDom = mmE.getLayout().getLocalNDIndex();
    Index II = lDom[0];
    Index JJ = lDom[1];
    Index mmII = mmlDom[0];
    Index mmJJ = mmlDom[1];

    _pbunch.push_particles(mmE, mmH, mmJ, _dt, 1);

    Timings::startTimer(_fieldUpdateTimer);

    assign(EFD[II][JJ](1), where(insideMask[II][JJ],
                                 EFD[II][JJ](1) - (etax * Z_0 * (HFD[II][JJ](0) - HFD[II-1][JJ](0) +
                                                                 HFD[II][JJ](1) - HFD[II-1][JJ](1)) +
                                                   dteps * JFD[II][JJ](1)),
                                 EFD[II][JJ](1)));

    assign(mmE[mmII][mmJJ](1), mmE[mmII][mmJJ](1) - (etax * Z_0 * (mmH[mmII][mmJJ](0) - mmH[mmII-1][mmJJ](0) +
                                                                   mmH[mmII][mmJJ](1) - mmH[mmII-1][mmJJ](1)) +
                                                     dteps * mmJ[mmII][mmJJ](1)));

    Timings::stopTimer(_fieldUpdateTimer);
    Timings::startTimer(_exchangeTimer);
    Communicator::exchangeBoundaries(EFD,
                                     mmE,
                                     Vektor<int,DIM>(ddim[0].first(), ddim[1].first()),
                                     EyFieldToDouble());
    Timings::stopTimer(_exchangeTimer);
    Timings::startTimer(_fieldUpdateTimer);

    bc.updateEy(EFD, HFD, JFD, *this, geo);
    EFD.fillGuardCells();

    Timings::stopTimer(_fieldUpdateTimer);

    _pbunch.drift_particles(-0.5 * _dt);
}

#endif
