/***************************************************************************
                            CalcEnergy.cpp
                         -------------------
    begin                : Tue Jun 16 2009
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

#include "CalcEnergy.hh"
#include "../Physics.hh"

extern std::ofstream dbg;
#define DBGOUT dbg << "CalcEnergy.cpp: " << __LINE__ << "\t"

CalcEnergyCommand::CalcEnergyCommand(const VField_Edge_t & EFD):
    _energy_EFD_old(0.0)
    // _file_name("Data/energy.txt")
{
    // if (Ippl::myNode() == 0) {
    //     _energy_out.open(_file_name.c_str());
    //     _energy_out.precision(10);
    // }
    _executeTimer = Timings::getTimer("calc energy");

    const NDIndex<DIM> lDom = EFD.getLayout().getLocalNDIndex();
    NDIndex<DIM> elem;
    double EE, energyEFD = 0.0;

    _area = EFD.get_mesh().get_meshSpacing(0) * EFD.get_mesh().get_meshSpacing(1);

    for (int j = lDom[1].first(); j <= lDom[1].last(); ++ j) {
        elem[1] = Index(j,j);
        for (int i = lDom[0].first(); i <= lDom[0].last(); ++ i) {
            elem[0] = Index(i,i);
            EE = dot(EFD.localElement(elem),EFD.localElement(elem));
            energyEFD += EE;
        }
    }

    DBGOUT << "using collective communication" << std::endl;
    MPI_Reduce(&energyEFD, &_energy_EFD_old, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    _energy_EFD_old *= _area;
}

CalcEnergyCommand::~CalcEnergyCommand()
{
    // _energy_out.close();
}

double CalcEnergyCommand::execute(const VField_Edge_t & EFD,
                                  const VField_Cell_t & HFD,
                                  const VField_Edge_t & mmEFD,
                                  const VField_Cell_t & mmHFD,
                                  const Vektor<int, DIM> & start,
                                  const double & t,
                                  Inform & msg)
{
    Timings::startTimer(_executeTimer);
    FieldLayout<DIM> & FL = EFD.getLayout();
    const NDIndex<DIM> & lDom = FL.getLocalNDIndex();
    NDIndex<DIM> elem;
    double energy_HFD = 0.0;
    double energy_EFD_new = 0.0;
    double totalEnergy = 0.0;

    double EE, EH;

    for (int j = lDom[1].first(); j <= lDom[1].last(); ++ j) {
        elem[1] = Index(j,j);
        for (int i = lDom[0].first(); i <= lDom[0].last(); ++ i) {
            elem[0] = Index(i,i);
            EE = dot(EFD.localElement(elem),EFD.localElement(elem));
            EH = HFD.localElement(elem)(0) + HFD.localElement(elem)(1);
            energy_EFD_new += EE;
            energy_HFD += EH * EH;
        }
    }

    const NDIndex<DIM> & mmlDom = mmEFD.getLayout().getLocalNDIndex();
    for (int j = mmlDom[1].first(); j <= mmlDom[1].last(); ++ j) {
        elem[1] = Index(j,j);
        for (int i = mmlDom[0].first(); i <= mmlDom[0].last(); ++ i) {
            elem[0] = Index(i,i);
            EE = dot(mmEFD.localElement(elem), mmEFD.localElement(elem));
            EH = mmHFD.localElement(elem)(0) + mmHFD.localElement(elem)(1);
            energy_EFD_new += EE;
            energy_HFD += EH * EH;
        }
    }

    NDIndex<DIM> mmgDom = mmEFD.getLayout().getDomain();
    bool inside = true;
    NDIndex<DIM> dom;
    for (unsigned int d = 0; d < DIM; ++ d) {
        mmgDom[d] = Index(mmgDom[d].first() + start[d], mmgDom[d].last() + start[d]);
        int lower = std::max(lDom[d].first(), mmgDom[d].first());
        int upper = std::min(lDom[d].last(),  mmgDom[d].last());
        if (lower > upper) {
            inside = false;
            break;
        }
        dom[d] = Index(lower, upper);
    }
    if (inside) {
        for (int j = dom[1].first(); j <= dom[1].last(); ++ j) {
            elem[1] = Index(j,j);
            for (int i = dom[0].first(); i <= dom[0].last(); ++ i) {
                elem[0] = Index(i,i);
                EE = dot(EFD.localElement(elem), EFD.localElement(elem));
                EH = HFD.localElement(elem)(0) + HFD.localElement(elem)(1);
                energy_EFD_new -= EE;
                energy_HFD -= EH * EH;
            }
        }
    }

    double energies[] = {energy_EFD_new, energy_HFD};
    double recv[2];
    DBGOUT << "using collective communication" << std::endl;
    MPI_Reduce(energies, recv, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    energy_EFD_new = recv[0] * _area;
    energy_HFD = recv[1] * _area;

    // if (Ippl::myNode() == 0) {
    //     _energy_out << t << "\t"
    //                 << 0.25 * Physics::epsilon_0 * (energy_EFD_new + _energy_EFD_old) + 0.5 * Physics::mu_0 * energy_HFD << "\t"
    //                 << 0.5 * Physics::epsilon_0 * energy_EFD_new << "\t"
    //                 << 0.5 * Physics::mu_0 * energy_HFD << std::endl;
    // }
    totalEnergy = (0.25 * Physics::epsilon_0 * (energy_EFD_new + _energy_EFD_old) + 0.5 * Physics::mu_0 * energy_HFD) * 1e9;

    msg << "\033[0;34m"
        << "Total field energy: "
        << totalEnergy << " nJ"
        << "\033[0m"
        << endl;

    _energy_EFD_old = energy_EFD_new;
    Timings::stopTimer(_executeTimer);

    return totalEnergy;
}
