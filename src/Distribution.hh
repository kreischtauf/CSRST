/***************************************************************************
                           Distribution.hh
                         -------------------
    begin                : Thu Feb 16 2012
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

#ifndef DISTRIBUTION_HH
#define DISTRIBUTION_HH

#include "defs.hh"
#include <vector>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_rng.h>

class LoadH5Part;

class Distribution {
public:
    Distribution(const std::vector<Vector_t> & r,
                 const std::vector<Vector_t> & oldr,
                 const std::vector<Vector_t> & p,
                 const std::vector<double> & q);

    Distribution(TCLAP::CmdLine & cmdl,
                 const VField_Edge_t & EFD,
                 LoadH5Part * fileHandler);

    void save(std::ostream & out);

    void purge();

    size_t size() const;

    const Vector_t & getR(const int & i) const;
    const Vector_t & getOldR(const int & i) const;
    const Vector_t & getP(const int & i) const;
    const double & getQ(const int & i) const;

    const double & getTotalCharge() const;
    const Vector_t & get_pmean() const;

    NDIndex<DIM> getDomain(const VField_Edge_t & EFD) const;

private:
    void loadFromFile(LoadH5Part & fileHandler);

    Vector_t gaussian2DSample(gsl_qrng * q, const double & sigma);
    Vector_t gaussian2DSample(gsl_rng * q, const double & sigma);

    std::vector<Vector_t> R;
    std::vector<Vector_t> oldR;
    std::vector<Vector_t> P;
    std::vector<double> Q;

    double _total_charge;
    Vector_t _momentum;
};

inline
void Distribution::purge()
{
    R.clear();
    oldR.clear();
    P.clear();
    Q.clear();
}

inline
size_t Distribution::size() const
{
    return R.size();
}

inline
const Vector_t & Distribution::getR(const int & i) const
{
    return R.at(i);
}

inline
const Vector_t & Distribution::getOldR(const int & i) const
{
    return oldR.at(i);
}

inline
const Vector_t & Distribution::getP(const int & i) const
{
    return P.at(i);
}

inline
const double & Distribution::getQ(const int & i) const
{
    return Q.at(i);
}

inline
const double & Distribution::getTotalCharge() const
{
    return _total_charge;
}

inline
const Vector_t & Distribution::get_pmean() const
{
    return _momentum;
}
#endif
