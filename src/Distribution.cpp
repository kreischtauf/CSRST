/***************************************************************************
                           Distribution.cpp
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

#include "Distribution.hh"
#include "Commands/LoadH5Part.hh"
#include "Physics.hh"
#include "utils.hh"

#include <gsl/gsl_randist.h>
#include <boost/mpi.hpp>

#define MASS (Physics::m_e * 1e3)

extern std::ostream dbg;
#define DBGOUT dbg << __FILE__ << ": " << __LINE__

Distribution::Distribution(const std::vector<Vector_t> & r,
                           const std::vector<Vector_t> & oldr,
                           const std::vector<Vector_t> & p,
                           const std::vector<double> & q):
    R(r.begin(), r.end()),
    oldR(oldr.begin(), oldr.end()),
    P(p.begin(), p.end()),
    Q(q.begin(), q.end())
{
    size_t numLocalParticles = std::min(std::min(R.size(), oldR.size()), std::min(P.size(), Q.size()));
    size_t numTotalParticles;
    R.resize(numLocalParticles);
    oldR.resize(numLocalParticles);
    P.resize(numLocalParticles);
    Q.resize(numLocalParticles);

    for (size_t i = 0; i < numLocalParticles; ++ i) {
        _momentum += P[i];
    }
    reduce(&_momentum[0], &_momentum[0] + DIM, &_momentum, OpAddAssign());
    reduce(numLocalParticles, numTotalParticles, OpAddAssign());
    _momentum /= numTotalParticles;
}

Distribution::Distribution(TCLAP::CmdLine & cmdl,
                           const VField_Edge_t & EFD,
                           LoadH5Part * fileHandler):
    _total_charge(0.0),
    _momentum(0.0)
{
    if (fileHandler != 0) {
        loadFromFile(*fileHandler);
    } else {
        std::list<TCLAP::Arg*> callArguments = cmdl.getArgList();
        int numTotalParticles = 0;
        Vect2D mean_x(0.0, 0.0);
        double sigma_x = 0.0;
        double ekin = 0.0;
        double sigma_p = 0.0;
        for (TCLAP::Arg* arg: callArguments) {
            if (arg->getName() == "Np") {
                TCLAP::ValueArg<int> & Varg = *dynamic_cast<TCLAP::ValueArg<int>*>(arg);
                numTotalParticles = Varg.getValue();
            } else if (arg->getName() == "mean_Rx") {
                TCLAP::ValueArg<double> & Varg = *dynamic_cast<TCLAP::ValueArg<double>*>(arg);
                mean_x.v[0] = Varg.getValue();
            } else if (arg->getName() == "mean_Ry") {
                TCLAP::ValueArg<double> & Varg = *dynamic_cast<TCLAP::ValueArg<double>*>(arg);
                mean_x.v[1] = Varg.getValue();
            } else if (arg->getName() == "sigma_x") {
                TCLAP::ValueArg<double> & Varg = *dynamic_cast<TCLAP::ValueArg<double>*>(arg);
                sigma_x = Varg.getValue();
            } else if (arg->getName() == "Ekin") {
                TCLAP::ValueArg<double> & Varg = *dynamic_cast<TCLAP::ValueArg<double>*>(arg);
                ekin = Varg.getValue();
            } else if (arg->getName() == "sigma_p") {
                TCLAP::ValueArg<double> & Varg = *dynamic_cast<TCLAP::ValueArg<double>*>(arg);
                sigma_p = Varg.getValue();
            } else if (arg->getName() == "Qtotal") {
                TCLAP::ValueArg<double> & Varg = *dynamic_cast<TCLAP::ValueArg<double>*>(arg);
                _total_charge = -Varg.getValue() * 1e-9;
            }
        }
        size_t numLocalParticles = std::floor(static_cast<double>(numTotalParticles) / Ippl::getNodes());
        size_t missing = numTotalParticles - numLocalParticles * Ippl::getNodes();
        if (static_cast<size_t>(Ippl::myNode()) < missing) {
            numLocalParticles ++;
        }

        R.resize(numLocalParticles);
        oldR.resize(numLocalParticles, 0.0);
        P.resize(numLocalParticles);
        Q.resize(numLocalParticles, _total_charge / numTotalParticles);

        double gamma = ekin / MASS + 1;
        double mean_p = std::sqrt(gamma*gamma - 1);

        gsl_qrng * p = gsl_qrng_alloc(gsl_qrng_halton, 2);
        gsl_qrng * q = gsl_qrng_alloc(gsl_qrng_halton, 2);
        size_t ii = 0;
        for (int i = 0; i < numTotalParticles; ++ i) {
            Vector_t x = gaussian2DSample(p, sigma_x);
            Vector_t p = gaussian2DSample(q, sigma_p);
            if (i % Ippl::getNodes() != Ippl::myNode()) continue;

            R[ii] = mean_x.toVector_t() + x;
            P[ii] = Vector_t(mean_p, 0.0) + p;
            _momentum += P[ii];
            ++ ii;
        }

        std::sort(R.begin(), R.end(), Utils::XFirst(EFD.get_mesh()));

        reduce(&_momentum[0], &_momentum[0] + DIM, &_momentum[0], OpAddAssign());
        _momentum /= numTotalParticles;

        gsl_qrng_free(p);
        gsl_qrng_free(q);
    }
}

void Distribution::loadFromFile(LoadH5Part & fileHandler)
{
    fileHandler.restoreBunch(R, oldR, P, Q);
    for (auto & p: P) {
        _momentum += p;
    }

    size_t numParticles = P.size();
    reduce(numParticles, numParticles, OpAddAssign());
    reduce(&_momentum[0], &_momentum[0] + DIM, &_momentum[0], OpAddAssign());

    _momentum /= numParticles;
}

Vector_t Distribution::gaussian2DSample(gsl_qrng * q, const double & sigma)
{
    double v[2];
    gsl_qrng_get(q, v);

    double dist = sigma * sqrt(- 2 * log(v[0]));
    double angle = Physics::two_pi * v[1];
    return Vector_t(dist * cos(angle), dist * sin(angle));
}

Vector_t Distribution::gaussian2DSample(gsl_rng * q, const double & sigma)
{
    // double dist = sigma * sqrt(- 2 * log(gsl_rng_uniform(q)));
    double dist = gsl_ran_gaussian_tail(q, 0.0, sigma);
    double angle = Physics::two_pi * gsl_rng_uniform(q);

    return Vector_t(dist * cos(angle), dist * sin(angle));
}

void Distribution::save(std::ostream & out)
{
    std::streamsize old_precision = out.precision();
    out.precision(10);
    for (size_t i = 0; i < size(); ++ i) {
        out << std::setw(20) << R[i](0)
            << std::setw(20) << P[i](0)
            << std::setw(20) << R[i](1)
            << std::setw(20) << P[i](1)
            << std::setw(20) << Q[i]
            << std::setw(20) << 0.0
            << "\n";
    }
    out.precision(old_precision);
}

NDIndex<DIM> Distribution::getDomain(const VField_Edge_t & EFD) const
{
    boost::mpi::communicator world;

    Vector_t minR(R[0]), maxR(R[0]);
    Vector_t origin = EFD.get_mesh().get_origin();
    Vector_t dx(EFD.get_mesh().get_meshSpacing(0),
                EFD.get_mesh().get_meshSpacing(1));
    for (size_t i = 0; i < size(); ++ i) {
        minR(0) = std::min(minR(0), R[i](0));
        minR(1) = std::min(minR(1), R[i](1));
        maxR(0) = std::max(maxR(0), R[i](0));
        maxR(1) = std::max(maxR(1), R[i](1));
    }
    int minima[2 * DIM], result[2 * DIM];
    for (unsigned int d = 0; d < DIM; ++ d) {
        minima[2 * d    ] = static_cast<int>(floor((minR(d) - origin(d)) / dx(d)));
        minima[2 * d + 1] = -static_cast<int>(floor((maxR(d) - origin(d)) / dx(d)));
    }
    boost::mpi::all_reduce(world, minima, 2 * DIM,  result, boost::mpi::minimum<int>());

    NDIndex<DIM> ret;
    for (unsigned int d = 0; d < DIM; ++ d) {
        ret[d] = Index(result[2 * d], -result[2 * d + 1]);
    }
    return ret;
}
