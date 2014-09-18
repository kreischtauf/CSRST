/***************************************************************************
                               Bunch.hh
                         -------------------
    begin                : Wed Jan 11 2012
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

#ifndef OPAL_BUNCH_HH
#define OPAL_BUNCH_HH

#include "defs.hh"

class Bunch: public ParticleBase<playout_t> {
public:
    ParticleAttrib< Vector_t > oldR;      // old position
    ParticleAttrib< Vector_t > P;      // particle momentum
    ParticleAttrib< double > Q;        // particle charge

    Bunch(playout_t & pl):
        ParticleBase<playout_t>(&pl)
    {
        addAttribute(oldR);
        addAttribute(P);
        addAttribute(Q);
    }

    Bunch()
    {
        addAttribute(oldR);
        addAttribute(P);
        addAttribute(Q);
    }

    /*
    virtual
    void makeReady(const double & dt)
    { }

    virtual
    void push_particles(VField_Edge_t& EFD, VField_Cell_t& BFD, VField_Edge_t& JFD, const double & dt) = 0;

    virtual
    void drift_particles(const double & dt) = 0;

    virtual
    void move_by(const Vector_t & dX) = 0;

    virtual
    void rotate(const Vector_t & X,
                const double& phi) = 0;

    virtual
    void create_particle(const Vector_t & X,
                         const Vector_t & PX,
                         const double & q)
    { }

    virtual
    void get_rmean(Vector_t & spos) const = 0;

    virtual
    void get_pmean(Vector_t & pmean) const = 0;

    virtual
    double get_qtotal() const = 0;

    virtual
    void scatterQ(SField_t & rho) const = 0;

    virtual
    size_t getLocalNP() const {
        return 0;
    }
    */
};
#endif
