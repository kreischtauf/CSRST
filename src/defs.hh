/***************************************************************************
                                  defs.hh
                             -------------------
    begin                : Tue Oct 6 2009
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

#ifndef DEFS_HH
#define DEFS_HH

#include "csrst_config.h"
#include "Ippl.h"
#include <fstream>
#include <boost/smart_ptr.hpp>
#include "tclap/CmdLine.h"
#include "Timings.hh"

#define DIM 2
#define GUARDCELLSIZE 1
#define SPACE_EPS 1e-10
#define MINAREA 1e-20

#define indent_l0 ""
#define indent_l1 "  "
#define indent_l2 "    "
#define indent_l3 "      "
#define indent_l4 "        "
#define indent_l5 "          "

#define dontGETUNIFORMLAYOUT

typedef unsigned int                                Idx_t;
typedef double                                      Scalar_t;
typedef double                                      SolverPrecision_t;
typedef std::pair<int,int>                          coord_t;

#ifdef GETUNIFORMLAYOUT
typedef ParticleUniformLayout<double, DIM>          playout_t;
#else
typedef ParticleSpatialLayout<double, DIM>          playout_t;
#endif
typedef playout_t::SingleParticlePos_t              Vector_t;
typedef Vektor<double, 3>                           Vector3D_t;

typedef IntCIC                                      IntOp_t;
// typedef IntS1                                       IntOp_t;
// typedef IntS2                                       IntOp_t;
typedef UniformCartesian<DIM,double>                Mesh_t;
typedef Vert                                        Center_t;
typedef CenteredFieldLayout<DIM, Mesh_t, Center_t>  FieldLayout_t;
typedef CenteredFieldLayout<DIM, Mesh_t, Edge>      FieldLayout_Edge_t;
typedef CenteredFieldLayout<DIM, Mesh_t, Cell>      FieldLayout_Cell_t;
typedef CenteredFieldLayout<DIM, Mesh_t, Vert>      FieldLayout_Vert_t;
typedef Field<Vector_t, DIM, Mesh_t, Center_t>      VField_t;
typedef Field<Vector_t, DIM, Mesh_t, Edge>          VField_Edge_t;
typedef Field<Vector_t, DIM, Mesh_t, Cell>          VField_Cell_t;
typedef Field<Vector_t, DIM, Mesh_t, Vert>          VField_Vert_t;
typedef Field<Scalar_t, DIM, Mesh_t, Center_t>      SField_t;
typedef Field<Scalar_t, DIM, Mesh_t, Cell>          SField_Cell_t;
typedef Field<Scalar_t, DIM, Mesh_t, Vert>          SField_Vert_t;
typedef GuardCellSizes<DIM>                         GCS_t;
typedef boost::shared_ptr<SField_Vert_t>            SField_Vert_ptr;
typedef boost::shared_ptr<SField_Cell_t>            SField_Cell_ptr;

enum BCType {PEC, ABC, PBC};

enum CellNeighbour {WEST = 0,
                    NORTHWEST,
                    NORTH,
                    NORTHEAST,
                    EAST,
                    SOUTHEAST,
                    SOUTH,
                    SOUTHWEST,
                    SELF};

// Define a simple 3D vector type
struct Vect3D {
    double v[3];

    Vect3D(const double & vx,
             const double & vy,
             const double & vz) {
        v[0] = vx;
        v[1] = vy;
        v[2] = vz;
    }

    // operator= will be used to assign to the vector
    Vect3D& operator=(const std::string &str)
    {
	std::istringstream iss(str);
	if (!(iss >> v[0] >> v[1] >> v[2]))
	    throw TCLAP::ArgParseException(str + " is not a 3D vector");

	return *this;
    }

    Vector3D_t toVector3D_t() const
    {
        return Vector3D_t(v[0], v[1], v[2]);
    }

    std::ostream& print(std::ostream &os) const
    {
	std::copy(v, v + 3, std::ostream_iterator<double>(os, " "));
	return os;
    }
};

struct Vect2D {
    double v[2];

    Vect2D(const double & vx,
           const double & vy) {
        v[0] = vx;
        v[1] = vy;
    }

    Vect2D(const Vect2D & o) {
        v[0] = o.v[0];
        v[1] = o.v[1];
    }

    // operator= will be used to assign to the vector
    Vect2D& operator=(const std::string &str)
    {
	std::istringstream iss(str);
	if (!(iss >> v[0] >> v[1]))
	    throw TCLAP::ArgParseException(str + " is not a 2D vector");

	return *this;
    }

    double operator()(const unsigned int & i) const
    {
        if (i >= 2) return 0.0;
        return v[i];
    }

    Vector_t toVector_t() const
    {
        return Vector_t(v[0], v[1]);
    }

    std::ostream& print(std::ostream &os) const
    {
	std::copy(v, v + 2, std::ostream_iterator<double>(os, " "));
	return os;
    }
};

inline
std::ostream& operator<<(std::ostream& out, const Vect2D& vec) {
    vec.print(out);
    return out;
}

struct Particle {
    Particle() {
        a[0] = 0.0;
        a[1] = 0.0;
        a[2] = 0.0;
        a[3] = 0.0;
        a[4] = 0.0;
        a[5] = 0.0;
        a[6] = 0.0;
    }
    Particle(const Vector_t & olR, const Vector_t & neR, const Vector_t & P, const double & Q) {
        a[0] = olR(0);
        a[1] = olR(1);
        a[2] = neR(0);
        a[3] = neR(1);
        a[4] = P(0);
        a[5] = P(1);
        a[6] = Q;
    }
    Vector_t getOldR() const {
        return Vector_t(a[0], a[1]);
    }
    Vector_t getR() const {
        return Vector_t(a[2], a[3]);
    }
    Vector_t getP() const {
        return Vector_t(a[4], a[5]);
    }
    double getQ() const {
        return a[6];
    }

    double a[7];
};

struct Coordinates {
    int idx_m;
    int idy_m;

    void print(std::ostream& out) const {
        out << idx_m << ", " << idy_m;
    }
};

inline
std::ostream& operator<<(std::ostream& out, const Coordinates& coord) {
    coord.print(out);
    return out;
}

#endif
