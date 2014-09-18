/***************************************************************************
                           BoundaryCell.hh
                         -------------------
    begin                : Tue Jun 5 2012
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

#ifndef BOUNDARYCELL_HH
#define BOUNDARYCELL_HH

#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

/** \brief stores information of a single cell that is
 *         cut by the boundary into a part inside and
 *         outside the domain.
 */
struct BoundaryCell {
    size_t localID_m;  /*!< the local id of the cell */
    Vector_t lambda_m; /*!< fraction inside the domain of the edges that belong to that cell */
    double area_m;     /*!< fraction inside the domain of the area */
    Vector_t lambdaToNeighbours_m; /*!< fraction inside the domain of the edges that belong to the neighbouring cells */

    BoundaryCell();

    BoundaryCell(const BoundaryCell& copy);

    BoundaryCell(const size_t& localID,
                 const Vector_t& lambda,
                 const double& area,
                 const Vector_t& lambdaToNeighbours = Vector_t(0.0));

    void print(std::ostream& out, const std::pair<int,int> & coord) const;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int version);
};

inline
BoundaryCell::BoundaryCell():
    localID_m(-1),
    lambda_m(0.0),
    area_m(0.0),
    lambdaToNeighbours_m(0.0)
{ }

inline
BoundaryCell::BoundaryCell(const BoundaryCell& copy):
    localID_m(copy.localID_m),
    lambda_m(copy.lambda_m),
    area_m(copy.area_m),
    lambdaToNeighbours_m(copy.lambdaToNeighbours_m)
{ }

inline
BoundaryCell::BoundaryCell(const size_t& localID,
                           const Vector_t& lambda,
                           const double& area,
                           const Vector_t& lambdaToNeighbours):
    localID_m(localID),
    lambda_m(lambda),
    area_m(area),
    lambdaToNeighbours_m(lambdaToNeighbours)
{ }

inline
void BoundaryCell::print(std::ostream& out, const std::pair<int,int> & coord) const {
    std::ios_base::fmtflags oldFormat = out.flags();
    out << std::setprecision(8) << std::fixed;
    out << "elem:      " << std::setw(8) << coord.first << ", " << coord.second << "\n"
        << "area:      " << std::setw(14) << area_m << "\n"
        << "lambda:    " << std::setw(14) << "(" << lambda_m(0) << " , " << lambda_m(1) << " , "
        << lambdaToNeighbours_m(0) << " , " << lambdaToNeighbours_m(1) << ")\n";

    out << std::endl;
    out.flags(oldFormat);
}

template<class Archive>
void BoundaryCell::serialize(Archive& ar, const unsigned int version) {
    ar & localID_m;
    ar & lambda_m[0];
    ar & lambda_m[1];
    ar & area_m;
    ar & lambdaToNeighbours_m[0];
    ar & lambdaToNeighbours_m[1];
}

#endif
