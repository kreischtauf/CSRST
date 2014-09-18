/***************************************************************************
                               Bend.hh
                         -------------------
    begin                : Wed Apr 27 2011
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

#ifndef BEND_HH
#define BEND_HH

#include "defs.hh"
#include "GenericGeometry.hh"

#include <array>
#include <vector>
#include <map>
#include <algorithm>

struct BoundaryCell;

/** \brief defines the geometry of the domain.
 */
class Bend {
public:
    Bend(const double & lengthBegin,
         const double & length,
         const double & width,
         const double & phi,
         const double & Ekin,
         const double & Bz);

    void Init(VField_t & alpha_e,
              VField_t & alpha_h);

    // void drawBoundary(const std::string & filename, const std::vector<boundaryCell> & bcs);
    void drawDomain(const std::string & filename,
                    const Field<char, DIM, Mesh_t, Vert> & values) const;

    void makeInsideCell(const size_t& localID);
    void makeOutsideCell(const size_t& localID);

    const std::vector<int>& getInverseMap() const;
    std::vector<int>& getInverseMap();
    const std::vector<int>& getInverseMapDualGrid() const;
    std::vector<int>& getInverseMapDualGrid();

    bool isInside(const Vector_t & R) const;
    bool isInsideBend(const Vector_t & R) const;

    const Field<bool, DIM, Mesh_t, Cell>& getInsideMask() const;
    Field<bool, DIM, Mesh_t, Cell>& getInsideMask();
    const Field<bool, DIM, Mesh_t, Vert>& getInsideMaskDualGrid() const;
    Field<bool, DIM, Mesh_t, Vert>& getInsideMaskDualGrid();

    size_t getSizeLocalDomain() const;

    std::pair<Vector_t, Vector_t> getLambdas(const int& localID,
                                             const std::vector<BoundaryCell>& boundaryCells) const;

    double getArea(const int& localID,
                   const std::vector<BoundaryCell>& boundaryCells) const;

    size_t getLocalID(const NDIndex<DIM>& elem) const;
    NDIndex<DIM> getNDIndex(const size_t& localID) const;
    Coordinates getCoordinates(const size_t& localID) const;

    size_t getNeighbourNorth(const size_t& localID) const;
    size_t getNeighbourNorthEast(const size_t& localID) const;
    size_t getNeighbourEast(const size_t& localID) const;
    size_t getNeighbourSouthEast(const size_t& localID) const;
    size_t getNeighbourSouth(const size_t& localID) const;
    size_t getNeighbourSouthWest(const size_t& localID) const;
    size_t getNeighbourWest(const size_t& localID) const;
    size_t getNeighbourNorthWest(const size_t& localID) const;
    size_t getNeighbour(const CellNeighbour& neighbour,
                        const size_t& localID) const;

    void fillCellProperties(std::vector<BoundaryCell>& boundaryCells,
                            std::vector<BoundaryCell>& boundaryCellsDualGrid);

    NDIndex<DIM> getGlobalDomain() const;
    NDIndex<DIM> getLocalDomain() const;

    size_t getLengthStraightSection() const;
    size_t getWidthStraightSection() const;

    void getPositionInCells(Vector_t & position) const;

    std::vector<double> getExternalField(const ParticleAttrib< Vector_t > & R) const;
private:
    double crossBoundary(const Vector_t & base, const Vector_t & dir, const int print = 0);

    void calcInternalVertices();

    void calcInternalVerticesDualGrid();

    void calculatePMLCoefficients(const double& domain_length,
                                  VField_t& alpha_e,
                                  VField_t& alpha_h);

    void getLambdaField(VField_Edge_t & lambdaField, std::vector<int> & bCellNrToLocal);
    void getLambdaFieldDualGrid(VField_Edge_t & lambdaField, std::vector<int> & bCellNrToLocal);
    void getBoundaryCellAreas(std::vector<BoundaryCell> & boundaryCells,
                              std::vector<size_t> & toBeDeleted,
                              const VField_Edge_t & lambdaField,
                              const std::vector<int> & bCellNrToLocal);
    void getBoundaryCellAreasDualGrid(std::vector<BoundaryCell> & boundaryCells,
                                      std::vector<size_t> & toBeDeleted,
                                      const VField_Edge_t & lambdaField,
                                      const std::vector<int> & bCellNrToLocal);
    void getBoundaryCellAreasGhostCells(const VField_Edge_t & lambdaField);
    void getBoundaryCellAreasGhostCellsImpl(const VField_Edge_t & lambda, const NDIndex<DIM> & dom);

    static
    double areaCases(const int & ccase, const Vector_t & lambda, const Vector_t & lambdaNeighbours);

    const double _mass; // rest mass in MeV
    double _domain_length;
    double _domain_width;
    double _dlengthAfter; // length of drift after bend
    double _dlengthBefore; // length of drift before bend
    double _vwidth; // width of vacuum chamber
    double _phi;  // bending angle
    double _ekin; // kinetic energy in MeV
    double _bz;  // magnetic field perpendicular to plane
    double _R; // radius of bend
    Vector_t _hr;
    Vector_t _center;
    Vector_t _mesh_origin;

    double _innerCutoffDistance;
    double _outerCutoffDistance;

    GenericGeometry _geo;
};



inline
void Bend::makeInsideCell(const size_t& localID) {
    _geo.makeInsideCell(localID);
}

inline
void Bend::makeOutsideCell(const size_t& localID) {
    _geo.makeOutsideCell(localID);
}

inline
const std::vector<int>& Bend::getInverseMap() const {
    return _geo.getInverseMap();
}

inline
std::vector<int>& Bend::getInverseMap() {
    return _geo.getInverseMap();
}

inline
const std::vector<int>& Bend::getInverseMapDualGrid() const {
    return _geo.getInverseMapDualGrid();
}

inline
std::vector<int>& Bend::getInverseMapDualGrid() {
    return _geo.getInverseMapDualGrid();
}

inline
const Field<bool, DIM, Mesh_t, Cell>& Bend::getInsideMask() const {
    return _geo.getInsideMask();
}

inline
Field<bool, DIM, Mesh_t, Cell>& Bend::getInsideMask() {
    return _geo.getInsideMask();
}

inline
const Field<bool, DIM, Mesh_t, Vert>& Bend::getInsideMaskDualGrid() const {
    return _geo.getInsideMaskDualGrid();
}

inline
Field<bool, DIM, Mesh_t, Vert>& Bend::getInsideMaskDualGrid() {
    return _geo.getInsideMaskDualGrid();
}

inline
size_t Bend::getSizeLocalDomain() const {
    return _geo.getSizeLocalDomain();
}

inline
size_t Bend::getNeighbourNorth(const size_t& localID) const {
    return _geo.getNeighbourNorth(localID);
}

inline
size_t Bend::getNeighbourNorthEast(const size_t& localID) const {
    return _geo.getNeighbourNorthEast(localID);
}

inline
size_t Bend::getNeighbourEast(const size_t& localID) const {
    return _geo.getNeighbourEast(localID);
}

inline
size_t Bend::getNeighbourSouthEast(const size_t& localID) const {
    return _geo.getNeighbourSouthEast(localID);
}

inline
size_t Bend::getNeighbourSouth(const size_t& localID) const {
    return _geo.getNeighbourSouth(localID);
}

inline
size_t Bend::getNeighbourSouthWest(const size_t& localID) const {
    return _geo.getNeighbourSouthWest(localID);
}

inline
size_t Bend::getNeighbourWest(const size_t& localID) const {
    return _geo.getNeighbourWest(localID);
}

inline
size_t Bend::getNeighbourNorthWest(const size_t& localID) const {
    return _geo.getNeighbourNorthWest(localID);
}

inline
size_t Bend::getNeighbour(const CellNeighbour& neighbour,
                          const size_t& localID) const {
    return _geo.getNeighbour(neighbour, localID);
}

inline
size_t Bend::getLocalID(const NDIndex<DIM>& elem) const {
    return _geo.getLocalID(elem);
}

inline
NDIndex<DIM> Bend::getNDIndex(const size_t& localID) const {
    return _geo.getNDIndex(localID);
}

inline
Coordinates Bend::getCoordinates(const size_t& localID) const {
    return _geo.getCoordinates(localID);
}

inline
NDIndex<DIM> Bend::getGlobalDomain() const {
    return _geo.getGlobalDomain();
}

inline
NDIndex<DIM> Bend::getLocalDomain() const {
    return _geo.getLocalDomain();
}

inline
std::pair<Vector_t, Vector_t> Bend::getLambdas(const int& localID,
                                               const std::vector<BoundaryCell>& boundaryCells) const {
    return _geo.getLambdas(localID, boundaryCells);
}

inline
double Bend::getArea(const int& localID,
                     const std::vector<BoundaryCell>& boundaryCells) const {
    return _geo.getArea(localID, boundaryCells);
}

inline
size_t Bend::getLengthStraightSection() const
{
    return _dlengthBefore / _hr(0);
}

inline
size_t Bend::getWidthStraightSection() const
{
    return _vwidth / _hr(1);
}

inline
void Bend::getPositionInCells(Vector_t & position) const
{
    position = (position + _mesh_origin) / _hr;

    return;
}

inline
bool Bend::isInsideBend(const Vector_t & R) const
{
    Tenzor<double, DIM> M(cos(_phi), -sin(_phi), sin(_phi), cos(_phi));

    Vector_t Dx = R - _center;
    Vector_t Dxrot = dot(M, Dx);

    if (Dx(0) < 0.0 || Dxrot(0) > 0.0)
        return false;

    if (Dx(0) > 0.0 && Dxrot(0) < 0.0)
        return true;

    return false;
}

#endif // BEND_HH
