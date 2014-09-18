/***************************************************************************
                          GenericGeometry.hh
                         -------------------
    begin                : Wed Jun 6 2012
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

#ifndef GENERICGEOMETRY_HH
#define GENERICGEOMETRY_HH

#include "defs.hh"
#include "BoundaryCell.hh"
#include <vector>

#define GEOMETRY_GCS 2

class GenericGeometry {
public:
    void init(FieldLayout<DIM>& FL, Mesh_t& mesh);

    void makeInsideCell(const size_t& localID);
    void makeOutsideCell(const size_t& localID);

    const std::vector<int>& getInverseMap() const;
    std::vector<int>& getInverseMap();
    const std::vector<int>& getInverseMapDualGrid() const;
    std::vector<int>& getInverseMapDualGrid();
    Field<bool, DIM, Mesh_t, Cell>& getInsideMask();
    const Field<bool, DIM, Mesh_t, Cell>& getInsideMask() const;
    Field<bool, DIM, Mesh_t, Vert>& getInsideMaskDualGrid();
    const Field<bool, DIM, Mesh_t, Vert>& getInsideMaskDualGrid() const;

    Field<char, DIM, Mesh_t, Vert>& getInternalVertices();
    Field<char, DIM, Mesh_t, Cell>& getInternalVerticesDualGrid();
    const Field<char, DIM, Mesh_t, Vert>& getInternalVertices() const;
    const Field<char, DIM, Mesh_t, Cell>& getInternalVerticesDualGrid() const;

    size_t getLocalID(const NDIndex<DIM>& elem) const;
    NDIndex<DIM> getNDIndex(const size_t& localID) const;
    Coordinates getCoordinates(const size_t& localID) const;

    size_t getSizeLocalDomain() const;

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

    std::pair<Vector_t, Vector_t> getLambdas(const int& localID,
                                             const std::vector<BoundaryCell>& boundaryCells) const;

    double getArea(const int& localID,
                   const std::vector<BoundaryCell>& boundaryCells) const;

    NDIndex<DIM> getGlobalDomain() const;
    NDIndex<DIM> getLocalDomain() const;

    std::vector<BoundaryCell>& getGCBoundaryCells();

private:
    size_t _sizeLocalDomain;

    NDIndex<DIM> _lDom;
    NDIndex<DIM> _gDom;

    Field<bool, DIM, Mesh_t, Cell> _insideMask;
    Field<bool, DIM, Mesh_t, Vert> _insideMaskDualGrid;
    Field<char, DIM, Mesh_t, Vert> _internalVertices;
    Field<char, DIM, Mesh_t, Cell> _internalVerticesDualGrid;
    std::vector<int> _localToBCellNr;
    std::vector<int> _localToBCellNrDualGrid;

    std::vector<BoundaryCell> _ghostCellBoundaryCells;
};

inline
void GenericGeometry::init(FieldLayout<DIM>& FL,
                           Mesh_t& mesh) {
    _lDom = FL.getLocalNDIndex();
    _gDom = FL.getDomain();
    _sizeLocalDomain = (_lDom[0].length() + 4) * (_lDom[1].length() + 4);
    _localToBCellNr.resize(_sizeLocalDomain, -1);
    _localToBCellNrDualGrid.resize(_sizeLocalDomain, -1);
    _internalVertices.initialize(mesh, FL, GuardCellSizes<DIM>(GEOMETRY_GCS));
    _internalVerticesDualGrid.initialize(mesh, FL, GuardCellSizes<DIM>(GEOMETRY_GCS));
    _insideMask.initialize(mesh, FL, GuardCellSizes<DIM>(GEOMETRY_GCS));
    _insideMaskDualGrid.initialize(mesh, FL, GuardCellSizes<DIM>(GEOMETRY_GCS));
}

inline
void GenericGeometry::makeInsideCell(const size_t& localID) {
    NDIndex<DIM> elem = getNDIndex(localID);
    _insideMask.localElement(elem) = 1;
}

inline
void GenericGeometry::makeOutsideCell(const size_t& localID) {
    NDIndex<DIM> elem = getNDIndex(localID);
    _insideMask.localElement(elem) = 0;
}

inline
const std::vector<int>& GenericGeometry::getInverseMap() const {
    return _localToBCellNr;
}

inline
std::vector<int>& GenericGeometry::getInverseMap() {
    return _localToBCellNr;
}

inline
const std::vector<int>& GenericGeometry::getInverseMapDualGrid() const {
    return _localToBCellNrDualGrid;
}

inline
std::vector<int>& GenericGeometry::getInverseMapDualGrid() {
    return _localToBCellNrDualGrid;
}

inline
Field<bool, DIM, Mesh_t, Cell>& GenericGeometry::getInsideMask() {
    return _insideMask;
}

inline
const Field<bool, DIM, Mesh_t, Cell>& GenericGeometry::getInsideMask() const {
    return _insideMask;
}

inline
Field<bool, DIM, Mesh_t, Vert>& GenericGeometry::getInsideMaskDualGrid() {
    return _insideMaskDualGrid;
}

inline
const Field<bool, DIM, Mesh_t, Vert>& GenericGeometry::getInsideMaskDualGrid() const {
    return _insideMaskDualGrid;
}

inline
Field<char, DIM, Mesh_t, Vert>& GenericGeometry::getInternalVertices() {
    return _internalVertices;
}

inline
const Field<char, DIM, Mesh_t, Vert>& GenericGeometry::getInternalVertices() const {
    return _internalVertices;
}

inline
Field<char, DIM, Mesh_t, Cell>& GenericGeometry::getInternalVerticesDualGrid() {
    return _internalVerticesDualGrid;
}

inline
const Field<char, DIM, Mesh_t, Cell>& GenericGeometry::getInternalVerticesDualGrid() const {
    return _internalVerticesDualGrid;
}

inline
size_t GenericGeometry::getSizeLocalDomain() const {
    return _sizeLocalDomain;
}

inline
size_t GenericGeometry::getNeighbourNorth(const size_t& localID) const {
    return localID + _lDom[0].length() + 4;
}

inline
size_t GenericGeometry::getNeighbourNorthEast(const size_t& localID) const {
    return localID + _lDom[0].length() + 5;
}

inline
size_t GenericGeometry::getNeighbourEast(const size_t& localID) const {
    return localID + 1;
}

inline
size_t GenericGeometry::getNeighbourSouthEast(const size_t& localID) const {
    return localID - _lDom[0].length() - 3;
}

inline
size_t GenericGeometry::getNeighbourSouth(const size_t& localID) const {
    return localID - _lDom[0].length() - 4;
}

inline
size_t GenericGeometry::getNeighbourSouthWest(const size_t& localID) const {
    return localID - _lDom[0].length() - 5;
}

inline
size_t GenericGeometry::getNeighbourWest(const size_t& localID) const {
    return localID - 1;
}

inline
size_t GenericGeometry::getNeighbourNorthWest(const size_t& localID) const {
    return localID + _lDom[0].length() + 3;
}

inline
size_t GenericGeometry::getNeighbour(const CellNeighbour& neighbour,
                                 const size_t& localID) const {
    switch (neighbour) {
    case WEST:
        return getNeighbourWest(localID);
    case NORTHWEST:
        return getNeighbourNorthWest(localID);
    case NORTH:
        return getNeighbourNorth(localID);
    case NORTHEAST:
        return getNeighbourNorthEast(localID);
    case EAST:
        return getNeighbourEast(localID);
    case SOUTHEAST:
        return getNeighbourSouthEast(localID);
    case SOUTH:
        return getNeighbourSouth(localID);
    case SOUTHWEST:
        return getNeighbourSouthWest(localID);
    default:
        return localID;
    }
}

inline
size_t GenericGeometry::getLocalID(const NDIndex<DIM>& elem) const {
    int idx = elem[0].first() + 2 - _lDom[0].first();
    int idy = elem[1].first() + 2 - _lDom[1].first();

    return idx + idy * (_lDom[0].length() + 4);
}

inline
NDIndex<DIM> GenericGeometry::getNDIndex(const size_t& localID) const {
    int idx = localID % (_lDom[0].length() + 4) - 2 + _lDom[0].first();
    int idy = localID / (_lDom[0].length() + 4) - 2 + _lDom[1].first();

    NDIndex<DIM> elem;
    elem[0] = Index(idx, idx);
    elem[1] = Index(idy, idy);

    return elem;
}

inline
Coordinates GenericGeometry::getCoordinates(const size_t& localID) const {
    Coordinates coords;

    coords.idx_m = localID % (_lDom[0].length() + 4) - 2 + _lDom[0].first();
    coords.idy_m = localID / (_lDom[0].length() + 4) - 2 + _lDom[1].first();

    return coords;
}

inline
NDIndex<DIM> GenericGeometry::getGlobalDomain() const {
    return _gDom;
}

inline
NDIndex<DIM> GenericGeometry::getLocalDomain() const {
    return _lDom;
}

inline
std::pair<Vector_t, Vector_t>
GenericGeometry::getLambdas(const int& localID,
                            const std::vector<BoundaryCell>& boundaryCells) const {
    NDIndex<DIM> elem = getNDIndex(localID);
    const int& idx = elem[0].first();
    const int& idy = elem[1].first();
    const GuardCellSizes<DIM> gcs = _insideMask.getGuardCellSizes();
    if (idx < _lDom[0].first() - (int)gcs.left(0) ||
        idx > _lDom[0].last() + (int)gcs.right(0) ||
        idy < _lDom[1].first() - (int)gcs.left(1) ||
              idy > _lDom[1].last() + (int)gcs.right(1))
        return std::make_pair(Vector_t(0.0), Vector_t(0.0));

    if (_insideMask.localElement(elem)) return std::make_pair(Vector_t(1.0), Vector_t(1.0));

    int bCellNr = _localToBCellNr[localID];
    if (bCellNr == -1) return std::make_pair(Vector_t(0.0), Vector_t(0.0));

    return std::make_pair(boundaryCells[bCellNr].lambda_m,
                          boundaryCells[bCellNr].lambdaToNeighbours_m);
}

inline
double GenericGeometry::getArea(const int& localID,
                                const std::vector<BoundaryCell>& boundaryCells) const {

    NDIndex<DIM> elem = getNDIndex(localID);
    const int& idx = elem[0].first();
    const int& idy = elem[1].first();
    const GuardCellSizes<DIM> gcs = _insideMask.getGuardCellSizes();
    if (idx < _lDom[0].first() - (int)gcs.left(0) ||
        idx > _lDom[0].last() + (int)gcs.right(0) ||
        idy < _lDom[1].first() - (int)gcs.left(1) ||
              idy > _lDom[1].last() + (int)gcs.right(1)) return 0.0;

    if (_insideMask.localElement(elem)) return 1.0;

    int bCellNr = _localToBCellNr[localID];
    if (bCellNr == -1) {
        for (const auto & cell: _ghostCellBoundaryCells)
            if (cell.localID_m == size_t(localID)) return cell.area_m;
        return 0.0;
    }

    return boundaryCells[bCellNr].area_m;
}

inline
std::vector<BoundaryCell>& GenericGeometry::getGCBoundaryCells() {
    return _ghostCellBoundaryCells;
}
#endif
