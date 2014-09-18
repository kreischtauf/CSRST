/***************************************************************************
                               Bend.cpp
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

#include "Bend.hh"
#include "Physics.hh"
#include "gauss_legendre.h"
#include "BoundaryCell.hh"
#include "utils.hh"

#include <algorithm>
#include <initializer_list>
#include <utility>

#define PML_WIDTH 10
// reflection of PML is 10^{-PML_RED}
#define PML_RED exp(-16.)
// sigma scales as \rho^PML_ORD
#define PML_ORD 4
extern std::ofstream dbg;
#define DBGOUT dbg << "Bend.cpp: " << __LINE__ << "\t"

double power(double x, int b) {
    double c = 1.0;

    while (b != 0) {
        if (b % 2 == 1) {
            c *= x;
        }
        x *= x;
        b /= 2;
    }
    return c;
}

double sigma(double x, double y, void* data) {
    const double sinphi = *(double*)data;
    const double cosphi = *((double*)data + 1);
    const double pmlwidth = *((double*)data + 2);

    double val = (cosphi * x + sinphi * y) / pmlwidth;
    val = val > 0.0? val: 0.0;
    val = val < 1.0? power(val, PML_ORD): 0.0;
    return val;
}

Bend::Bend(const double & lengthBefore,
           const double & lengthAfter,
           const double & width,
           const double & phi,
           const double & Ekin,
           const double & Bz):
    _mass(0.511e6),
    _dlengthAfter(lengthAfter),
    _dlengthBefore(lengthBefore),
    _vwidth(width),
    _ekin(Ekin*1e6),
    _bz(Bz),
    _center(0.0),
    _innerCutoffDistance(0.0),
    _outerCutoffDistance(0.0)
{
    if (std::abs(_bz) > 1e-8) {
        _R = sqrt(_ekin*_ekin + 2. * _ekin * _mass) / (Physics::c * std::abs(_bz));
        _phi = -Bz / std::abs(Bz) * std::abs(phi) * Physics::pi / 180.;
    } else {
        _R = 9999999.9;
        _phi = 0.0;
    }
}

void Bend::Init(VField_t& alpha_e,
                VField_t& alpha_h)
{
    Mesh_t& mesh = alpha_e.get_mesh();
    FieldLayout<DIM>& FL = alpha_e.getLayout();

    _geo.init(FL, mesh);

    const NDIndex<DIM> gDom = _geo.getGlobalDomain();

    const double cosphi = cos(_phi);
    const double sinphi = std::abs(sin(_phi));
    if (std::abs(_bz) > 1e-8) {
        _domain_length = (_dlengthBefore + _dlengthAfter * cosphi + (_R + 0.5 * _vwidth) * sinphi);
        _domain_width = ((_R - 0.5 * _vwidth) * (1 - cosphi) + _vwidth + _dlengthAfter * sinphi);
    } else {
        _domain_length = _dlengthAfter + _dlengthBefore;
        _domain_width = _vwidth;
    }
    _hr = Vector_t(_domain_length / (gDom[0].length() - 2),
                   _domain_width / (gDom[1].length() - 2));

    if (_bz <= 0.0) {
        _mesh_origin = Vector_t(0, _vwidth - _domain_width) - _hr;
        _center = Vector_t(_dlengthBefore, 0.5 *_vwidth - _R);
    } else {
        _mesh_origin = -_hr;
        _center = Vector_t(_dlengthBefore, 0.5 *_vwidth + _R);
    }

    _domain_length += _hr(0);
    _domain_width += _hr(1);

    mesh.set_meshSpacing(&(_hr(0)));
    mesh.set_origin(_mesh_origin);

    calcInternalVertices();
    calcInternalVerticesDualGrid();
    calculatePMLCoefficients(_domain_length, alpha_e, alpha_h);
}

double Bend::crossBoundary(const Vector_t & base, const Vector_t & dir, const int print) {
    const double cosphi = cos(_phi);
    const double sinphi = sin(_phi);

    const Vector_t x = base + _mesh_origin;
    const Vector_t Dx = x - _center;

    if (Dx(0) < 0.0) { // in first straight section
        double t = 99.9;
        if (x(1) < _vwidth + SPACE_EPS && x(1) + dir(1) > _vwidth - SPACE_EPS) {
            t = std::max(1.0, std::min(0.0, (_vwidth - x(1)) / dir(1)));

        }
        if (x(1) > -SPACE_EPS && x(1) + dir(1) < SPACE_EPS) {
            t = std::min(t, -x(1) / dir(1));
        }
        if (x(0) > -SPACE_EPS && x(0) + dir(0) < SPACE_EPS) {
            t = std::min(t, -x(0) / dir(0));
        }
        return t;
    } else if (Dx(1) * sinphi > Dx(0) * cosphi) { // in bend section
        double rmin_sq = (_R - 0.5 * _vwidth) * (_R - 0.5 * _vwidth);
        double rmax_sq = (_R + 0.5 * _vwidth) * (_R + 0.5 * _vwidth);

        double dist1_sq = dot(Dx, Dx);
        double dist2_sq = dot(Dx + dir, Dx + dir);
        double a = dot(dir, dir);
        double b = dot(Dx, dir);
        double c = 0.0;
        if (dist1_sq < rmax_sq && dist2_sq > rmax_sq) {
            c = dot(Dx, Dx) - rmax_sq;
        } else if (dist1_sq > rmin_sq && dist2_sq < rmin_sq) {
            c = dot(Dx, Dx) - rmin_sq;
        } else {
            return 99.9;
        }
        double t = 99.9;
        double tau = (-b + sqrt(b*b - a * c)) / a;
        if (std::abs(tau - 0.5) < 0.5)
            t = tau;
        tau = -tau - 2 * b / a;
        if (std::abs(tau - 0.5) < 0.5)
            t = std::min(t, tau);
        return t;
    } else { // in second straight section
        Vector_t xrot(0.0);
        Vector_t dirrot(0.0);

        if (std::abs(_phi) > SPACE_EPS) {
            xrot(0) = cosphi * Dx(0) - sinphi * Dx(1) + _center(0);
            xrot(1) = sinphi * Dx(0) + cosphi * Dx(1) + _center(1);

            dirrot(0) = cosphi * dir(0) - sinphi * dir(1);
            dirrot(1) = sinphi * dir(0) + cosphi * dir(1);
        } else {
            xrot = x;
            dirrot = dir;
        }

        double t0 = 99.9, t1 = 99.9, t2 = 99.9;
        if (std::abs(dirrot(1)) > SPACE_EPS) {
            t0 = -xrot(1) / dirrot(1);
            if (t0 < 0.0) t0 = 99.9;
            t1 = (_vwidth - xrot(1)) / dirrot(1);
            if (t1 < 0.0) t1 = 99.9;
        }
        if (std::abs(dirrot(0)) > SPACE_EPS) {
            t2 = (_dlengthBefore + _dlengthAfter - xrot(0)) / dirrot(0);
            if (t2 < 0.0) t2 = 99.9;
        }

        double t = std::min(t0, std::min(t1, t2));
        if (t > 1 + SPACE_EPS && t < 1.001) t = 1.0;

        PAssert(t < 1 + SPACE_EPS);
        PAssert(t >= 0);
        return std::min(1.0, std::max(0.0, t));
    }
}

bool Bend::isInside(const Vector_t & x) const
{
    static const double cosphi = cos(_phi);
    static const double sinphi = sin(_phi);
    static const double rmin_sq = (_R - 0.5 * _vwidth) * (_R - 0.5 * _vwidth);
    static const double rmax_sq = (_R + 0.5 * _vwidth) * (_R + 0.5 * _vwidth);
    static const Tenzor<double, DIM> rotation(cosphi, -sinphi, sinphi, cosphi);

    Vector_t Dx = x - _center;

    if (Dx(0) < 0.0) {
        if (x(0) < SPACE_EPS ||
            x(1) < SPACE_EPS ||
            x(1) > _vwidth - SPACE_EPS) return false;
    } else if (Dx(1) * sinphi > Dx(0) * cosphi) {
        double length_sq = dot(Dx, Dx);
        if (length_sq <= rmin_sq || length_sq >= rmax_sq) return false;
    } else {
        Vector_t dxp;
        if (std::abs(_phi) > SPACE_EPS) {
            dxp = dot(rotation, Dx) + _center;
        } else {
            dxp = x;
        }
        if (dxp(1) < SPACE_EPS ||
            dxp(1) > _vwidth - SPACE_EPS ||
            dxp(0) > _dlengthBefore + _dlengthAfter - SPACE_EPS) return false;
    }
    return true;
}

void Bend::calcInternalVertices() {
    const NDIndex<DIM> lDom = _geo.getLocalDomain();

    auto& internalVertices = _geo.getInternalVertices();
    internalVertices = 0;

    Vector_t x;
    NDIndex<DIM> elem;

    for (int cy = lDom[1].first(); cy <= lDom[1].last(); cy++) {
        elem[1] = Index(cy,cy);
        x(1) = _mesh_origin(1) + cy * _hr(1);
        for (int cx = lDom[0].first(); cx <= lDom[0].last(); cx++) {
            elem[0] = Index(cx,cx);
            x(0) = _mesh_origin(0) + cx * _hr(0);

            internalVertices.localElement(elem) = isInside(x) ? 1: 0;
        }
    }

    internalVertices.fillGuardCells();
}

void Bend::calcInternalVerticesDualGrid() {
    const NDIndex<DIM> lDom = _geo.getLocalDomain();

    auto& internalVertices = _geo.getInternalVerticesDualGrid();
    internalVertices = 0;

    Vector_t x;
    NDIndex<DIM> elem;

    for (int cy = lDom[1].first(); cy <= lDom[1].last(); cy++) {
        elem[1] = Index(cy,cy);
        x(1) = _mesh_origin(1) + (cy + 0.5) * _hr(1);
        for (int cx = lDom[0].first(); cx <= lDom[0].last(); cx++) {
            elem[0] = Index(cx,cx);
            x(0) = _mesh_origin(0) + (cx + 0.5) * _hr(0);

            internalVertices.localElement(elem) = isInside(x) ? 1: 0;
        }
    }
    internalVertices.fillGuardCells();
}

void Bend::fillCellProperties(std::vector<BoundaryCell>& boundaryCells,
                              std::vector<BoundaryCell>& boundaryCellsDualGrid) {
    const auto gDom = _geo.getGlobalDomain();

    std::vector<int> bCellNrToLocal;
    std::vector<int> bCellNrToLocalDualGrid;

    auto& localToBCellNr = _geo.getInverseMap();
    auto& localToBCellNrDualGrid = _geo.getInverseMapDualGrid();

    e_dim_tag decomp[] = {PARALLEL, PARALLEL, SERIAL};
    double spacing[] = {1.0, 1.0};

    Mesh_t mesh(gDom[0],
                gDom[1],
                spacing,
                Vector_t(0.0));

    FieldLayout_Edge_t FL_edge(mesh, decomp);

    VField_Edge_t lambdaField(mesh,
                              FL_edge,
                              GuardCellSizes<DIM>(GEOMETRY_GCS));
    lambdaField = 0.0;
    getLambdaField(lambdaField, bCellNrToLocal);

    std::vector<size_t> toBeDeleted;
    getBoundaryCellAreas(boundaryCells, toBeDeleted, lambdaField, bCellNrToLocal);

    if (toBeDeleted.size() > 0) {
        size_t j = 0;
        for (size_t i = 0; i < bCellNrToLocal.size(); ++ i) {
            if (bCellNrToLocal[i] == static_cast<long>(toBeDeleted[j])) {
                ++ j;
                continue;
            }

            localToBCellNr[bCellNrToLocal[i]] = i - j;
        }
    } else {
        for (size_t i = 0; i < bCellNrToLocal.size(); ++ i) {
            localToBCellNr[bCellNrToLocal[i]] = i;
        }
    }

    toBeDeleted.clear();
    getBoundaryCellAreasGhostCells(lambdaField);

    lambdaField = 0.0;
    getLambdaFieldDualGrid(lambdaField, bCellNrToLocalDualGrid);

    getBoundaryCellAreasDualGrid(boundaryCellsDualGrid,
                                 toBeDeleted,
                                 lambdaField,
                                 bCellNrToLocalDualGrid);

    if (toBeDeleted.size() > 0) {
        size_t j = 0;
        for (size_t i = 0; i < bCellNrToLocalDualGrid.size(); ++ i) {
            if (bCellNrToLocalDualGrid[i] == static_cast<long>(toBeDeleted[j])) {
                ++ j;
                continue;
            }

            localToBCellNrDualGrid[bCellNrToLocalDualGrid[i]] = i - j;
        }
    } else {
        for (size_t i = 0; i < bCellNrToLocalDualGrid.size(); ++ i) {
            localToBCellNrDualGrid[bCellNrToLocalDualGrid[i]] = i;
        }
    }
}

void Bend::calculatePMLCoefficients(const double& domain_length,
                                    VField_t& alpha_e,
                                    VField_t& alpha_h) {
    const auto lDom = _geo.getLocalDomain();
    const auto gDom = _geo.getGlobalDomain();
    const double sinphi{0.5 * _phi};
    const double cosphi{0.5 * _phi};
    const double Z_0 = sqrt(Physics::mu_0 / Physics::epsilon_0);

    alpha_e = Vector_t(0.0);
    alpha_h = Vector_t(0.0);

    NDIndex<DIM> elem;
    {
        double data[3] = {sinphi, cosphi, (cosphi * _hr(0) + sinphi * _hr(1)) * PML_WIDTH};
        double P[2] = {sinphi * _vwidth + cosphi * data[2], sinphi * data[2]}; // add here _hr(0) due to pec
        double sigma_m = -(PML_ORD + 1) * log(PML_RED) / (2.0 * Z_0 * data[2] * _hr(0) * _hr(1));
        int BBx1 = 0;
        int BBx2 = (int)floor(P[0] / _hr(0) + 0.99);
        int BBy1 = 0;
        int BBy2 = (int)floor((cosphi * _vwidth + sinphi * data[2]) / _hr(1) + 0.99);
        int X1 = lDom[0].first() > BBx1? lDom[0].first(): BBx1;
        int X2 = lDom[0].last() < BBx2? lDom[0].last(): BBx2;
        int Y1 = lDom[1].first() > BBy1? lDom[1].first(): BBy1;
        int Y2 = lDom[1].last() < BBy2? lDom[1].last(): BBy2;


        for (int X = X1; X <= X2; X ++) {
            elem[0] = Index(X, X);
            for (int Y = Y1; Y <= Y2; Y ++) {
                elem[1] = Index(Y, Y);
                double dx = P[0] - sinphi * (Y * _hr(1) - P[1]) / cosphi - X * _hr(0);
                double sigex = sigma_m * gauss_legendre_2D_cube((PML_ORD + 1), sigma, (void*)data, dx              , dx       + _hr(0), -0.5 * _hr(1), 0.5 * _hr(1));
                double sigey = sigma_m * gauss_legendre_2D_cube((PML_ORD + 1), sigma, (void*)data, dx - 0.5 * _hr(0), dx + 0.5 * _hr(0),       -_hr(1),         0.0);
                double sig   = sigma_m * gauss_legendre_2D_cube((PML_ORD + 1), sigma, (void*)data, dx       - _hr(0), dx              ,       -_hr(1),         0.0);
                alpha_e.localElement(elem) = Vector_t(cosphi * sigey, sinphi * sigex);
                alpha_h.localElement(elem) = Vector_t(cosphi * sig  , sinphi * sig  );
            }
        }
    }
    {
        const double data[3] = {-sinphi, cosphi, (cosphi * _hr(0) + sinphi * _hr(1)) * PML_WIDTH};
        const double P[2] = {domain_length - sinphi * _vwidth - cosphi * data[2], sinphi * data[2]};
        const double sigma_m = -(PML_ORD + 1) * log(PML_RED) / (2.0 * Z_0 * data[2] * _hr(0) * _hr(1));
        int BBx1 = (int)floor(P[0] / _hr(0) + 0.99);
        int BBx2 = gDom[0].last();
        int BBy1 = 0;
        int BBy2 = (int)floor((cosphi * _vwidth + sinphi * data[2]) / _hr(1) + 0.99);
        int X1 = lDom[0].first() > BBx1? lDom[0].first(): BBx1;
        int X2 = lDom[0].last() < BBx2? lDom[0].last(): BBx2;
        int Y1 = lDom[1].first() > BBy1? lDom[1].first(): BBy1;
        int Y2 = lDom[1].last() < BBy2? lDom[1].last(): BBy2;

        for (int X = X1; X <= X2; X ++) {
            elem[0] = Index(X, X);
            for (int Y = Y1; Y <= Y2; Y ++) {
                elem[1] = Index(Y, Y);
                double dx = sinphi * (P[1] - Y * _hr(1)) / cosphi - (P[0] - X * _hr(0));
                double sigex = sigma_m * gauss_legendre_2D_cube((PML_ORD + 1), sigma, (void*)data, dx              , dx       + _hr(0), -0.5 * _hr(1), 0.5 * _hr(1));
                double sigey = sigma_m * gauss_legendre_2D_cube((PML_ORD + 1), sigma, (void*)data, dx - 0.5 * _hr(0), dx + 0.5 * _hr(0),       -_hr(1),         0.0);
                double sig   = sigma_m * gauss_legendre_2D_cube((PML_ORD + 1), sigma, (void*)data, dx              , dx       + _hr(0),       -_hr(1),         0.0);

                alpha_e.localElement(elem) = Vector_t(cosphi * sigey, sinphi * sigex);
                alpha_h.localElement(elem) = Vector_t(cosphi * sig  , sinphi * sig  );
            }
        }
    }

    alpha_e.fillGuardCells();
    alpha_h.fillGuardCells();
}

void Bend::getLambdaField(VField_Edge_t & lambdaField, std::vector<int> & bCellNrToLocal) {
    const NDIndex<DIM> & lDom = _geo.getLocalDomain();
    const auto& internalVertices = _geo.getInternalVertices();
    auto& insideMask = _geo.getInsideMask();
    insideMask = false;

    // drawDomain("Data/domain.ppm",
    //            internalVertices);

    NDIndex<DIM> elem, elempx, elempy, elempxpy;
    for (int cy = lDom[1].first(); cy <= lDom[1].last(); cy++) {
        elem[1] = elempx[1] = Index(cy,cy);
        elempy[1] = elempxpy[1] = Index(cy+1,cy+1);
        for (int cx = lDom[0].first(); cx <= lDom[0].last(); cx++) {
            elem[0] = elempy[0] = Index(cx,cx);
            elempx[0] = elempxpy[0] = Index(cx+1,cx+1);

            int ccase = internalVertices.localElement(elem) +
                (internalVertices.localElement(elempx) << 1) +
                (internalVertices.localElement(elempy) << 2) +
                (internalVertices.localElement(elempxpy) << 3);

            auto& lambda = lambdaField.localElement(elem);

            if (ccase == 15) {
                lambdaField.localElement(elem) = Vector_t(1.0);
                insideMask.localElement(elem) = true;
                continue;
            }
            if (ccase == 0) {
                continue;
            }

            size_t localID = getLocalID(elem);
            bCellNrToLocal.push_back(localID);

            switch (ccase) {
            case 1:
                {
                    lambda[0] = crossBoundary(Vector_t(cx * _hr(0), cy * _hr(1)), Vector_t(_hr(0), 0.0));
                    lambda[1] = crossBoundary(Vector_t(cx * _hr(0), cy * _hr(1)), Vector_t(0.0, _hr(1)));
                    break;
                }
            case 2:
                {
                    lambda[0] = crossBoundary(Vector_t((cx + 1) * _hr(0), cy * _hr(1)), Vector_t(-_hr(0), 0.0));
                    lambda[1] = 0.0;
                    break;
                }
            case 3:
                {
                    lambda[0] = 1.0;
                    lambda[1] = crossBoundary(Vector_t(cx * _hr(0), cy * _hr(1)), Vector_t(0.0, _hr(1)));
                    break;
                }
            case 4:
                {
                    lambda[0] = 0.0;
                    lambda[1] = crossBoundary(Vector_t(cx * _hr(0), (cy + 1) * _hr(1)), Vector_t(0.0, -_hr(1)));
                    break;
                }
            case 5:
                {
                    lambda[0] = crossBoundary(Vector_t(cx * _hr(0), cy * _hr(1)), Vector_t(_hr(0), 0.0));
                    lambda[1] = 1.0;
                    break;
                }
            case 7:
                {
                    lambda[0] = 1.0;
                    lambda[1] = 1.0;
                    break;
                }
            case 8:
                {
                    lambda[0] = 0.0;
                    lambda[1] = 0.0;
                    break;
                }
            case 10:
                {
                    lambda[0] = crossBoundary(Vector_t((cx + 1) * _hr(0), cy * _hr(1)), Vector_t(-_hr(0), 0.0));
                    lambda[1] = 0.0;
                    break;
                }
            case 11:
                {
                    lambda[0] = 1.0;
                    lambda[1] = crossBoundary(Vector_t(cx * _hr(0), cy * _hr(1)), Vector_t(0.0, _hr(1)));
                    break;
                }
            case 12:
                {
                    lambda[0] = 0.0;
                    lambda[1] = crossBoundary(Vector_t(cx * _hr(0), (cy + 1) * _hr(1)), Vector_t(0.0, -_hr(1)));
                    break;
                }
            case 13:
                {
                    lambda[0] = crossBoundary(Vector_t(cx * _hr(0), cy * _hr(1)), Vector_t(_hr(0), 0.0));
                    lambda[1] = 1.0;
                    break;
                }
            case 14:
                {
                    lambda[0] = crossBoundary(Vector_t((cx + 1) * _hr(0), cy * _hr(1)), Vector_t(-_hr(0), 0.0));
                    lambda[1] = crossBoundary(Vector_t(cx * _hr(0), (cy + 1) * _hr(1)), Vector_t(0.0, -_hr(1)));
                    break;
                }
            default:
                // nothing to do since eigther all vertices are outside
                continue;
            }
        }
    }

    lambdaField.fillGuardCells();
    insideMask.fillGuardCells();
}

void Bend::getLambdaFieldDualGrid(VField_Edge_t & lambdaField, std::vector<int> & bCellNrToLocal) {
    const NDIndex<DIM> & lDom = _geo.getLocalDomain();
    const auto& internalVertices = _geo.getInternalVerticesDualGrid();
    auto& insideMask = _geo.getInsideMaskDualGrid();
    insideMask = false;

    NDIndex<DIM> elem, elempx, elempy, elempxpy;
    for (int cy = lDom[1].first(); cy <= lDom[1].last(); cy++) {
        elem[1] = elempx[1] = Index(cy,cy);
        elempy[1] = elempxpy[1] = Index(cy+1,cy+1);
        for (int cx = lDom[0].first(); cx <= lDom[0].last(); cx++) {
            elem[0] = elempy[0] = Index(cx,cx);
            elempx[0] = elempxpy[0] = Index(cx+1,cx+1);

            int ccase = internalVertices.localElement(elem) +
                (internalVertices.localElement(elempx) << 1) +
                (internalVertices.localElement(elempy) << 2) +
                (internalVertices.localElement(elempxpy) << 3);

            auto& lambda = lambdaField.localElement(elem);

            if (ccase == 15) {
                lambdaField.localElement(elem) = Vector_t(1.0);
                insideMask.localElement(elem) = true;
                continue;
            }
            if (ccase == 0) {
                continue;
            }

            size_t localID = getLocalID(elem);
            bCellNrToLocal.push_back(localID);

            switch (ccase) {
            case 1:
                {
                    lambda[0] = crossBoundary(Vector_t((cx + 0.5) * _hr(0), (cy + 0.5) * _hr(1)), Vector_t(_hr(0), 0.0));
                    lambda[1] = crossBoundary(Vector_t((cx + 0.5) * _hr(0), (cy + 0.5) * _hr(1)), Vector_t(0.0, _hr(1)));
                    break;
                }
            case 2:
                {
                    lambda[0] = crossBoundary(Vector_t((cx + 1.5) * _hr(0), (cy + 0.5) * _hr(1)), Vector_t(-_hr(0), 0.0));
                    lambda[1] = 0.0;
                    break;
                }
            case 3:
                {
                    lambda[0] = 1.0;
                    lambda[1] = crossBoundary(Vector_t((cx + 0.5) * _hr(0), (cy + 0.5) * _hr(1)), Vector_t(0.0, _hr(1)));
                    break;
                }
            case 4:
                {
                    lambda[0] = 0.0;
                    lambda[1] = crossBoundary(Vector_t((cx + 0.5) * _hr(0), (cy + 1.5) * _hr(1)), Vector_t(0.0, -_hr(1)));
                    break;
                }
            case 5:
                {
                    lambda[0] = crossBoundary(Vector_t((cx + 0.5) * _hr(0), (cy + 0.5) * _hr(1)), Vector_t(_hr(0), 0.0));
                    lambda[1] = 1.0;
                    break;
                }
            case 7:
                {
                    lambda[0] = 1.0;
                    lambda[1] = 1.0;
                    break;
                }
            case 8:
                {
                    lambda[0] = 0.0;
                    lambda[1] = 0.0;
                    break;
                }
            case 10:
                {
                    lambda[0] = crossBoundary(Vector_t((cx + 1.5) * _hr(0), (cy + 0.5) * _hr(1)), Vector_t(-_hr(0), 0.0));
                    lambda[1] = 0.0;
                    break;
                }
            case 11:
                {
                    lambda[0] = 1.0;
                    lambda[1] = crossBoundary(Vector_t((cx + 0.5) * _hr(0), (cy + 0.5) * _hr(1)), Vector_t(0.0, _hr(1)));
                    break;
                }
            case 12:
                {
                    lambda[0] = 0.0;
                    lambda[1] = crossBoundary(Vector_t((cx + 0.5) * _hr(0), (cy + 1.5) * _hr(1)), Vector_t(0.0, -_hr(1)));
                    break;
                }
            case 13:
                {
                    lambda[0] = crossBoundary(Vector_t((cx + 0.5) * _hr(0), (cy + 0.5) * _hr(1)), Vector_t(_hr(0), 0.0));
                    lambda[1] = 1.0;
                    break;
                }
            case 14:
                {
                    lambda[0] = crossBoundary(Vector_t((cx + 1.5) * _hr(0), (cy + 0.5) * _hr(1)), Vector_t(-_hr(0), 0.0));
                    lambda[1] = crossBoundary(Vector_t((cx + 0.5) * _hr(0), (cy + 1.5) * _hr(1)), Vector_t(0.0, -_hr(1)));
                    break;
                }
            default:
                // nothing to do since eigther all vertices are outside
                continue;
            }
        }
    }

    lambdaField.fillGuardCells();
    insideMask.fillGuardCells();
}

void Bend::getBoundaryCellAreas(std::vector<BoundaryCell> & boundaryCells,
                                std::vector<size_t> & toBeDeleted,
                                const VField_Edge_t & lambdaField,
                                const std::vector<int> & bCellNrToLocal) {
    const auto lDom = _geo.getLocalDomain();
    const auto& internalVertices = _geo.getInternalVertices();

    NDIndex<DIM> elem, elempx, elempy, elempxpy;
    for (const size_t& localID: bCellNrToLocal) {
        elem = elempx = elempy = elempxpy = getNDIndex(localID);
        elempx[0] = elempxpy[0] = elem[0] + 1;
        elempy[1] = elempxpy[1] = elem[1] + 1;

        int ccase = internalVertices.localElement(elem) +
            (internalVertices.localElement(elempx) << 1) +
            (internalVertices.localElement(elempy) << 2) +
            (internalVertices.localElement(elempxpy) << 3);

        Vector_t lambda = lambdaField.localElement(elem);
        Vector_t lambdaPx = lambdaField.localElement(elempx);
        Vector_t lambdaPy = lambdaField.localElement(elempy);
        Vector_t lambdaNeighbours(lambdaPy[0], lambdaPx[1]);

        double area = areaCases(ccase, lambda, lambdaNeighbours);
        PAssert(area < 1.0 + SPACE_EPS);
        PAssert(area >= 0.0);
        if (area < SPACE_EPS) {
            toBeDeleted.push_back(localID);
        } else {
            boundaryCells.push_back({localID, lambda, area, lambdaNeighbours});
        }
    }
}

void Bend::getBoundaryCellAreasDualGrid(std::vector<BoundaryCell> & boundaryCells,
                                        std::vector<size_t> & toBeDeleted,
                                        const VField_Edge_t & lambdaField,
                                        const std::vector<int> & bCellNrToLocal) {
    const auto lDom = _geo.getLocalDomain();
    const auto& internalVertices = _geo.getInternalVerticesDualGrid();

    NDIndex<DIM> elem, elempx, elempy, elempxpy;
    for (const size_t& localID: bCellNrToLocal) {
        elem = elempx = elempy = elempxpy = getNDIndex(localID);
        elempx[0] = elempxpy[0] = elem[0] + 1;
        elempy[1] = elempxpy[1] = elem[1] + 1;

        int ccase = internalVertices.localElement(elem) +
            (internalVertices.localElement(elempx) << 1) +
            (internalVertices.localElement(elempy) << 2) +
            (internalVertices.localElement(elempxpy) << 3);

        Vector_t lambda = lambdaField.localElement(elem);
        Vector_t lambdaPx = lambdaField.localElement(elempx);
        Vector_t lambdaPy = lambdaField.localElement(elempy);
        Vector_t lambdaNeighbours(lambdaPy[0], lambdaPx[1]);

        double area = areaCases(ccase, lambda, lambdaNeighbours);

        if (area < SPACE_EPS * SPACE_EPS) {
            toBeDeleted.push_back(localID);
        } else {
            boundaryCells.push_back({localID, lambda, area, lambdaNeighbours});
        }
    }
}

double Bend::areaCases(const int & ccase, const Vector_t & lambda, const Vector_t & lambdaNeighbours) {
    double area = 0.0;
    switch (ccase) {
    case 1:
        {
            area = lambda[0] * lambda[1] / 2;
            break;
        }
    case 2:
        {
            area = lambda[0] * lambdaNeighbours[1] / 2;
            break;
        }
    case 3:
        {
            area = (lambda[1] + lambdaNeighbours[1]) / 2;
            break;
        }
    case 4:
        {
            area = lambdaNeighbours[0] * lambda[1] / 2;
            break;
        }
    case 5:
        {
            area = (lambda[0] + lambdaNeighbours[0]) / 2;
            break;
        }
    case 7:
        {
            area = 1 - (1 - lambdaNeighbours[0]) * (1 - lambdaNeighbours[1]) / 2;
            break;
        }
    case 8:
        {
            area = lambdaNeighbours[0] * lambdaNeighbours[1] / 2;
            break;
        }
    case 10:
        {
            area = (lambda[0] + lambdaNeighbours[0]) / 2;
            break;
        }
    case 11:
        {
            area = 1 - (1 - lambda[1]) * (1 - lambdaNeighbours[0]) / 2;
            break;
        }
    case 12:
        {
            area = (lambda[1] + lambdaNeighbours[1]) / 2;
            break;
        }
    case 13:
        {
            area = 1 - (1 - lambda[0]) * (1 - lambdaNeighbours[1]) / 2;
            break;
        }
    case 14:
        {
            area = 1 - (1 - lambda[0]) * (1 - lambda[1]) / 2;
            break;
        }
    default:
        // nothing to do since eigther all vertices are outside
        ;
    }
    return area;
}

void Bend::getBoundaryCellAreasGhostCells(const VField_Edge_t & lambda) {
    const NDIndex<DIM> & lDom = _geo.getLocalDomain();
    const NDIndex<DIM> & gDom = _geo.getGlobalDomain();

    NDIndex<DIM> elem;
    if (lDom[0].first() > gDom[0].first()) {
        elem[0] = Index(lDom[0].first() - 1,lDom[0].first() - 1);
        elem[1] = lDom[1];
        getBoundaryCellAreasGhostCellsImpl(lambda,elem);
    }
    if (lDom[0].last() < gDom[0].last()) {
        elem[0] = Index(lDom[0].last() + 1,lDom[0].last() + 1);
        elem[1] = lDom[1];
        getBoundaryCellAreasGhostCellsImpl(lambda,elem);
    }
    if (lDom[1].first() > gDom[1].first()) {
        elem[0] = lDom[0];
        elem[1] = Index(lDom[1].first() - 1,lDom[1].first() - 1);
        getBoundaryCellAreasGhostCellsImpl(lambda,elem);
    }
    if (lDom[1].last() < gDom[1].last()) {
        elem[0] = lDom[0];
        elem[1] = Index(lDom[1].last() + 1,lDom[1].last() + 1);
        getBoundaryCellAreasGhostCellsImpl(lambda,elem);
    }
}

void Bend::getBoundaryCellAreasGhostCellsImpl(const VField_Edge_t & lambdaField, const NDIndex<DIM> & dom) {
    const auto& insideMask = _geo.getInsideMask();
    const auto& internalVertices = _geo.getInternalVertices();
    auto& ghostCellBoundaryCells = _geo.getGCBoundaryCells();

    NDIndex<DIM> elem, elempx, elempy, elempxpy;
    for (int cy = dom[1].first(); cy <= dom[1].last(); ++ cy) {
        elem[1] = elempx[1] = Index(cy,cy);
        elempy[1] = elempxpy[1] = Index(cy+1,cy+1);
        for (int cx = dom[0].first(); cx <= dom[0].last(); ++ cx) {
            elem[0] = Index(cx,cx);

            if (insideMask.localElement(elem)) continue;

            elempx[0] = elempxpy[0] = Index(cx+1,cx+1);
            elempy[0] = elem[0];

            int ccase = internalVertices.localElement(elem) +
                (internalVertices.localElement(elempx) << 1) +
                (internalVertices.localElement(elempy) << 2) +
                (internalVertices.localElement(elempxpy) << 3);

            if (ccase == 0) continue;

            Vector_t lambda = lambdaField.localElement(elem);
            Vector_t lambdaPx = lambdaField.localElement(elempx);
            Vector_t lambdaPy = lambdaField.localElement(elempy);
            Vector_t lambdaNeighbours(lambdaPy[0], lambdaPx[1]);

            double area = areaCases(ccase, lambda, lambdaNeighbours);
            size_t localID = getLocalID(elem);

            ghostCellBoundaryCells.push_back({localID, lambda, area, lambdaNeighbours});
        }
    }
}

std::vector<double> Bend::getExternalField(const ParticleAttrib< Vector_t > & R) const
{
    typedef ParticleAttrib< Vector_t >::const_iterator const_iterator;
    Tenzor<double, DIM> M(cos(_phi), -sin(_phi), sin(_phi), cos(_phi));
    std::vector<double> bf;
    const double cosphihalf = cos(0.5 * _phi);
    const double sinphihalf = sin(0.5 * _phi);
    const double halfWidth = 0.5 * _vwidth;
    const double bendWidth = 2 * _R * std::abs(sinphihalf);

    // DBGOUT << R[0](0) - _dlengthBefore << "\t" << R[0](1) - halfWidth << "\t"
    //        << cosphihalf * (R[0](0) - _dlengthBefore) + sinphihalf * (R[0](1) - halfWidth) << std::endl;
    for (const_iterator r = R.cbegin(); r != R.cend(); ++ r) {
        Vector_t Dx = *r - Vector_t(_dlengthBefore, halfWidth);
        double distance = cosphihalf * Dx(0) - sinphihalf * Dx(1);
        if (distance > 0.0 && distance < bendWidth) {
            bf.push_back(_bz);
        } else {
            bf.push_back(0.0);
        }

        // Vector_t Dx = *r - _center;
        // Vector_t Dxrot = dot(M, Dx);
        // if (Dx(0) < -_outerCutoffDistance ||
        //     Dxrot(0) > _outerCutoffDistance) {
        //     bf.push_back(0.0);
        //     continue;
        // }

        // if (Dx(0) > _innerCutoffDistance &&
        //     Dxrot(0) < -_innerCutoffDistance) {
        //     bf.push_back(_bz);
        //     continue;
        // }

        // if (Dx(0) < _innerCutoffDistance) {
        //     bf.push_back(_bz);
        //     continue;
        // }
    }

    return bf;
}

void Bend::drawDomain(const std::string & filename,
                      const Field<char, DIM, Mesh_t, Vert> & values) const {
    // This routine is based on the official solution of an assignment in Informatik I
    // at ETH Zurich written by Cyril Flaig
    typedef Utils::color_rgb color_rgb;
    const auto & insideMask = values;

    size_t pos = filename.find_last_of('.');
    char myFilename[100];
    sprintf(myFilename, "%s%03d.%s", filename.substr(0, pos + 1).c_str(), Ippl::myNode(), filename.substr(pos + 1).c_str());

    const NDIndex<DIM>& lDom = values.getLayout().getLocalNDIndex();
    int Nx = lDom[0].length();
    int Ny = lDom[1].length();

    std::vector<color_rgb> pic(Ny*Nx);
    color_rgb white(1,1,1);
    color_rgb black(0,0,0);
    color_rgb red(1,0,0);
    color_rgb green(0,1,0);
    color_rgb yellow(1,1,0);

    std::fill(pic.begin(), pic.end(), black);

    NDIndex<DIM> elem;
    for (int cy = lDom[1].first(); cy <= lDom[1].last(); ++ cy) {
        int ny = cy - lDom[1].first();
        elem[1] = Index(cy, cy);
        for (int cx = lDom[0].first(); cx <= lDom[0].last(); ++ cx) {
            int nx = cx - lDom[0].first();
            elem[0] = Index(cx, cx);
            if (insideMask.localElement(elem))
                pic[ny * Nx + nx] = yellow;
        }
    }

    std::ofstream outfile(myFilename, std::ios::binary);

    std::stringstream header;
    // Definiert den Header der PPM-Datei
    header << "P6" << std::endl;
    header << Nx << std::endl;
    header << Ny << std::endl;
    header << 255;
    outfile << header.str() << std::endl;

    // Fuer jeden Pixel werden die Farben geschrieben
    for (int i = 0; i < Ny; i++)
        for (int j = 0; j < Nx; j++)
            for (int l = 0; l < 3; l++)
                outfile << (unsigned char)(255*pic[(Ny - i-1)*Nx +j][l]);
}
