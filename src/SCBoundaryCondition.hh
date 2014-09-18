/***************************************************************************
                        SCBoundaryCondition.hh
                         -------------------
    begin                : Thu Jun 7 2012
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

#ifndef SCBOUNDARYCONDITION_H
#define SCBOUNDARYCONDITION_H

#define MINEDGELENGTH 0.00001

#include "defs.hh"
#include "utils.hh"

#include <vector>
#include <map>
#include <iostream>
#include "BoundaryCell.hh"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
/** \brief implements the SC boundary condition as in
 *  <a href="http://dx.doi.org/10.1016/j.jcp.2007.02.002">Conformal FDTD-methods to
 *  avoid time step reduction with and without cell enlargement</a> by
 *  I. Zagorodnov, R. Schuhmann and T. Weiland
 */
template<class G, class S>
class SCBoundaryCondition {
public:
    SCBoundaryCondition(const double & sratio, G & geo);

    void drawBoundary(const std::string& filename, const G& geo);
    void drawBoundaryDualGrid(const std::string & filename, const G& geo);

    void updateHz(VField_Cell_t& Hz,
                  const VField_Edge_t& EFD,
                  const S& solver,
                  const G& geo) const;

    void updateEx(VField_Edge_t& EFD,
                  const VField_Cell_t& HFD,
                  const VField_Edge_t& JFD,
                  const S& solver,
                  const G& geo) const;

    void updateEy(VField_Edge_t& EFD,
                  const VField_Cell_t& HFD,
                  const VField_Edge_t& JFD,
                  const S& solver,
                  const G& geo) const;

    const std::vector<BoundaryCell>& getBoundaryCells() const;
    const std::vector<BoundaryCell>& getBoundaryCellsDualGrid() const;

    void resetOutsideDomain(VField_Edge_t&EFD,
                            const G& geo) const;
private:
    void fillCellProperties(G& geo);

    void collectBoundaryCells(const G& geo,
                              std::vector<BoundaryCell>& bcells);

    void limitArea(BoundaryCell& cell) const;
    void limitEdges(BoundaryCell& cell,
                    const G& geo) const;

    std::vector<BoundaryCell> boundaryCells_m;
    std::vector<BoundaryCell> boundaryCellsDualGrid_m;
    std::vector<double> sumColumns_m;

    double stableRatio_m;
};

template<class G, class S>
SCBoundaryCondition<G,S>::SCBoundaryCondition(const double & sratio, G & geo):
    stableRatio_m(2.00)
{
    fillCellProperties(geo);

    // drawBoundary("Data/Boundary.ppm", geo);
}

template<class G, class S>
void SCBoundaryCondition<G,S>::fillCellProperties(G& geo) {
    // extern std::ofstream dbg;
    geo.fillCellProperties(boundaryCells_m,
                           boundaryCellsDualGrid_m);

    for(auto& cell: boundaryCells_m) {
        // limitEdges(cell, geo);
        limitArea(cell);
    }

    // std::string filename("Data/BoundaryCells.dat");
    // size_t pos = filename.find_last_of('.');
    // char myFilename[100];
    // sprintf(myFilename, "%s%03d.%s", filename.substr(0, pos + 1).c_str(), Ippl::myNode(), filename.substr(pos + 1).c_str());

    // std::ofstream out(myFilename);
    // for(auto& cell: boundaryCells_m) {
    //     Coordinates coord = geo.getCoordinates(cell.localID_m);
    //     cell.print(out, std::make_pair(coord.idx_m, coord.idy_m));
    // }
    // out.close();

    // filename = "Data/BoundaryCellsDualGrid.dat";
    // pos = filename.find_last_of('.');
    // sprintf(myFilename, "%s%03d.%s", filename.substr(0, pos + 1).c_str(), Ippl::myNode(), filename.substr(pos + 1).c_str());

    // out.open(myFilename);
    // for(auto& cell: boundaryCellsDualGrid_m) {
    //     Coordinates coord = geo.getCoordinates(cell.localID_m);
    //     cell.print(out, std::make_pair(coord.idx_m, coord.idy_m));
    // }
    // out.close();
}

template<class G, class S>
void SCBoundaryCondition<G,S>::limitArea(BoundaryCell& cell) const
{
    cell.area_m = std::max(cell.area_m, cell.lambda_m[0] / 2);
    cell.area_m = std::max(cell.area_m, cell.lambda_m[1]  / 2);
    cell.area_m = std::max(cell.area_m, cell.lambdaToNeighbours_m[0] / 2);
    cell.area_m = std::max(cell.area_m, cell.lambdaToNeighbours_m[1] / 2);
}

template<class G, class S>
void SCBoundaryCondition<G,S>::limitEdges(BoundaryCell& cell,
                                          const G& geo) const
{
    double areaSouth = std::min(cell.area_m,
                                geo.getArea(geo.getNeighbourSouth(cell.localID_m),
                                            boundaryCells_m));
    double areaWest = std::min(cell.area_m,
                               geo.getArea(geo.getNeighbourWest(cell.localID_m),
                                           boundaryCells_m));
    double areaNorth = std::min(cell.area_m,
                                geo.getArea(geo.getNeighbourNorth(cell.localID_m),
                                            boundaryCells_m));
    double areaEast = std::min(cell.area_m,
                               geo.getArea(geo.getNeighbourEast(cell.localID_m),
                                           boundaryCells_m));
    cell.lambda_m[0] = std::min(cell.lambda_m[0], 2 * areaSouth);
    cell.lambda_m[1] = std::min(cell.lambda_m[1], 2 * areaWest);
    cell.lambdaToNeighbours_m[0] = std::min(cell.lambdaToNeighbours_m[0],
                                            2 * areaNorth);
    cell.lambdaToNeighbours_m[1] = std::min(cell.lambdaToNeighbours_m[1],
                                            2 * areaEast);
}

template<class G, class S>
void SCBoundaryCondition<G,S>::resetOutsideDomain(VField_Edge_t&EFD,
                                                  const G& geo) const
{
    const Field<bool, DIM>& insideMask = geo.getInsideMask();
    const std::vector<int>& localToBCellNr = geo.getInverseMap();
    const NDIndex<DIM>& lDom = geo.getLocalDomain();
    NDIndex<DIM> elem;

    for (int cy = lDom[1].first(); cy <= lDom[1].last(); ++ cy) {
        elem[1] = Index(cy, cy);
        for (int cx = lDom[0].first(); cx <= lDom[0].last(); ++ cx) {
            elem[0] = Index(cx, cx);
            if (insideMask.localElement(elem)) continue;

            size_t localID = geo.getLocalID(elem);
            int bCellNr = localToBCellNr[localID];
            if (bCellNr == -1) {
                EFD.localElement(elem) = Vector_t(0.0);
            } else {
                if (boundaryCells_m[bCellNr].lambda_m[0] < SPACE_EPS) {
                    EFD.localElement(elem)[0] = 0.0;
                }
                if (boundaryCells_m[bCellNr].lambda_m[1] < SPACE_EPS) {
                    EFD.localElement(elem)[1] = 0.0;
                }
                if (boundaryCells_m[bCellNr].area_m < MINAREA) {
                    EFD.localElement(elem) = Vector_t(0.0);
                }
            }
        }
    }
}

template<class G, class S>
void SCBoundaryCondition<G,S>::collectBoundaryCells(const G& geo,
                                                    std::vector<BoundaryCell>& bcells) {
    extern std::ofstream dbg;

    char fname[100];
    sprintf(fname, "Data/serializedBoundaryCells.%03d.dat", Ippl::myNode());
    {
        std::ofstream outf(fname);
        boost::archive::text_oarchive oa(outf);
        oa << boundaryCells_m;
    }

    int* limits = new int[4 * Ippl::getNodes()];
    NDIndex<DIM> lDom = geo.getLocalDomain();
    limits[4*Ippl::myNode()] = lDom[0].first();
    limits[4*Ippl::myNode() + 1] = lDom[0].last();
    limits[4*Ippl::myNode() + 2] = lDom[1].first();
    limits[4*Ippl::myNode() + 3] = lDom[1].last();

    if (Ippl::getNodes() > 1) {
        dbg << "SCBoundaryCondition.hh: " << __LINE__ << "using collective communication" << std::endl;
        MPI_Gather(&(limits[4*Ippl::myNode()]),
                   4,
                   MPI_INT,
                   limits,
                   4,
                   MPI_INT,
                   0,
                   MPI_COMM_WORLD);
    }

    std::vector<BoundaryCell> tmpCells;
    if (Ippl::myNode() == 0) {
        NDIndex<DIM> gDom = geo.getGlobalDomain();
        tmpCells.insert(bcells.end(), boundaryCells_m.begin(), boundaryCells_m.end());
        for (auto& cell: tmpCells) {
            int idx = cell.localID_m % (lDom[0].length() + 4) + lDom[0].first() - gDom[0].first();
            int idy = cell.localID_m / (lDom[0].length() + 4) + lDom[1].first() - gDom[1].first();
            cell.localID_m = idx + idy * (gDom[0].length() + 4);
        }
        bcells.insert(bcells.end(), tmpCells.begin(), tmpCells.end());

        for (int i = 1; i < Ippl::getNodes(); ++ i) {
            sprintf(fname, "Data/serializedBoundaryCells.%03d.dat", i);
            std::ifstream inf(fname);
            if (inf.good()) {
                boost::archive::text_iarchive ia(inf);
                ia >> tmpCells;
            }
            NDIndex<DIM> rDom(Index(limits[4 * i], limits[4 * i + 1]),
                              Index(limits[4 * i + 2], limits[4 * i + 3]));
            for (auto& cell: tmpCells) {
                int idx = cell.localID_m % (rDom[0].length() + 4) + rDom[0].first() - gDom[0].first();
                int idy = cell.localID_m / (rDom[0].length() + 4) + rDom[1].first() - gDom[1].first();
                cell.localID_m = idx + idy * (gDom[0].length() + 4);
            }
            bcells.insert(bcells.end(), tmpCells.begin(), tmpCells.end());
        }

        std::sort(bcells.begin(), bcells.end(),
                  [](const BoundaryCell& a, const BoundaryCell& b) ->bool {
                      return a.localID_m < b.localID_m;
                  });
    }

    delete[] limits;
}

template<class G, class S>
void SCBoundaryCondition<G,S>::drawBoundary(const std::string & filename,
                                            const G& geo) {
    // This routine is based on the official solution of an assignment in Informatik I
    // at ETH Zurich written by Cyril Flaig
    typedef Utils::color_rgb color_rgb;
    const Field<bool, DIM>& insideMask = geo.getInsideMask();

    size_t pos = filename.find_last_of('.');
    char myFilename[100];
    sprintf(myFilename, "%s%03d.%s", filename.substr(0, pos + 1).c_str(), Ippl::myNode(), filename.substr(pos + 1).c_str());

    const NDIndex<DIM>& lDom = geo.getLocalDomain();
    int Nx = lDom[0].length();
    int Ny = lDom[1].length();

    std::vector<color_rgb> pic(Ny*Nx);
    color_rgb white(1,1,1);
    color_rgb black(0,0,0);
    color_rgb red(1,0,0);
    color_rgb green(0,1,0);
    color_rgb yellow(1,1,0);

    std::fill(pic.begin(), pic.end(), white);

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

    for (const auto& cell: boundaryCells_m) {
        Coordinates coord = geo.getCoordinates(cell.localID_m);
        int nx = coord.idx_m - lDom[0].first();
        int ny = coord.idy_m - lDom[1].first();

        pic[ny * Nx + nx] += green;
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

template<class G, class S>
void SCBoundaryCondition<G,S>::drawBoundaryDualGrid(const std::string & filename,
                                                    const G& geo) {
    // This routine is based on the official solution of an assignment in Informatik I
    // at ETH Zurich written by Cyril Flaig
    typedef Utils::color_rgb color_rgb;

    const auto & insideMask = geo.getInsideMaskDualGrid();

    size_t pos = filename.find_last_of('.');
    char myFilename[100];
    sprintf(myFilename, "%s%03d.%s", filename.substr(0, pos + 1).c_str(), Ippl::myNode(), filename.substr(pos + 1).c_str());

    const NDIndex<DIM>& lDom = geo.getLocalDomain();
    int Nx = lDom[0].length();
    int Ny = lDom[1].length();

    std::vector<color_rgb> pic(Ny*Nx);
    color_rgb white(1,1,1);
    color_rgb black(0,0,0);
    color_rgb red(1,0,0);
    color_rgb green(0,1,0);
    color_rgb yellow(1,1,0);

    std::fill(pic.begin(), pic.end(), white);

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

    for (const auto& cell: boundaryCellsDualGrid_m) {
        Coordinates coord = geo.getCoordinates(cell.localID_m);
        int nx = coord.idx_m - lDom[0].first();
        int ny = coord.idy_m - lDom[1].first();

        pic[ny * Nx + nx] += green;
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

template<class G, class S>
void SCBoundaryCondition<G,S>::updateHz(VField_Cell_t& Hz,
                                        const VField_Edge_t& EFD,
                                        const S& solver,
                                        const G& geo) const {
    for (const auto& cell: boundaryCells_m) {
        NDIndex<DIM> elem = geo.getNDIndex(cell.localID_m);
        Hz.localElement(elem) += solver.updateHz(EFD, elem, cell);
    }
}

template<class G, class S>
void SCBoundaryCondition<G,S>::updateEx(VField_Edge_t& EFD,
                                        const VField_Cell_t& HFD,
                                        const VField_Edge_t& JFD,
                                        const S& solver,
                                        const G& geo) const {
    for (const auto& cell: boundaryCells_m) {
        if (cell.area_m < MINAREA) continue;

        NDIndex<DIM> elem = geo.getNDIndex(cell.localID_m);
        if (cell.lambda_m(0) > SPACE_EPS)
            EFD.localElement(elem)(0) += solver.updateEx(HFD, JFD, elem);
    }
}

template<class G, class S>
void SCBoundaryCondition<G,S>::updateEy(VField_Edge_t& EFD,
                                        const VField_Cell_t& HFD,
                                        const VField_Edge_t& JFD,
                                        const S& solver,
                                        const G& geo) const {

    for (const auto& cell: boundaryCells_m) {
        if (cell.area_m < MINAREA) continue;

        NDIndex<DIM> elem = geo.getNDIndex(cell.localID_m);
        if (cell.lambda_m(1) > SPACE_EPS)
            EFD.localElement(elem)(1) += solver.updateEy(HFD, JFD, elem);
    }
}

template<class G, class S>
const std::vector<BoundaryCell>& SCBoundaryCondition<G,S>::getBoundaryCells() const {
    return boundaryCells_m;
}

template<class G, class S>
const std::vector<BoundaryCell>& SCBoundaryCondition<G,S>::getBoundaryCellsDualGrid() const {
    return boundaryCellsDualGrid_m;
}

#endif
