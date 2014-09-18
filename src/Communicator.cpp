/***************************************************************************
                           Communicator.cpp
                         -------------------
    begin                : Sat Sep 21 2013
    copyright            : (C) 2013 by Christof Kraus
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

#include "Communicator.hh"
#include "utils.hh"

extern std::ofstream dbg;

#define DBGOUT dbg << "Communicator.cpp: " << __LINE__ << "\t"

std::vector<unsigned int> Communicator::_index2node;
std::vector<unsigned int> Communicator::_node2index;
std::vector<unsigned int> Communicator::_neighbours;
unsigned int Communicator::_nNodeCols = 0;
unsigned int Communicator::_nNodeRows = 0;
bool Communicator::_initialized = false;

void Communicator::initialize()
{
    if (_initialized) return;

    _index2node.resize(Ippl::getNodes());
    _node2index.resize(Ippl::getNodes());

    unsigned int depth = static_cast<unsigned int>(floor(log(1.0 * Ippl::getNodes()) / log(2.0) + 0.5));
    unsigned int numNodes = Utils::ipow(2, depth);
    if (numNodes != static_cast<unsigned int>(Ippl::getNodes()) && Ippl::myNode() == 0) {
        Ippl::exitAllNodes("\033[31;1mnumber of cores not a power of 2!!\033[0m", true);
    }
    unsigned int coldepth = static_cast<unsigned int>(floor(0.5 * (depth + 1.0)));
    unsigned int rowdepth = static_cast<unsigned int>(floor(0.5 * depth));
    _nNodeCols = Utils::ipow(2, coldepth);
    _nNodeRows = Utils::ipow(2, rowdepth);
    std::vector<std::pair<unsigned int, unsigned int> > n2mi(numNodes, std::make_pair(0,0));
    unsigned int nNodes = 1;
    for (unsigned int m = 1; m <= depth; ++ m) {
        for (unsigned int l = nNodes; l > 0; -- l) {
            if (m % 2 == 0) {
                coldepth = static_cast<unsigned int>(floor(0.5 * (m + 1.0)));
                unsigned int curcol = Utils::ipow(2, coldepth);
                n2mi[2*(l-1)].second = n2mi[l-1].second + (n2mi[l-1].second / curcol) * curcol;
                n2mi[2*(l-1)+1].second = n2mi[2*(l-1)].second + curcol;
            } else {
                n2mi[2*(l-1)].second = 2 * n2mi[l-1].second;
                n2mi[2*(l-1)+1].second = n2mi[2*(l-1)].second + 1;
            }
        }
        nNodes *= 2;
    }
    for (unsigned int i = 0; i < numNodes; ++ i) {
        n2mi[i].first = i;
        _node2index[i] = n2mi[i].second;
    }
    std::sort(n2mi.begin(), n2mi.end(),
              [](const std::pair<unsigned int, unsigned int> & a,
                 const std::pair<unsigned int, unsigned int> & b) { return a.second < b.second;});

    for (unsigned int i = 0; i < numNodes; ++ i) {
        _index2node[i] = n2mi[i].first;
    }

    unsigned int index = _node2index[Ippl::myNode()];
    unsigned int myCol = index % _nNodeCols;
    unsigned int myRow = index / _nNodeCols;
    unsigned int startCol = std::min(myCol, myCol - 1);
    unsigned int endCol   = std::min(_nNodeCols, myCol + 2);
    unsigned int startRow = std::min(myRow, myRow - 1);
    unsigned int endRow   = std::min(_nNodeRows, myRow + 2);

    for (unsigned int j = startRow; j < endRow; ++ j) {
        for (unsigned int i = startCol; i < endCol; ++ i) {
            if (i == myCol && j == myRow) continue;
            index = j * _nNodeCols + i;
            _neighbours.push_back(_index2node[index]);
        }
    }
    std::sort(_neighbours.begin(), _neighbours.end());

    _initialized = true;
}

std::vector<NDIndex<DIM> > Communicator::getAllLocalDomains(NDIndex<DIM> lDom)
{
    unsigned int myNode = Ippl::myNode();
    unsigned int numNodes = Ippl::getNodes();
    int *values = new int[2 * DIM * numNodes];
    for (unsigned int d = 0; d < DIM; ++ d) {
        values[2 * (myNode * DIM + d)    ] = lDom[d].first();
        values[2 * (myNode * DIM + d) + 1] = lDom[d].last();
    }
    DBGOUT << "using collective communication" << std::endl;
    MPI_Allgather(MPI_IN_PLACE, 2 * DIM, MPI_INT, values, 2 * DIM, MPI_INT, MPI_COMM_WORLD);

    std::vector<NDIndex<DIM> > ret(numNodes);
    for (unsigned int k = 0; k < numNodes; ++ k) {
        for (unsigned int d = 0; d < DIM; ++ d) {
            int lower = values[2 * (k * DIM + d)    ];
            int upper = values[2 * (k * DIM + d) + 1];
            int sign = lower <= upper? 1: upper - lower;
            ret[k][d] = Index(lower, upper, sign);
        }
    }

    delete[] values;

    return ret;
}

void Communicator::exchangeBoundariesMovingMesh(VField_Edge_t & surroundingE,
                                                VField_Cell_t & surroundingH,
                                                VField_Edge_t & embeddedE,
                                                VField_Cell_t & embeddedH,
                                                const Vektor<int, DIM> & positionEmbedded,
                                                const Vektor<int, DIM> & move)
{
    const unsigned int myNode = Ippl::myNode();
    const unsigned int numNodes = Ippl::getNodes();

    NDIndex<DIM> elem;
    NDIndex<DIM> gSurroundingEDomain = surroundingE.getLayout().getDomain();
    NDIndex<DIM> gSurroundingHDomain = surroundingH.getLayout().getDomain();
    NDIndex<DIM> gEmbeddedEDomain = embeddedE.getLayout().getDomain();
    NDIndex<DIM> gEmbeddedHDomain = embeddedH.getLayout().getDomain();

    std::vector<NDIndex<DIM> > ExS2EDoms, ExE2SDoms;
    std::vector<NDIndex<DIM> > EyS2EDoms, EyE2SDoms;
    std::vector<NDIndex<DIM> > HzS2EDoms, HzE2SDoms;

    std::vector<NDIndex<DIM> > localEmbeddedEDomains;
    std::vector<NDIndex<DIM> > localSurroundingEDomains;
    Utils::getLocalDomains(embeddedE.getLayout(), localEmbeddedEDomains);
    Utils::getLocalDomains(surroundingE.getLayout(), localSurroundingEDomains);
    for (unsigned int k = 0; k < numNodes; ++ k) {
        NDIndex<DIM> & dom = localEmbeddedEDomains[k];
        for (unsigned int d = 0; d < DIM; ++ d) {
            if (dom[d].first() == gEmbeddedEDomain[d].first())
                dom[d] = Index(dom[d].first() - 1, dom[d].last());
            if (dom[d].last() == gEmbeddedEDomain[d].last())
                dom[d] = Index(dom[d].first(), dom[d].last() + 1);
        }
    }

    std::vector<NDIndex<DIM> > localEmbeddedHDomains;
    std::vector<NDIndex<DIM> > localSurroundingHDomains;
    Utils::getLocalDomains(embeddedH.getLayout(), localEmbeddedHDomains);
    Utils::getLocalDomains(surroundingH.getLayout(), localSurroundingHDomains);
    for (unsigned int k = 0; k < numNodes; ++ k) {
        NDIndex<DIM> & dom = localEmbeddedHDomains[k];
        for (unsigned int d = 0; d < DIM; ++ d) {
            if (dom[d].first() == gEmbeddedHDomain[d].first())
                dom[d] = Index(dom[d].first() - 1, dom[d].last());
            if (dom[d].last() == gEmbeddedHDomain[d].last())
                dom[d] = Index(dom[d].first(), dom[d].last() + 1);
        }
    }

    if (move[0] == 1) {
        NDIndex<DIM> dom = gEmbeddedHDomain;
        dom[0] = Index(dom[0].first() - move[0], dom[0].first() - move[0]);
        for (unsigned int d = 1; d < DIM; ++ d) {
            dom[d] = Index(dom[d].first() - move[d], dom[d].last() - move[d] - 1);
        }
        HzE2SDoms.push_back(dom);

        dom = gEmbeddedEDomain;
        dom[0] = Index(dom[0].first() - move[0], dom[0].first() - move[0]);
        for (unsigned int d = 1; d < DIM; ++ d) {
            dom[d] = Index(dom[d].first() - move[d] + 1, dom[d].last() - move[d]);
        }
        ExE2SDoms.push_back(dom);

        dom = gEmbeddedEDomain;
        dom[0] = Index(dom[0].first() + 1 - move[0], dom[0].first() + 1 - move[0]);
        for (unsigned int d = 1; d < DIM; ++ d) {
            dom[d] = Index(dom[d].first() - move[d], dom[d].last() - move[d]);
        }
        EyE2SDoms.push_back(dom);

        dom = gEmbeddedHDomain;
        dom[0] = Index(dom[0].last() + positionEmbedded[0],
                       dom[0].last() + positionEmbedded[0]);
        for (unsigned int d = 1; d < DIM; ++ d) {
            dom[d] = Index(dom[d].first() + positionEmbedded[d],
                           dom[d].last() + positionEmbedded[d]);
        }
        HzS2EDoms.push_back(dom);

        dom = gEmbeddedEDomain;
        dom[0] = Index(dom[0].last() + positionEmbedded[0],
                       dom[0].last() + positionEmbedded[0]);
        for (unsigned int d = 1; d < DIM; ++ d) {
            dom[d] = Index(dom[d].first() + positionEmbedded[d],
                           dom[d].last() + positionEmbedded[d] + 1);
        }
        ExS2EDoms.push_back(dom);

        dom = gEmbeddedEDomain;
        dom[0] = Index(dom[0].last() + positionEmbedded[0] + 1,
                       dom[0].last() + positionEmbedded[0] + 1);
        for (unsigned int d = 1; d < DIM; ++ d) {
            dom[d] = Index(dom[d].first() + positionEmbedded[d],
                           dom[d].last() + positionEmbedded[d]);
        }
        EyS2EDoms.push_back(dom);
    }

    if (move[1] == 1) {
        NDIndex<DIM> dom = gEmbeddedHDomain;
        dom[1] = Index(dom[1].first() - move[1], dom[1].first() - move[1]);
        for (unsigned int d = 1; d < DIM; ++ d) {
            unsigned int d2 = (d + 1) % DIM;
            dom[d2] = Index(dom[d2].first() - move[d2], dom[d2].last() - 1 - move[d2]);
        }
        HzE2SDoms.push_back(dom);

        dom = gEmbeddedEDomain;
        dom[1] = Index(dom[1].first() + 1 - move[1], dom[1].first() + 1 - move[1]);
        for (unsigned int d = 1; d < DIM; ++ d) {
            unsigned int d2 = (d + 1) % DIM;
            dom[d2] = Index(dom[d2].first() - move[d2], dom[d2].last() - move[d2]);
        }
        ExE2SDoms.push_back(dom);

        dom = gEmbeddedEDomain;
        dom[1] = Index(dom[1].first() - move[1], dom[1].first() - move[1]);
        for (unsigned int d = 1; d < DIM; ++ d) {
            unsigned int d2 = (d + 1) % DIM;
            dom[d2] = Index(dom[d2].first() + 1 - move[d2], dom[d2].last() - move[d2]);
        }
        EyE2SDoms.push_back(dom);

        dom = gEmbeddedHDomain;
        dom[1] = Index(dom[1].last() + positionEmbedded[1],
                       dom[1].last() + positionEmbedded[1]);
        for (unsigned int d = 1; d < DIM; ++ d) {
            unsigned int d2 = (d + 1) % DIM;
            dom[d2] = Index(dom[d2].first() + positionEmbedded[d2],
                            dom[d2].last() + positionEmbedded[d2]);
        }
        HzS2EDoms.push_back(dom);

        dom = gEmbeddedEDomain;
        dom[1] = Index(dom[1].last() + 1 + positionEmbedded[1],
                       dom[1].last() + 1 + positionEmbedded[1]);
        for (unsigned int d = 1; d < DIM; ++ d) {
            unsigned int d2 = (d + 1) % DIM;
            dom[d2] = Index(dom[d2].first() + positionEmbedded[d2],
                            dom[d2].last() + positionEmbedded[d2]);
        }
        ExS2EDoms.push_back(dom);

        dom = gEmbeddedEDomain;
        dom[1] = Index(dom[1].last() + positionEmbedded[1],
                       dom[1].last() + positionEmbedded[1]);
        for (unsigned int d = 1; d < DIM; ++ d) {
            unsigned int d2 = (d + 1) % DIM;
            dom[d2] = Index(dom[d2].first() + positionEmbedded[d2],
                            dom[d2].last() + positionEmbedded[d2] + 1);
        }
        EyS2EDoms.push_back(dom);

    } else if (move[1] == -1) {
        NDIndex<DIM> dom = gEmbeddedHDomain;
        dom[1] = Index(dom[1].last() - 1 - move[1], dom[1].last() - move[1] - 1);
        for (unsigned int d = 1; d < DIM; ++ d) {
            unsigned int d2 = (d + 1) % DIM;
            dom[d2] = Index(dom[d2].first() - move[d2], dom[d2].last() - move[d2] - 1);
        }
        HzE2SDoms.push_back(dom);

        dom = gEmbeddedEDomain;
        dom[1] = Index(dom[1].last() - move[1], dom[1].last() - move[1]);
        for (unsigned int d = 1; d < DIM; ++ d) {
            unsigned int d2 = (d + 1) % DIM;
            dom[d2] = Index(dom[d2].first() - move[d2], dom[d2].last() - move[d2]);
        }
        ExE2SDoms.push_back(dom);

        dom = gEmbeddedEDomain;
        dom[1] = Index(dom[1].last() - move[1], dom[1].last() - move[1]);
        for (unsigned int d = 1; d < DIM; ++ d) {
            unsigned int d2 = (d + 1) % DIM;
            dom[d2] = Index(dom[d2].first() - move[d2] + 1, dom[d2].last() - move[d2]);
        }
        EyE2SDoms.push_back(dom);

        dom = gEmbeddedHDomain;
        dom[1] = Index(dom[1].first() - 1 + positionEmbedded[1],
                       dom[1].first() - 1 + positionEmbedded[1]);
        for (unsigned int d = 1; d < DIM; ++ d) {
            unsigned int d2 = (d + 1) % DIM;
            dom[d2] = Index(dom[d2].first() + positionEmbedded[d2],
                            dom[d2].last() + positionEmbedded[d2]);
        }
        HzS2EDoms.push_back(dom);

        dom = gEmbeddedEDomain;
        dom[1] = Index(dom[1].first() + positionEmbedded[1],
                       dom[1].first() + positionEmbedded[1]);
        for (unsigned int d = 1; d < DIM; ++ d) {
            unsigned int d2 = (d + 1) % DIM;
            dom[d2] = Index(dom[d2].first() + positionEmbedded[d2],
                            dom[d2].last() + positionEmbedded[d2]);
        }
        ExS2EDoms.push_back(dom);

        dom = gEmbeddedEDomain;
        dom[1] = Index(dom[1].first() + positionEmbedded[1],
                       dom[1].first() + positionEmbedded[1]);
        for (unsigned int d = 1; d < DIM; ++ d) {
            unsigned int d2 = (d + 1) % DIM;
            dom[d2] = Index(dom[d2].first() + positionEmbedded[d2],
                            dom[d2].last() + positionEmbedded[d2] + 1);
        }
        EyS2EDoms.push_back(dom);
    }

    const std::vector<std::vector<NDIndex<DIM> > > allLayers = {HzE2SDoms, ExE2SDoms, EyE2SDoms, HzS2EDoms, ExS2EDoms, EyS2EDoms};
    const std::vector<std::vector<NDIndex<DIM> > > allLocalDomains = {localEmbeddedHDomains,
                                                                      localEmbeddedEDomains,
                                                                      localSurroundingHDomains,
                                                                      localSurroundingEDomains};
    Vektor<int,DIM>  diffs[] = {-positionEmbedded, positionEmbedded};

    std::vector<MPI_Request> requests;
    int tag = Ippl::Comm->next_tag(F_GUARD_CELLS_TAG, F_TAG_CYCLE);

    std::vector<std::vector<char> > sendData(numNodes);
    std::vector<std::vector<char> > receiveData(numNodes);
    // determine the local contribution to the boundary layers and fill in the values
    for (unsigned int k = 0; k < numNodes; ++ k) {
        if (k == myNode) continue;

        auto & data = receiveData[k];

        for (unsigned int layerType = 0; layerType < 6; ++ layerType) {
            const std::vector<NDIndex<DIM> > & currentLayers = allLayers[layerType];
            unsigned int didx;
            {
                unsigned int mul = layerType / 3;
                unsigned int add = (layerType % 3 == 0? 0: 1);
                didx = 2 * mul + add;
            }
            const NDIndex<DIM> & domain = allLocalDomains[didx][k];

            for (unsigned int lidx = 0; lidx < currentLayers.size(); ++ lidx) {
                const NDIndex<DIM> & layer = currentLayers[lidx];

                NDIndex<DIM> dom;
                bool hasSize = true;
                for (unsigned int d = 0; d < DIM; ++ d) {
                    int lower = std::max(domain[d].first(), layer[d].first());
                    int upper = std::min(domain[d].last(), layer[d].last());
                    if (lower > upper) {
                        hasSize = false;
                        break;
                    }
                    dom[d] = Index(lower, upper);
                }
                if (!hasSize) continue;

                unsigned int idx = (layerType < 3? 0: 1);

                bool inside = true;
                const NDIndex<DIM> & flddom = allLocalDomains[(didx + 2) % 4][myNode];
                NDIndex<DIM> overlap, adjdom;
                for (unsigned int d2 = 0; d2 < DIM; ++ d2) {
                    adjdom[d2] = Index(flddom[d2].first() + diffs[idx][d2],
                                       flddom[d2].last() + diffs[idx][d2]);
                }
                // DBGOUT << "k = " << k << "\tlayerType = " << layerType << "\tlidx = " << lidx << "\t didx compl = " << (didx + 2) % 4 << std::endl;
                // DBGOUT << "dom = " << dom << "\tflddom = " << flddom << "\tdiff = " << diffs[idx] << std::endl;
                for (unsigned int d2 = 0; d2 < DIM; ++ d2) {
                    int lower = std::max(dom[d2].first(), flddom[d2].first() + diffs[idx][d2]);
                    int upper = std::min(dom[d2].last(),  flddom[d2].last() + diffs[idx][d2]);
                    if (lower > upper) {
                        inside = false;
                        break;
                    }
                    overlap[d2] = Index(lower, upper);
                    adjdom[d2] = Index(lower - diffs[idx][d2], upper - diffs[idx][d2]);
                }
                if (!inside) continue;

                unsigned int sizeData = overlap.size();
                data.resize(data.size() + 1 + 2 * DIM * sizeof(int) + MD5_DIGEST_LENGTH + sizeData * sizeof(double));

            }
        }
        if (receiveData[k].size() == 0) continue;

        MPI_Request req;
        MPI_Irecv(&(receiveData[k][0]), receiveData[k].size(), MPI_BYTE, k, tag, MPI_COMM_WORLD, &req);
        requests.push_back(req);
    }

    // determine the local contribution to the boundary layers and fill in the values
    for (unsigned int k = 0; k < numNodes; ++ k) {
        unsigned int start = sendData[k].size();
        auto & data = (k != myNode? sendData[k]: receiveData[k]);

        for (unsigned int layerType = 0; layerType < 6; ++ layerType) {
            const std::vector<NDIndex<DIM> > & currentLayers = allLayers[layerType];
            unsigned int didx;
            {
                unsigned int mul = layerType / 3;
                unsigned int add = (layerType % 3 == 0? 0: 1);
                didx = 2 * mul + add;
            }
            const NDIndex<DIM> & domain = allLocalDomains[didx][myNode];

            for (unsigned int lidx = 0; lidx < currentLayers.size(); ++ lidx) {
                const NDIndex<DIM> & layer = currentLayers[lidx];

                NDIndex<DIM> dom;
                bool hasSize = true;
                for (unsigned int d = 0; d < DIM; ++ d) {
                    int lower = std::max(domain[d].first(), layer[d].first());
                    int upper = std::min(domain[d].last(), layer[d].last());
                    if (lower > upper) {
                        hasSize = false;
                        break;
                    }
                    dom[d] = Index(lower, upper);
                }
                if (!hasSize) continue;

                unsigned int idx = (layerType < 3? 0: 1);

                bool inside = true;
                const NDIndex<DIM> & flddom = allLocalDomains[(didx + 2) % 4][k];
                NDIndex<DIM> overlap, adjdom;
                for (unsigned int d2 = 0; d2 < DIM; ++ d2) {
                    adjdom[d2] = Index(flddom[d2].first() + diffs[idx][d2],
                                       flddom[d2].last() + diffs[idx][d2]);
                }

                for (unsigned int d2 = 0; d2 < DIM; ++ d2) {
                    int lower = std::max(dom[d2].first(), flddom[d2].first() + diffs[idx][d2]);
                    int upper = std::min(dom[d2].last(),  flddom[d2].last() + diffs[idx][d2]);
                    if (lower > upper) {
                        inside = false;
                        break;
                    }
                    overlap[d2] = Index(lower, upper);
                    adjdom[d2] = Index(lower - diffs[idx][d2], upper - diffs[idx][d2]);
                }
                if (!inside) continue;

                unsigned int sizeData = overlap.size();
                data.reserve(data.size() + 1 + 2 * DIM * sizeof(int) + MD5_DIGEST_LENGTH + sizeData * sizeof(double));
                const char* buffer;

                char FType = layerType;
                data.insert(data.end(), &FType, &FType + 1);

                for (unsigned int d2 = 0; d2 < DIM; ++ d2) {
                    int lower = overlap[d2].first() - diffs[idx][d2];
                    buffer = reinterpret_cast<const char*>(&lower);
                    data.insert(data.end(), buffer, buffer + sizeof(int));

                    int upper = overlap[d2].last() - diffs[idx][d2];
                    buffer = reinterpret_cast<const char*>(&upper);
                    data.insert(data.end(), buffer, buffer + sizeof(int));
                }

                double * rawData = new double[sizeData];

                switch(layerType) {
                case 0:
                    {
                        Field<Vector_t, DIM, Mesh_t, Cell> & fld = embeddedH;
                        serializeField(rawData, fld, overlap, HFieldToDouble());
                        break;
                    }
                case 1:
                    {
                        Field<Vector_t, DIM, Mesh_t, Edge> & fld = embeddedE;
                        serializeField(rawData, fld, overlap, ExFieldToDouble());
                        break;
                    }
                case 2:
                    {
                        Field<Vector_t, DIM, Mesh_t, Edge> & fld = embeddedE;
                        serializeField(rawData, fld, overlap, EyFieldToDouble());
                        break;
                    }
                case 3:
                    {
                        Field<Vector_t, DIM, Mesh_t, Cell> & fld = surroundingH;
                        serializeField(rawData, fld, overlap, HFieldToDouble());
                        break;
                    }
                case 4:
                    {
                        Field<Vector_t, DIM, Mesh_t, Edge> & fld = surroundingE;
                        serializeField(rawData, fld, overlap, ExFieldToDouble());
                        break;
                    }
                case 5:
                    {
                        Field<Vector_t, DIM, Mesh_t, Edge> & fld = surroundingE;
                        serializeField(rawData, fld, overlap, EyFieldToDouble());
                        break;
                    }

                }

                buffer = reinterpret_cast<const char*>(rawData);
#ifdef OpenSSL_FOUND
                unsigned char md5Hash[MD5_DIGEST_LENGTH];
                MD5(reinterpret_cast<const unsigned char*>(buffer), sizeData * sizeof(double), md5Hash);
                data.insert(data.end(), md5Hash, md5Hash + MD5_DIGEST_LENGTH);
#endif
                data.insert(data.end(), buffer, buffer + sizeData * sizeof(double));

                delete[] rawData;

            }
        }

        if (sendData[k].size() - start == 0) continue;

        MPI_Send(&(sendData[k][0]), sendData[k].size(), MPI_BYTE, k, tag, MPI_COMM_WORLD);
    }


    std::vector<MPI_Status> status(requests.size());
    MPI_Waitall(requests.size(), &requests[0], &status[0]);

    for (unsigned int k = 0; k < numNodes; ++ k) {
        std::vector<char> & data = receiveData[k];

        unsigned int size = data.size();
        if (size == 0) continue;

        NDIndex<DIM> overlap;
        unsigned int it = 0;
        while (it < size) {
            unsigned int FType = data[it ++];

            for (unsigned int d2 = 0; d2 < DIM; ++ d2) {
                int lower = *reinterpret_cast<const int*>(&data[it]);
                it += sizeof(int);
                int upper = *reinterpret_cast<const int*>(&data[it]);
                it += sizeof(int);
                overlap[d2] = Index(lower, upper);
            }
            unsigned int sizeData = overlap.size();

#ifdef OpenSSL_FOUND
            unsigned char sentMd5Hash[MD5_DIGEST_LENGTH], receivedMd5Hash[MD5_DIGEST_LENGTH];
            char rcMD5Buffer[2 * MD5_DIGEST_LENGTH], snMD5Buffer[2 * MD5_DIGEST_LENGTH];
            memcpy(sentMd5Hash, &data[it], MD5_DIGEST_LENGTH);
            it += MD5_DIGEST_LENGTH;

            MD5(reinterpret_cast<const unsigned char*>(&data[it]), sizeData * sizeof(double), receivedMd5Hash);
            unsigned int ll = 0;
            for (unsigned int l = 0; l < MD5_DIGEST_LENGTH; ++ l) {
                sprintf(&rcMD5Buffer[ll], "%02x", receivedMd5Hash[l]);
                sprintf(&snMD5Buffer[ll], "%02x", sentMd5Hash[l]);
                ll += 2;
            }
            if (strcmp(rcMD5Buffer, snMD5Buffer) != 0) {
                dbg << "Communicator.hh: " << __LINE__ << "\t"
                    << "MD5 hashes are not equal" << "\n"
                    << "    sentMd5Hash = " << snMD5Buffer << "\n"
                    << "recievedMd5Hash = " << rcMD5Buffer << std::endl;
            }
#endif

            double * rawData = new double[sizeData];
            for (unsigned int l = 0; l < sizeData; ++ l) {
                rawData[l] = *reinterpret_cast<const double*>(&data[it]);
                it += sizeof(double);
            }

            switch (FType) {
            case 0:
                {
                    Field<Vector_t, DIM, Mesh_t, Cell> & fld = surroundingH;
                    deserializeField(fld, rawData, overlap, HFieldToDouble());
                    break;
                }
            case 1:
                {
                    Field<Vector_t, DIM, Mesh_t, Edge> & fld = surroundingE;
                    deserializeField(fld, rawData, overlap, ExFieldToDouble());
                    break;
                }
            case 2:
                {
                    Field<Vector_t, DIM, Mesh_t, Edge> & fld = surroundingE;
                    deserializeField(fld, rawData, overlap, EyFieldToDouble());
                    break;
                }
            case 3:
                {
                    Field<Vector_t, DIM, Mesh_t, Cell> & fld = embeddedH;
                    deserializeField(fld, rawData, overlap, HFieldToDouble());
                    break;
                }
            case 4:
                {
                    Field<Vector_t, DIM, Mesh_t, Edge> & fld = embeddedE;
                    deserializeField(fld, rawData, overlap, ExFieldToDouble());
                    break;
                }
            case 5:
                {
                    Field<Vector_t, DIM, Mesh_t, Edge> & fld = embeddedE;
                    deserializeField(fld, rawData, overlap, EyFieldToDouble());
                    break;
                }
            }
            delete[] rawData;
        }
    }

    embeddedE.fillGuardCells();
    embeddedH.fillGuardCells();
    surroundingE.fillGuardCells();
    surroundingH.fillGuardCells();
}

const unsigned int HFieldToDouble::sides = 1 << (2 * DIM)| 1 << (1 + 2 * DIM) | 1 << 2 | 1 << 3;
const unsigned int ExFieldToDouble::sides = 1 << (3 + 2 * DIM) | 1 << 1;
const unsigned int EyFieldToDouble::sides = 1 << (2 + 2 * DIM) | 1;
