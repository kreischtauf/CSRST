/***************************************************************************
                           Communicator.hh
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

#ifndef COMMUNICATOR_HH
#define COMMUNICATOR_HH

#include "defs.hh"
#include "utils.hh"
#include "FieldPatch.hh"
// #include "BoundaryLayer.hh"

#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/split_member.hpp>

#ifdef OpenSSL_FOUND
#include <openssl/md5.h>
#else
#define MD5_DIGEST_LENGTH 0
#endif

#include <numeric>

#define P2PCOMM

struct ExFieldToDouble;
struct EyFieldToDouble;
struct HFieldToDouble;

class Communicator {
public:
    template<class T, class C>
    static void copy(const Field<T, DIM, Mesh_t, C> & from,
                     Field<T, DIM, Mesh_t, C> & to,
                     const Vektor<int,DIM> & start,
                     const bool all2all = false);

    template<class T, class C, class Op>
    static void exchangeBoundariesCC(Field<T, DIM, Mesh_t, C> & surrounding,
                                     Field<T, DIM, Mesh_t, C> & embedded,
                                     const Vektor<int,DIM> & positionEmbedded,
                                     const Op & op,
                                     const bool forwardDir = false);

    template<class T, class C, class Op>
    static void exchangeBoundaries(Field<T, DIM, Mesh_t, C> & surrounding,
                                   Field<T, DIM, Mesh_t, C> & embedded,
                                   const Vektor<int,DIM> & positionEmbedded,
                                   const Op & op);

    static void exchangeBoundariesMovingMesh(VField_Edge_t & surroundingE,
                                             VField_Cell_t & surroundingH,
                                             VField_Edge_t & embeddedE,
                                             VField_Cell_t & embeddedH,
                                             const Vektor<int, DIM> & positionEmbedded,
                                             const Vektor<int, DIM> & move);

    template<class T, class C, class Op>
    static void serializeField(double * rawData,
                               const Field<T, DIM, Mesh_t, C> & fld,
                               const NDIndex<DIM> & dom,
                               const Op & op);

    template<class T, class C, class Op>
    static void deserializeField(Field<T, DIM, Mesh_t, C> & fld,
                                 const double * rawData,
                                 const NDIndex<DIM> & dom,
                                 const Op & op);

    static std::vector<NDIndex<DIM> > getAllLocalDomains(NDIndex<DIM> lDom);

    template <class T>
    static void communicateFields(FieldPatch<T> & fp,
                                  const NDIndex<DIM> & dom,
                                  const Timings::TimerRef & commTimer,
                                  const bool  all2all = false);

    template <class T>
    static void communicateFields(FieldPatch<T> & fp,
                                  const std::vector<NDIndex<DIM> >  & localPDomains,
                                  const std::vector<NDIndex<DIM> >  & localFDomains,
                                  const Timings::TimerRef & commTimer,
                                  const bool all2all = false);

    template <class T>
    static void collectivelyCommunicateFields(FieldPatch<T> & fp,
                                              const NDIndex<DIM> & dom,
                                              const Timings::TimerRef & commTimer);

    template <class T>
    static void p2pCommunicateFields(FieldPatch<T> & fp,
                                     const std::vector<NDIndex<DIM> >  & localPDomains,
                                     const std::vector<NDIndex<DIM> >  & localFDomains,
                                     const Timings::TimerRef & commTimer,
                                     const bool all2all = false);

    static void initialize();
    static const std::vector<unsigned int> & getNeighbours();
    static const std::vector<unsigned int> & getNode2Index();
    static const std::vector<unsigned int> & getIndex2Node();
    static std::pair<unsigned int, unsigned int> getNodeMeshSize();

private:
    static std::vector<unsigned int> _index2node;
    static std::vector<unsigned int> _node2index;
    static std::vector<unsigned int> _neighbours;
    static unsigned int _nNodeCols;
    static unsigned int _nNodeRows;
    static bool _initialized;
};

class WrappedVector: public Vector_t {
public:
    WrappedVector():
        Vector_t()
    { }

    WrappedVector(const Vector_t &rhs):
        Vector_t(rhs)
    { }

    WrappedVector(const double& x00):
        Vector_t(x00)
    { }

    WrappedVector(const double& x00, const double& x01):
        Vector_t(x00, x01)
    { }

    WrappedVector(const double& x00, const double& x01, const double& x02):
        Vector_t(x00, x01, x02)
    { }

private:
    friend class boost::serialization::access;

    template<class Archive>
    void save(Archive& ar, const unsigned int version) const
    {
        for (unsigned int d = 0; d < DIM; ++ d) {
            double a = (*this)[d];
            ar & a;
        }
    }

    template<class Archive>
    void load(Archive& ar, const unsigned int version)
    {
        for (unsigned int d = 0; d < DIM; ++ d) {
            double a;
            ar & a;
            (*this)[d] = a;
        }
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()

};

template <class T>
class Wrapper
{
public:
    typedef T W_t;
};

template<>
struct Wrapper<Vector_t>
{
    typedef WrappedVector W_t;
};

template<class T>
class BoundaryLayer: public FieldPatch<T> {
public:
    enum BLType {EMBEDDED, SURROUNDING};
    BoundaryLayer();
    BoundaryLayer(BLType type, NDIndex<DIM> dom);
    BoundaryLayer(const BoundaryLayer<T>& bl);

    const BLType & getType() const { return _bltype;}

private:
    friend class boost::serialization::access;

    template<class Archive>
    void save(Archive& ar, const unsigned int version) const;

    template<class Archive>
    void load(Archive& ar, const unsigned int version);

    BOOST_SERIALIZATION_SPLIT_MEMBER()

    BLType _bltype;
};

template<class T, class C>
void Communicator::copy(const Field<T, DIM, Mesh_t, C> & from,
                        Field<T, DIM, Mesh_t, C> & to,
                        const Vektor<int,DIM> & start,
                        const bool all2all)
{
    typedef T T_t;
    typedef typename Wrapper<T_t>::W_t W_t;

    NDIndex<DIM> lFromDomain = from.getLayout().getLocalNDIndex();
    NDIndex<DIM> gFromDomain = from.getLayout().getDomain();
    NDIndex<DIM> lToDomain = to.getLayout().getLocalNDIndex();
    NDIndex<DIM> gToDomain = to.getLayout().getDomain();
    int lowerI = std::max(lFromDomain[0].first(), gToDomain[0].first() + start[0]);
    int upperI = std::min(lFromDomain[0].last(),  gToDomain[0].last() + start[0]);
    int signI = lowerI <= upperI ? 1 : upperI - lowerI;

    int lowerJ = std::max(lFromDomain[1].first(), gToDomain[1].first() + start[1]);
    int upperJ = std::min(lFromDomain[1].last(),  gToDomain[1].last() + start[1]);
    int signJ = lowerJ <= upperJ ? 1 : upperJ - lowerJ;

    NDIndex<DIM> overlap(Index(lowerI, upperI, signI),
                         Index(lowerJ, upperJ, signJ));
    FieldPatch<W_t> helper(overlap);

    NDIndex<DIM> elem;
    for (int j = lowerJ; j <= upperJ; ++ j) {
        elem[1] = Index(j,j);
        for (int i = lowerI; i <= upperI; ++ i) {
            elem[0] = Index(i,i);
            helper(i,j) = from.localElement(elem);
        }
    }

    lowerI = std::max(gFromDomain[0].first(), lToDomain[0].first() + start[0]);
    upperI = std::min(gFromDomain[0].last(), lToDomain[0].last() + start[0]);
    signI = lowerI <= upperI ? 1: upperI - lowerI;

    lowerJ = std::max(gFromDomain[1].first(), lToDomain[1].first() + start[1]);
    upperJ = std::min(gFromDomain[1].last(), lToDomain[1].last() + start[1]);
    signJ = lowerJ <= upperJ ? 1 : upperJ - lowerJ;

    NDIndex<DIM> lHelperDomain(Index(lowerI, upperI, signI),
                               Index(lowerJ, upperJ, signJ));


    Timings::TimerRef comm_timer;
    comm_timer = Timings::getTimer("field exchange");
    communicateFields(helper, lHelperDomain, comm_timer, all2all);

    int jj = lowerJ - start[1];
    for (int j = lowerJ; j <= upperJ; ++ j, ++ jj) {
        elem[1] = Index(jj,jj);
        int ii = lowerI - start[0];
        for (int i = lowerI; i <= upperI; ++ i, ++ ii) {
            elem[0] = Index(ii,ii);
            to.localElement(elem) = helper(i,j);
        }
    }

    to.fillGuardCells();
}


template<class T, class C, class Op>
void Communicator::exchangeBoundariesCC(Field<T, DIM, Mesh_t, C> & surrounding,
                                        Field<T, DIM, Mesh_t, C> & embedded,
                                        const Vektor<int,DIM> & positionEmbedded,
                                        const Op & op,
                                        const bool forwardDir)
{
    extern std::ofstream dbg;

    typedef T T_t;
    typedef typename Wrapper<T_t>::W_t W_t;

    // unsigned int myNode = Ippl::myNode();
    // unsigned int numNodes = Ippl::getNodes();

    NDIndex<DIM> elem;
    NDIndex<DIM> gSurroundingDomain = surrounding.getLayout().getDomain();
    NDIndex<DIM> gEmbeddedDomain = embedded.getLayout().getDomain();
    NDIndex<DIM> incrGEmbeddedDomain; // one cell thick layer around and outside of the embedded domain
    for (unsigned int d = 0; d < DIM; ++ d) {
        incrGEmbeddedDomain[d] = Index(gEmbeddedDomain[d].first() - 1 + positionEmbedded[d],
                                       gEmbeddedDomain[d].last() + 1 + positionEmbedded[d]);
    }
    std::vector<Vektor<int, DIM> > startEmbeddedDomain(2 * DIM);
    std::vector<Vektor<int, DIM> > startSurroundingDomain(2 * DIM);
    std::vector<int> startInnerLayer(2 * DIM, 0);
    std::vector<int> startOuterLayer(2 * DIM, 0);

    unsigned int sizeInnerLayer = 0, sizeOuterLayer = 0;
    std::vector<int> LayerSizes(2 * DIM);
    for (unsigned int d = 0; d < DIM; ++ d) {
        NDIndex<DIM> dom = gEmbeddedDomain;
        for (unsigned int d2 = 0; d2 < d; ++ d2) {
            dom[d2] = Index(dom[d2].first() + 1, dom[d2].last() - 1);
        }
        dom[d] = Index(dom[d].first(), dom[d].first());
        sizeInnerLayer += 2 * dom.size();
        LayerSizes[d] = dom.size();
        LayerSizes[d + DIM] = dom.size();

        dom = incrGEmbeddedDomain;
        for (unsigned int d2 = 0; d2 < d; ++ d2) {
            dom[d2] = Index(dom[d2].first() + 1, dom[d2].last() - 1);
        }
        dom[d] = Index(dom[d].first(), dom[d].first());
        sizeOuterLayer += 2 * dom.size();
    }

    for (unsigned int d = 1; d < 2 * DIM; ++ d) {
        startInnerLayer[d] = LayerSizes[d-1] + startInnerLayer[d-1];
        startOuterLayer[d] = LayerSizes[d-1] + startOuterLayer[d-1] + 2;
    }

    double * innerLayer = new double[sizeInnerLayer + sizeOuterLayer];
    double * outerLayer = innerLayer + sizeInnerLayer;
    std::fill(innerLayer, innerLayer + sizeInnerLayer + sizeOuterLayer, 0.0);

    NDIndex<DIM> lEDom = embedded.getLayout().getLocalNDIndex();
    NDIndex<DIM> lSDom = surrounding.getLayout().getLocalNDIndex();

    // the boundary layer is treated as 2 * DIM one-dimensional domains
    // compute the global domains of these one-dimensional fields
    NDIndex<DIM> gDoms[4 * DIM]; // global boundary layers
    NDIndex<DIM> exDoms[4]; // local contribution to global boundary layer
    for (unsigned int d = 0; d < DIM; ++ d) {
        gDoms[d] = gEmbeddedDomain;
        if (forwardDir && d == 0) { // in the direction of propagation a layer one cell further inside has to be communicated
            // due to the shift of the embedded domain after each step
            gDoms[d][0] = Index(gDoms[0][0].first() + 1, gDoms[0][0].first() + 1);
        } else {
            gDoms[d][d] = Index(gDoms[d][d].first(), gDoms[d][d].first());
        }

        gDoms[d + 2] = gEmbeddedDomain;
        gDoms[d + 2][d] = Index(gDoms[d + 2][d].last(), gDoms[d + 2][d].last());

        gDoms[d + 4] = incrGEmbeddedDomain;
        gDoms[d + 4][d] = Index(gDoms[d + 4][d].first(), gDoms[d + 4][d].first());

        gDoms[d + 6] = incrGEmbeddedDomain;
        gDoms[d + 6][d] = Index(gDoms[d + 6][d].last(), gDoms[d + 6][d].last());

        for (unsigned int d2 = 0; d2 < d; ++ d2) { // make sure that there is no overlap within gS2EDoms and within gE2SDoms
            gDoms[d    ][d2] = Index(gDoms[d    ][d2].first() + 1, gDoms[d    ][d2].last() - 1);
            gDoms[d + 2][d2] = Index(gDoms[d + 2][d2].first() + 1, gDoms[d + 2][d2].last() - 1);
            gDoms[d + 4][d2] = Index(gDoms[d + 4][d2].first() + 1, gDoms[d + 4][d2].last() - 1);
            gDoms[d + 6][d2] = Index(gDoms[d + 6][d2].first() + 1, gDoms[d + 6][d2].last() - 1);
        }

        for (unsigned int d2 = 0; d2 < DIM; ++ d2) {
            startEmbeddedDomain[d](d2)          = gDoms[d][d2].min();
            startEmbeddedDomain[d + DIM](d2)    = gDoms[d + DIM][d2].min();
            startSurroundingDomain[d](d2)       = gDoms[d + 2 * DIM][d2].min();
            startSurroundingDomain[d + DIM](d2) = gDoms[d + 3 * DIM][d2].min();
        }

        startEmbeddedDomain[d + DIM](d) = gDoms[d + DIM][d].min();
        startSurroundingDomain[d + DIM](d) = gDoms[d + 3 * DIM][d].min();
    }

    // determine the local contribution to the boundary layers and fill in the values
    for (unsigned int d = 0; d < DIM; ++ d) {
        // calculate overlap between the one-dimensional domains and the local DIM-dimensional domains
        for (unsigned int d2 = 0; d2 < DIM; ++ d2) {
            int lower = std::max(lEDom[d2].first(), gDoms[d][d2].first());
            int upper = std::min(lEDom[d2].last(),  gDoms[d][d2].last());
            int sign = lower <= upper? 1: upper - lower;
            exDoms[0][d2] = Index(lower, upper, sign);

            lower = std::max(lEDom[d2].first(), gDoms[d + 2][d2].first());
            upper = std::min(lEDom[d2].last(),  gDoms[d + 2][d2].last());
            sign = lower <= upper? 1: upper - lower;
            exDoms[1][d2] = Index(lower, upper, sign);

            lower = std::max(lSDom[d2].first(), gDoms[d + 4][d2].first());
            upper = std::min(lSDom[d2].last(),  gDoms[d + 4][d2].last());
            sign = lower <= upper? 1: upper - lower;
            exDoms[2][d2] = Index(lower, upper, sign);

            lower = std::max(lSDom[d2].first(), gDoms[d + 6][d2].first());
            upper = std::min(lSDom[d2].last(),  gDoms[d + 6][d2].last());
            sign = lower <= upper? 1: upper - lower;
            exDoms[3][d2] = Index(lower, upper, sign);
        }

        unsigned int side = 0;
        for (const auto & bld: {std::make_pair(BoundaryLayer<W_t>::EMBEDDED, exDoms[0]),
                    std::make_pair(BoundaryLayer<W_t>::EMBEDDED, exDoms[1]),
                    std::make_pair(BoundaryLayer<W_t>::SURROUNDING, exDoms[2]),
                    std::make_pair(BoundaryLayer<W_t>::SURROUNDING, exDoms[3])}) {
            ++ side;
            const NDIndex<DIM> & dom = bld.second;

            bool hasSize = true;
            for (unsigned int d2 = 0; d2 < DIM; ++ d2) {
                if (dom[d2].first() > dom[d2].last()) {
                    hasSize = false;
                    break;
                }
            }
            if (!hasSize) continue;

            double *bl = NULL;
            unsigned int startLayer = 0;
            Vektor<int, DIM> startDomain;
            Field<T, DIM, Mesh_t, C> * fld = NULL;
            switch(side) {
            case 1:
                bl = innerLayer;
                startLayer = startInnerLayer[d];
                startDomain = startEmbeddedDomain[d];
                fld = &embedded;
                break;
            case 2:
                bl = innerLayer;
                startLayer = startInnerLayer[d + DIM];
                startDomain = startEmbeddedDomain[d + DIM];
                fld = &embedded;
                break;
            case 3:
                bl = outerLayer;
                startLayer = startOuterLayer[d];
                startDomain = startSurroundingDomain[d];
                fld = &surrounding;
                break;
            case 4:
                bl = outerLayer;
                startLayer = startOuterLayer[d + DIM];
                startDomain = startSurroundingDomain[d + DIM];
                fld = &surrounding;
                break;
            }

            for (int j = dom[1].first(); j <= dom[1].last(); ++ j) {
                elem[1] = Index(j,j);
                int jj = j - startDomain[1];
                for (int i = dom[0].first(); i <= dom[0].last(); ++ i) {
                    elem[0] = Index(i,i);
                    int ii = i - startDomain[0] + jj + startLayer;
                    op.execute(*fld, elem, bl[ii]);
                }
            }
        }
    }

    dbg << "Communicator.hh: " << __LINE__  << "using collective communication" << std::endl;
    MPI_Allreduce(MPI_IN_PLACE, innerLayer, sizeInnerLayer + sizeOuterLayer, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (unsigned int d = 0; d < DIM; ++ d) {
        for (unsigned int d2 = 0; d2 < DIM; ++ d2) {
            gDoms[d          ][d2] = Index(gDoms[d          ][d2].first() + positionEmbedded[d2],
                                           gDoms[d          ][d2].last()  + positionEmbedded[d2]);
            gDoms[d +     DIM][d2] = Index(gDoms[d +     DIM][d2].first() + positionEmbedded[d2],
                                           gDoms[d +     DIM][d2].last()  + positionEmbedded[d2]);
            gDoms[d + 2 * DIM][d2] = Index(gDoms[d + 2 * DIM][d2].first() - positionEmbedded[d2],
                                           gDoms[d + 2 * DIM][d2].last()  - positionEmbedded[d2]);
            gDoms[d + 3 * DIM][d2] = Index(gDoms[d + 3 * DIM][d2].first() - positionEmbedded[d2],
                                           gDoms[d + 3 * DIM][d2].last()  - positionEmbedded[d2]);

            startEmbeddedDomain[d](d2)          = gDoms[d][d2].min();
            startEmbeddedDomain[d + DIM](d2)    = gDoms[d + DIM][d2].min();
            startSurroundingDomain[d](d2)       = gDoms[d + 2 * DIM][d2].min();
            startSurroundingDomain[d + DIM](d2) = gDoms[d + 3 * DIM][d2].min();
        }

        startEmbeddedDomain[d + DIM](d) = gDoms[d + DIM][d].min();
        startSurroundingDomain[d + DIM](d) = gDoms[d + 3 * DIM][d].min();
    }

    for (unsigned int d = 0; d < DIM; ++ d) {
        // calculate overlap between the one-dimensional domains and the local DIM-dimensional domains
        for (unsigned int d2 = 0; d2 < DIM; ++ d2) {
            int lower = std::max(lSDom[d2].first(), gDoms[d][d2].first());
            int upper = std::min(lSDom[d2].last(),  gDoms[d][d2].last());
            int sign = lower <= upper? 1: upper - lower;
            exDoms[0][d2] = Index(lower, upper, sign);

            lower = std::max(lSDom[d2].first(), gDoms[d + 2][d2].first());
            upper = std::min(lSDom[d2].last(),  gDoms[d + 2][d2].last());
            sign = lower <= upper? 1: upper - lower;
            exDoms[1][d2] = Index(lower, upper, sign);

            lower = std::max(lEDom[d2].first() - 1, gDoms[d + 4][d2].first());
            upper = std::min(lEDom[d2].last() + 1,  gDoms[d + 4][d2].last());
            sign = lower <= upper? 1: upper - lower;
            exDoms[2][d2] = Index(lower, upper, sign);

            lower = std::max(lEDom[d2].first() - 1, gDoms[d + 6][d2].first());
            upper = std::min(lEDom[d2].last() + 1,  gDoms[d + 6][d2].last());
            sign = lower <= upper? 1: upper - lower;
            exDoms[3][d2] = Index(lower, upper, sign);
        }

        unsigned int side = 0;
        for (const auto & bld: {std::make_pair(BoundaryLayer<W_t>::EMBEDDED, exDoms[0]),
                    std::make_pair(BoundaryLayer<W_t>::EMBEDDED, exDoms[1]),
                    std::make_pair(BoundaryLayer<W_t>::SURROUNDING, exDoms[2]),
                    std::make_pair(BoundaryLayer<W_t>::SURROUNDING, exDoms[3])}) {
            ++ side;
            const NDIndex<DIM> & dom = bld.second;

            bool hasSize = true;
            for (unsigned int d2 = 0; d2 < DIM; ++ d2) {
                if (dom[d2].first() > dom[d2].last()) {
                    hasSize = false;
                    break;
                }
            }
            if (!hasSize) continue;

            double *bl = NULL;
            unsigned int startLayer = 0;
            Vektor<int, DIM> startDomain;
            Field<T, DIM, Mesh_t, C> * fld = NULL;
            switch(side) {
            case 1:
                bl = innerLayer;
                startLayer = startInnerLayer[d];
                startDomain = startEmbeddedDomain[d];
                fld = &surrounding;
                break;
            case 2:
                bl = innerLayer;
                startLayer = startInnerLayer[d + DIM];
                startDomain = startEmbeddedDomain[d + DIM];
                fld = &surrounding;
                break;
            case 3:
                bl = outerLayer;
                startLayer = startOuterLayer[d];
                startDomain = startSurroundingDomain[d];
                fld = &embedded;
                break;
            case 4:
                bl = outerLayer;
                startLayer = startOuterLayer[d + DIM];
                startDomain = startSurroundingDomain[d + DIM];
                fld = &embedded;
                break;
            }

            for (int j = dom[1].first(); j <= dom[1].last(); ++ j) {
                elem[1] = Index(j,j);
                int jj = j - startDomain[1];
                for (int i = dom[0].first(); i <= dom[0].last(); ++ i) {
                    elem[0] = Index(i,i);
                    int ii = i - startDomain[0] + jj + startLayer;
                    op.invert(*fld, elem, bl[ii]);
                }
            }
        }
    }

    embedded.fillGuardCells();
    surrounding.fillGuardCells();

    delete[] innerLayer;
}


template<class T, class C, class Op>
void Communicator::exchangeBoundaries(Field<T, DIM, Mesh_t, C> & surrounding,
                                      Field<T, DIM, Mesh_t, C> & embedded,
                                      const Vektor<int,DIM> & positionEmbedded,
                                      const Op & op)
{
    typedef T T_t;
    typedef typename Wrapper<T_t>::W_t W_t;
    typedef C C_t;

    const unsigned int myNode = Ippl::myNode();
    const unsigned int numNodes = Ippl::getNodes();

    NDIndex<DIM> elem;
    NDIndex<DIM> gSurroundingDomain = surrounding.getLayout().getDomain();
    NDIndex<DIM> gEmbeddedDomain = embedded.getLayout().getDomain();
    NDIndex<DIM> incrGEmbeddedDomain; // one cell thick layer around and outside of the embedded domain
    for (unsigned int d = 0; d < DIM; ++ d) {
        incrGEmbeddedDomain[d] = Index(gEmbeddedDomain[d].first() - 1 + positionEmbedded[d],
                                       gEmbeddedDomain[d].last() + 1 + positionEmbedded[d]);
    }
    std::vector<NDIndex<DIM> > localEmbeddedDomains;
    std::vector<NDIndex<DIM> > localSurroundingDomains;
    Utils::getLocalDomains(embedded.getLayout(), localEmbeddedDomains);
    Utils::getLocalDomains(surrounding.getLayout(), localSurroundingDomains);

    std::vector<NDIndex<DIM> > lDoms[] = {localEmbeddedDomains, localSurroundingDomains};
    Vektor<int,DIM>  diffs[] = {-positionEmbedded, positionEmbedded};

    std::vector<MPI_Request> requests;
    int tag = Ippl::Comm->next_tag(F_GUARD_CELLS_TAG, F_TAG_CYCLE);

    // the boundary layer is treated as 2 * DIM one-dimensional domains
    // compute the global domains of these one-dimensional fields
    NDIndex<DIM> gDoms[4 * DIM]; // global boundary layers
    for (unsigned int d = 0; d < DIM; ++ d) {
        gDoms[d] = gEmbeddedDomain;
        gDoms[d][d] = Index(gDoms[d][d].first(), gDoms[d][d].first());

        gDoms[d + DIM] = gEmbeddedDomain;
        gDoms[d + DIM][d] = Index(gDoms[d + DIM][d].last(), gDoms[d + DIM][d].last());

        gDoms[d + 2 * DIM] = incrGEmbeddedDomain;
        gDoms[d + 2 * DIM][d] = Index(gDoms[d + 2 * DIM][d].first(), gDoms[d + 2 * DIM][d].first());

        gDoms[d + 3 * DIM] = incrGEmbeddedDomain;
        gDoms[d + 3 * DIM][d] = Index(gDoms[d + 3 * DIM][d].last(), gDoms[d + 3 * DIM][d].last());
    }

    std::vector<std::vector<char> > sendData(numNodes);
    std::vector<std::vector<char> > receiveData(numNodes);
    // determine the local contribution to the boundary layers and fill in the values
    for (unsigned int k = 0; k < numNodes; ++ k) {
        if (k == myNode) continue;

        auto & data = receiveData[k];

        for (unsigned int side = 0; side < 4 * DIM; ++ side) {
            unsigned int test = 1 << side;
            if ((op.sides & test) > 0) {
                NDIndex<DIM> dom;
                if (side >= 2 * DIM) {
                    // from surrounding to embedded
                    for (unsigned int d2 = 0; d2 < DIM; ++ d2) {
                        int lower = std::max(localSurroundingDomains[k][d2].first(), gDoms[side][d2].first());
                        int upper = std::min(localSurroundingDomains[k][d2].last(),  gDoms[side][d2].last());
                        int sign = lower <= upper? 1: upper - lower;
                        dom[d2] = Index(lower, upper, sign);
                    }
                } else {
                    // from embedded to surrounding
                    for (unsigned int d2 = 0; d2 < DIM; ++ d2) {
                        int lower = std::max(localEmbeddedDomains[k][d2].first(), gDoms[side][d2].first());
                        int upper = std::min(localEmbeddedDomains[k][d2].last(),  gDoms[side][d2].last());
                        int sign = lower <= upper? 1: upper - lower;
                        dom[d2] = Index(lower, upper, sign);
                    }
                }

                bool hasSize = true;
                for (unsigned int d2 = 0; d2 < DIM; ++ d2) {
                    if (dom[d2].first() > dom[d2].last()) {
                        hasSize = false;
                        break;
                    }
                }
                if (!hasSize) continue;

                // if from embedded to surrounding then idx = 0 otherwise idx = 1
                unsigned int idx = (side < 2 * DIM? 0: 1);

                bool inside = true;
                const NDIndex<DIM> & flddom = lDoms[1-idx][myNode];
                NDIndex<DIM> overlap;

                for (unsigned int d2 = 0; d2 < DIM; ++ d2) {
                    int lower = std::max(dom[d2].first(), flddom[d2].first() + diffs[idx][d2] - 1);
                    int upper = std::min(dom[d2].last(),  flddom[d2].last() + diffs[idx][d2] + 1);
                    if (lower > upper) {
                        inside = false;
                        break;
                    }
                    overlap[d2] = Index(lower, upper);
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

        for (unsigned int side = 0; side < 4 * DIM; ++ side) {
            unsigned int test = 1 << side;
            if ((op.sides & test) == 0) continue;
            NDIndex<DIM> dom;
            if (side >= 2 * DIM) {
                // from surrounding to embedded
                for (unsigned int d2 = 0; d2 < DIM; ++ d2) {
                    int lower = std::max(localSurroundingDomains[myNode][d2].first(), gDoms[side][d2].first());
                    int upper = std::min(localSurroundingDomains[myNode][d2].last(),  gDoms[side][d2].last());
                    int sign = lower <= upper? 1: upper - lower;
                    dom[d2] = Index(lower, upper, sign);
                }
            } else {
                // from embedded to surrounding
                for (unsigned int d2 = 0; d2 < DIM; ++ d2) {
                    int lower = std::max(localEmbeddedDomains[myNode][d2].first(), gDoms[side][d2].first());
                    int upper = std::min(localEmbeddedDomains[myNode][d2].last(),  gDoms[side][d2].last());
                    int sign = lower <= upper? 1: upper - lower;
                    dom[d2] = Index(lower, upper, sign);
                }
            }

            bool hasSize = true;
            for (unsigned int d2 = 0; d2 < DIM; ++ d2) {
                if (dom[d2].first() > dom[d2].last()) {
                    hasSize = false;
                    break;
                }
            }
            if (!hasSize) continue;

            unsigned int idx = (side < 2 * DIM? 0: 1);

            bool inside = true;
            const NDIndex<DIM> & flddom = lDoms[1-idx][k];
            NDIndex<DIM> overlap;

            for (unsigned int d2 = 0; d2 < DIM; ++ d2) {
                int lower = std::max(dom[d2].first(), flddom[d2].first() + diffs[idx][d2] - 1);
                int upper = std::min(dom[d2].last(),  flddom[d2].last() + diffs[idx][d2] + 1);
                if (lower > upper) {
                    inside = false;
                    break;
                }
                overlap[d2] = Index(lower, upper);
            }
            if (!inside) continue;

            unsigned int sizeData = overlap.size();
            data.reserve(data.size() + 1 + 2 * DIM * sizeof(int) + MD5_DIGEST_LENGTH + sizeData * sizeof(double));
            const char* buffer;
            Field<T_t, DIM, Mesh_t, C_t> * fld = NULL;
            if (side < 2 * DIM) {
                char FType = 0;
                data.insert(data.end(), &FType, &FType + 1);

                for (unsigned int d2 = 0; d2 < DIM; ++ d2) {
                    int lower = overlap[d2].first() + positionEmbedded[d2];
                    buffer = reinterpret_cast<const char*>(&lower);
                    data.insert(data.end(), buffer, buffer + sizeof(int));

                    int upper = overlap[d2].last() + positionEmbedded[d2];
                    buffer = reinterpret_cast<const char*>(&upper);
                    data.insert(data.end(), buffer, buffer + sizeof(int));
                }

                fld = &embedded;
            } else {
                char FType = 1;
                data.insert(data.end(), &FType, &FType + 1);

                for (unsigned int d2 = 0; d2 < DIM; ++ d2) {
                    int lower = overlap[d2].first() - positionEmbedded[d2];
                    buffer = reinterpret_cast<const char*>(&lower);
                    data.insert(data.end(), buffer, buffer + sizeof(int));

                    int upper = overlap[d2].last() - positionEmbedded[d2];
                    buffer = reinterpret_cast<const char*>(&upper);
                    data.insert(data.end(), buffer, buffer + sizeof(int));
                }

                fld = &surrounding;
            }

            double * rawData = new double[sizeData];
            unsigned int ii = 0;
            for (int j = overlap[1].first(); j <= overlap[1].last(); ++ j) {
                elem[1] = Index(j,j);
                for (int i = overlap[0].first(); i <= overlap[0].last(); ++ i) {
                    elem[0] = Index(i,i);
                    op.execute(*fld, elem, rawData[ii++]);
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
            char FType = data[it ++];
            Field<T_t, DIM, Mesh_t, C_t> * fld = (FType == 0? &surrounding: &embedded);

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
            memcpy(sentMd5Hash, &data[it], MD5_DIGEST_LENGTH);
            it += MD5_DIGEST_LENGTH;

            MD5(reinterpret_cast<const unsigned char*>(&data[it]), sizeData * sizeof(double), receivedMd5Hash);
            char rcMD5Buffer[2 * MD5_DIGEST_LENGTH];
            char snMD5Buffer[2 * MD5_DIGEST_LENGTH];
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

            int ii = 0;
            for (int j = overlap[1].first(); j <= overlap[1].last(); ++ j) {
                elem[1] = Index(j,j);
                for (int i = overlap[0].first(); i <= overlap[0].last(); ++ i) {
                    elem[0] = Index(i,i);

                    op.invert(*fld, elem, rawData[ii++]);
                }
            }

            delete[] rawData;
        }
    }

    embedded.fillGuardCells();
    surrounding.fillGuardCells();
}

template<class T, class C, class Op>
void Communicator::serializeField(double * rawData,
                                  const Field<T, DIM, Mesh_t, C> & fld,
                                  const NDIndex<DIM> & dom,
                                  const Op & op)
{
    NDIndex<DIM> elem;
    unsigned int ii = 0;
    for (int j = dom[1].first(); j <= dom[1].last(); ++ j) {
        elem[1] = Index(j,j);
        for (int i = dom[0].first(); i <= dom[0].last(); ++ i) {
            elem[0] = Index(i,i);
            op.execute(fld, elem, rawData[ii++]);
        }
    }
}

template<class T, class C, class Op>
void Communicator::deserializeField(Field<T, DIM, Mesh_t, C> & fld,
                                    const double * rawData,
                                    const NDIndex<DIM> & dom,
                                    const Op & op)
{
    NDIndex<DIM> elem;
    unsigned int ii = 0;
    for (int j = dom[1].first(); j <= dom[1].last(); ++ j) {
        elem[1] = Index(j,j);
        for (int i = dom[0].first(); i <= dom[0].last(); ++ i) {
            elem[0] = Index(i,i);
            op.invert(fld, elem, rawData[ii++]);
        }
    }
}

inline
const std::vector<unsigned int> & Communicator::getNeighbours()
{
    initialize();
    return _neighbours;
}

inline
const std::vector<unsigned int> & Communicator::getNode2Index()
{
    initialize();
    return _node2index;
}

inline
const std::vector<unsigned int> & Communicator::getIndex2Node()
{
    initialize();
    return _index2node;
}

inline
std::pair<unsigned int, unsigned int> Communicator::getNodeMeshSize()
{
    initialize();
    return std::make_pair(_nNodeCols, _nNodeRows);
}

template <class T>
void Communicator::communicateFields(FieldPatch<T> & fp,
                                     const NDIndex<DIM> & dom,
                                     const Timings::TimerRef & commTimer,
                                     const bool all2all)
{
    // extern std::ofstream dbg;
#ifdef P2PCOMM
    NDIndex<DIM> lPDom = fp.getDomain();
    // dbg << "Communicator.hh: " << __LINE__ << "\t" << lPDom << "\t" << dom << std::endl;
    auto localPDomains = getAllLocalDomains(lPDom);
    auto localFDomains = getAllLocalDomains(dom);

    p2pCommunicateFields(fp, localPDomains, localFDomains, commTimer, all2all);
#else
    collectivelyCommunicateFields(fp, dom, commTimer);
#endif
}

template <class T>
void Communicator::communicateFields(FieldPatch<T> & fp,
                                     const std::vector<NDIndex<DIM> >  & localPDomains,
                                     const std::vector<NDIndex<DIM> >  & localFDomains,
                                     const Timings::TimerRef & commTimer,
                                     const bool all2all)
{
    p2pCommunicateFields(fp, localPDomains, localFDomains, commTimer, all2all);
}

template <class T>
void Communicator::collectivelyCommunicateFields(FieldPatch<T> & fp,
                                                 const NDIndex<DIM> & dom,
                                                 const Timings::TimerRef & commTimer)
{
    extern std::ofstream dbg;
    boost::mpi::communicator world;

    Timings::startTimer(commTimer);
    std::vector<FieldPatch<T> > patches;
    dbg << "Communicator.hh: " << __LINE__ << std::endl;
    boost::mpi::all_gather(world, fp, patches);

    fp.resize(dom);

    for (size_t i = 0; i < patches.size(); ++ i) {
        if (i == (size_t) Ippl::myNode() || patches[i].size() == 0) {
            patches[i].clear();
            continue;
        }
        fp.addWithoutResize(patches[i]);
        patches[i].clear();
    }
    Timings::stopTimer(commTimer);
}

template <class T>
void Communicator::p2pCommunicateFields(FieldPatch<T> & fp,
                                        const std::vector<NDIndex<DIM> >  & localPDomains,
                                        const std::vector<NDIndex<DIM> >  & localFDomains,
                                        const Timings::TimerRef & commTimer,
                                        const bool all2all)
{
    boost::mpi::communicator world;
    initialize();

    const unsigned int myNode = Ippl::myNode();
    const unsigned int numNodes = Ippl::getNodes();
    Timings::startTimer(commTimer);

    std::vector<unsigned int> partners;
    if (all2all) {
        partners.resize(numNodes - 1);
        for (unsigned int k = 0; k < numNodes - 1; ++ k)
            partners[k] = k < myNode? k: k + 1;
    } else {
        partners.assign(_neighbours.begin(), _neighbours.end());
    }

    std::vector<FieldPatch<T> > mySplitPatches(partners.size());
    std::vector<FieldPatch<T> > othersSplitPatches(partners.size());

    int totalmsgsend = 0;
    std::vector<boost::mpi::request> requests;
    int tag = Ippl::Comm->next_tag(F_GUARD_CELLS_TAG, F_TAG_CYCLE);

    size_t totalmsgrecv = 0;
    for (unsigned int k: partners) {
        bool inside = true;
        NDIndex<DIM> dom;
        for (unsigned int d = 0; d < DIM; ++ d) {
            int lower = std::max(localPDomains[k][d].first(), localFDomains[myNode][d].first());
            int upper = std::min(localPDomains[k][d].last(),  localFDomains[myNode][d].last());
            if (lower > upper) {
                inside = false;
                break;
            }
            dom[d] = Index(lower, upper);
        }
        if (!inside) continue;

        boost::mpi::request req = world.irecv(k, tag, &othersSplitPatches[totalmsgrecv], 1);
        ++ totalmsgrecv;
        requests.push_back(req);
    }

    for (unsigned int k: partners){
        NDIndex<DIM> dom;
        bool inside = true;
        for (unsigned int d = 0; d < DIM; ++ d) {
            int lower = std::max(localPDomains[myNode][d].first(), localFDomains[k][d].first());
            int upper = std::min(localPDomains[myNode][d].last(),  localFDomains[k][d].last());
            if (lower > upper) {
                inside = false;
                break;
            }
            dom[d] = Index(lower, upper);
        }
        if (!inside) continue;

        mySplitPatches[totalmsgsend].resize(dom);
        mySplitPatches[totalmsgsend].addWithoutResize(fp);

        // world.send(k, tag, &mySplitPatches[totalmsgsend], 1);
        boost::mpi::request req = world.isend(k, tag, &mySplitPatches[totalmsgsend], 1);
        requests.push_back(req);
        ++ totalmsgsend;
    }

    fp.resize(localFDomains[myNode]);

    boost::mpi::wait_all(&requests[0], &requests[0] + requests.size());

    for (size_t i = 0; i < totalmsgrecv; ++ i) {

        fp.addWithoutResize(othersSplitPatches[i]);
    }

    Timings::stopTimer(commTimer);
}

template<class T>
BoundaryLayer<T>::BoundaryLayer(BLType type, NDIndex<DIM> dom):
    FieldPatch<T>(dom),
    _bltype(type)
{ }

template<class T>
BoundaryLayer<T>::BoundaryLayer():
    FieldPatch<T>()
{ }

template<class T>
BoundaryLayer<T>::BoundaryLayer(const BoundaryLayer<T>& bl):
    FieldPatch<T>(bl),
    _bltype(bl._bltype)
{ }

template<class T>
template<class Archive>
void BoundaryLayer<T>::save(Archive& ar, const unsigned int version) const
{
    FieldPatch<T>::save(ar, version);

    ar & _bltype;
}

template<class T>
template<class Archive>
void BoundaryLayer<T>::load(Archive& ar, const unsigned int version)
{
    FieldPatch<T>::load(ar, version);

    ar & _bltype;
}

struct HFieldToDouble {
    static
    const unsigned int sides;

    template<class T, class C>
    void execute(const Field<T, DIM, Mesh_t, C> & fld,
                 const NDIndex<DIM> & elem,
                 double & array) const
    {
        const Vector_t & value = fld.localElement(elem);
        array = value[0] + value[1];
    }

    template<class T, class C>
    void invert(Field<T, DIM, Mesh_t, C> & fld,
                const NDIndex<DIM> & elem,
                const double & array) const
    {
        Vector_t & value = fld.localElement(elem);
        value[0] = 0.5 * array;
        value[1] = 0.5 * array;
    }
};

struct ExFieldToDouble {
    static
    const unsigned int sides;

    template<class T, class C>
    void execute(const Field<T, DIM, Mesh_t, C> & fld,
                 const NDIndex<DIM> & elem,
                 double & array) const
    {
        const Vector_t & value = fld.localElement(elem);
        array = value[0];
    }

    template<class T, class C>
    void invert(Field<T, DIM, Mesh_t, C> & fld,
                const NDIndex<DIM> & elem,
                const double & array) const
    {
        Vector_t & value = fld.localElement(elem);
        value[0] = array;
    }
};

struct EyFieldToDouble {
    static
    const unsigned int sides;

    template<class T, class C>
    void execute(const Field<T, DIM, Mesh_t, C> & fld,
                 const NDIndex<DIM> & elem,
                 double & array) const
    {
        const Vector_t & value = fld.localElement(elem);
        array = value[1];
    }

    template<class T, class C>
    void invert(Field<T, DIM, Mesh_t, C> & fld,
                const NDIndex<DIM> & elem,
                const double & array) const
    {
        Vector_t & value = fld.localElement(elem);
        value[1] = array;
    }
};

#endif
