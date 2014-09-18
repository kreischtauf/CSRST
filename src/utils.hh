#ifndef UTILS_HH
#define UTILS_HH

#include "defs.hh"

namespace Utils {

    void getLocalDomains(const FieldLayout<DIM> & FLayout,
                         std::vector<NDIndex<DIM> > & vectorLocalDomains);

    void getNodeIndices(const FieldLayout<DIM> & cellCenteredLayout,
                        std::vector<int> & vectorNodeIndices);

    void increaseLocalDomainsLast(const FieldLayout<DIM> & cellCenteredLayout,
                                  std::vector<NDIndex<DIM> > & vectorLocalDomains);

    void decreaseLocalDomainsLast(const FieldLayout<DIM> & cellCenteredLayout,
                                  std::vector<NDIndex<DIM> > & vectorLocalDomains);

    void increaseLocalDomainsFirst(const FieldLayout<DIM> & cellCenteredLayout,
                                   std::vector<NDIndex<DIM> > & vectorLocalDomains);

    void decreaseLocalDomainsFirst(const FieldLayout<DIM> & cellCenteredLayout,
                                   std::vector<NDIndex<DIM> > & vectorLocalDomains);

    void addGostCellToLocalDomains(const FieldLayout<DIM> & cellCenteredLayout,
                                   std::vector<NDIndex<DIM> > & vectorLocalDomains);

    FieldLayout_Cell_t *
    getCellFromVertCenteredLayout(UniformCartesian<DIM> & mesh, FieldLayout<DIM> & CenteredLayout);
    FieldLayout_Cell_t *
    getCellFromEdgeCenteredLayout(UniformCartesian<DIM> & mesh, FieldLayout<DIM> & CenteredLayout);
    FieldLayout_Vert_t *
    getVertFromCellCenteredLayout(UniformCartesian<DIM> & mesh, FieldLayout<DIM> & CenteredLayout);
    FieldLayout_Vert_t *
    getVertFromEdgeCenteredLayout(UniformCartesian<DIM> & mesh, FieldLayout<DIM> & CenteredLayout);
    FieldLayout_Edge_t *
    getEdgeFromVertCenteredLayout(UniformCartesian<DIM> & mesh, FieldLayout<DIM> & CenteredLayout);
    FieldLayout_Edge_t *
    getEdgeFromCellCenteredLayout(UniformCartesian<DIM> & mesh, FieldLayout<DIM> & CenteredLayout);

    SField_Vert_ptr getScalarVertField(const VField_Vert_t & vectorField);
    SField_Vert_ptr getScalarVertField(const VField_Cell_t & vectorField);
    SField_Vert_ptr getScalarVertField(const VField_Edge_t & vectorField);
    SField_Cell_ptr getScalarCellField(const VField_Vert_t & vectorField);
    SField_Cell_ptr getScalarCellField(const VField_Cell_t & vectorField);
    SField_Cell_ptr getScalarCellField(const VField_Edge_t & vectorField);

    template<class M, class T>
    class CompanionDeleter {
    public:
        CompanionDeleter(T* t);
        void operator()(M* m);
    private:
        T* companion_;
    };

    struct CompanionSorterItem {
        double fraction;
        Vector_t subStep;
    };

    inline
    bool CompanionSorter(const CompanionSorterItem & a,
                         const CompanionSorterItem & b)
    {
        return a.fraction < b.fraction;
    }

    class samePos {
    public:
        samePos(const double & ref):
            _ref(ref) { }

        bool operator() (const Vector_t & p1, const Vector_t & p2) {
            return (std::abs(p1(0) - p2(0)) < 1e-6 * _ref &&
                    std::abs(p1(1) - p2(1)) < 1e-6 * _ref);
        }

        bool operator() (const CompanionSorterItem & a,
                         const CompanionSorterItem & b) {
            return (std::abs(a.fraction - b.fraction) < 1e-6);
        }

    private:
        double _ref;
    };

    class ascDistance {
    public:
        ascDistance(Vector_t ref):
            _ref(ref) { }

        bool operator()(const Vector_t & p1, const Vector_t & p2) {
            const bool retval = VecDistance(_ref, p1) < VecDistance(_ref, p2);
            return retval;
        }
    private:
        double VecDistance(const Vector_t & a, const Vector_t & b)
        {
#if DIM == 2
            return sqrt((a(0) - b(0))*(a(0) - b(0)) +
                        (a(1) - b(1))*(a(1) - b(1)));
#elif DIM == 3
            return sqrt((a(0) - b(0))*(a(0) - b(0)) +
                        (a(1) - b(1))*(a(1) - b(1)) +
                        (a(2) - b(2))*(a(2) - b(2)));
#else
            return 0.0;
#endif
        }

        Vector_t _ref;
    };

    class XFirst {
    public:
        XFirst(const Mesh_t & mesh):
            _dx(mesh.get_meshSpacing(0),
                mesh.get_meshSpacing(1)),
            _origin(mesh.get_origin())
        { }

        bool operator() (const Vector_t & v1, const Vector_t & v2) {
#if DIM == 3
            int K1 = (v1(2) - _origin(2)) / _dx(2);
            int K2 = (v2(2) - _origin(2)) / _dx(2);
            if (K2 > K1) return true;
            if (K2 < K1) return false;
#else
            int K1, K2;
#endif
            K1 = (v1(1) - _origin(1)) / _dx(1);
            K2 = (v2(1) - _origin(1)) / _dx(1);
            if (K2 > K1) return true;
            if (K2 < K1) return false;

            K1 = (v1(0) - _origin(0)) / _dx(0);
            K2 = (v2(0) - _origin(0)) / _dx(0);
            if (v2(0) >= v1(0)) return true;
            return false;
        }
    private:
        Vector_t _dx;
        Vector_t _origin;
    };


    size_t calcDualCellBorderCrossings(CompanionSorterItem xy_crossings[],
                                       Vector_t oldR,
                                       Vector_t newR,
                                       const Vector_t & spacing,
                                       const Vector_t & origin,
                                       const Utils::samePos & uni);
    size_t calcCellBorderCrossings(CompanionSorterItem xy_crossings[],
                                   Vector_t oldR,
                                   Vector_t newR,
                                   const Vector_t & spacing,
                                   const Vector_t & origin,
                                   const Utils::samePos & uni);

    int ipow(int base, unsigned int exp);

    void process_mem_usage(double& vm_usage, double& resident_set);

    float Q_rsqrt( float x );

    struct color_rgb {
        color_rgb() {
            vec[0] = vec[1] = vec[2] = 0.0;
        }

        color_rgb(float x, float y, float z){
            vec[0] = x;
            vec[1] = y;
            vec[2] = z;
        }

        float & operator[](int i) {
            return vec[i];
        }
        color_rgb operator*(float t) {
            return color_rgb(t*vec[0],t*vec[1],t*vec[2]);
        }

        void operator+=(const color_rgb& v2) {
            vec[0] += v2.vec[0];
            double maxEntry = vec[0];
            vec[1] += v2.vec[1];
            maxEntry = vec[1] > maxEntry ? vec[1]: maxEntry;
            vec[2] += v2.vec[2];
            maxEntry = vec[2] > maxEntry ? vec[2]: maxEntry;

            for (int i = 0; i < 3; ++ i) {
                vec[i] /= maxEntry;
            }
        }

        float vec[3];
    };

    struct KahanAccumulation
    {
        long double sum;
        long double correction;
        KahanAccumulation();
    };

    KahanAccumulation KahanSum(KahanAccumulation accumulation, double value);
}

template<class M, class T>
Utils::CompanionDeleter<M,T>::CompanionDeleter(T* t):
companion_(t)
{ }

template<class M, class T>
void
Utils::CompanionDeleter<M,T>::operator()(M* master)
{
    delete master;
    delete companion_;
}

inline
void Utils::getLocalDomains(const FieldLayout<DIM> & FLayout,
                            std::vector<NDIndex<DIM> > & vectorLocalDomains)
{
    vectorLocalDomains.resize(Ippl::getNodes());
    vectorLocalDomains[Ippl::myNode()] = FLayout.getLocalNDIndex();
    for (FieldLayout<DIM>::const_iterator_dv domainMapIterator = FLayout.begin_rdv();
         domainMapIterator != FLayout.end_rdv();
         ++ domainMapIterator) {
        RefCountedP<Vnode<DIM> > virtualNode = (*domainMapIterator).second;
        vectorLocalDomains[virtualNode->getNode()] = virtualNode->getDomain();
    }
}

inline
void Utils::getNodeIndices(const FieldLayout<DIM> & cellCenteredLayout,
                           std::vector<int> & vectorNodeIndices)
{
    vectorNodeIndices.resize(Ippl::getNodes());
    vectorNodeIndices[Ippl::myNode()] = Ippl::myNode();
    for (FieldLayout<DIM>::const_iterator_dv domainMapIterator = cellCenteredLayout.begin_rdv();
         domainMapIterator != cellCenteredLayout.end_rdv();
         ++ domainMapIterator)
        {
            RefCountedP<Vnode<DIM> > virtualNode = (*domainMapIterator).second;
            vectorNodeIndices[virtualNode->getNode()] = virtualNode->getNode();
        }
}

inline
void Utils::increaseLocalDomainsLast(const FieldLayout<DIM> & cellCenteredLayout,
                                     std::vector<NDIndex<DIM> > & vectorLocalDomains)
{
    const NDIndex<DIM> & globalDomain = cellCenteredLayout.getDomain();

    for (std::vector<NDIndex<DIM> >::iterator NDIndexIterator = vectorLocalDomains.begin();
         NDIndexIterator != vectorLocalDomains.end();
         ++ NDIndexIterator) {
        NDIndex<DIM> & domain = *NDIndexIterator;
        if (domain[0].last() == globalDomain[0].last()) {
            domain[0] = Index(domain[0].first(), globalDomain[0].last() + 1);
        }

        if (domain[1].last() == globalDomain[1].last()) {
            domain[1] = Index(domain[1].first(), globalDomain[1].last() + 1);
        }
    }
}

inline
void Utils::decreaseLocalDomainsLast(const FieldLayout<DIM> & cellCenteredLayout,
                                     std::vector<NDIndex<DIM> > & vectorLocalDomains)
{
    const NDIndex<DIM> & globalDomain = cellCenteredLayout.getDomain();

    for (std::vector<NDIndex<DIM> >::iterator NDIndexIterator = vectorLocalDomains.begin();
         NDIndexIterator != vectorLocalDomains.end();
         ++ NDIndexIterator) {
        NDIndex<DIM> & domain = *NDIndexIterator;
        if (domain[0].last() == globalDomain[0].last()) {
            domain[0] = Index(domain[0].first(), globalDomain[0].last() - 1);
        }

        if (domain[1].last() == globalDomain[1].last()) {
            domain[1] = Index(domain[1].first(), globalDomain[1].last() - 1);
        }
    }
}

inline
void Utils::increaseLocalDomainsFirst(const FieldLayout<DIM> & cellCenteredLayout,
                                      std::vector<NDIndex<DIM> > & vectorLocalDomains)
{
    const NDIndex<DIM> & globalDomain = cellCenteredLayout.getDomain();

    for (std::vector<NDIndex<DIM> >::iterator NDIndexIterator = vectorLocalDomains.begin();
         NDIndexIterator != vectorLocalDomains.end();
         ++ NDIndexIterator) {
        NDIndex<DIM> & domain = *NDIndexIterator;
        if (domain[0].first() == globalDomain[0].first()) {
            domain[0] = Index(globalDomain[0].first() - 1, domain[0].last());
        }

        if (domain[1].first() == globalDomain[1].first()) {
            domain[1] = Index(globalDomain[1].first() - 1, domain[1].last());
        }
    }
}

inline
void Utils::decreaseLocalDomainsFirst(const FieldLayout<DIM> & cellCenteredLayout,
                                      std::vector<NDIndex<DIM> > & vectorLocalDomains)
{
    const NDIndex<DIM> & globalDomain = cellCenteredLayout.getDomain();

    for (std::vector<NDIndex<DIM> >::iterator NDIndexIterator = vectorLocalDomains.begin();
         NDIndexIterator != vectorLocalDomains.end();
         ++ NDIndexIterator) {
        NDIndex<DIM> & domain = *NDIndexIterator;
        if (domain[0].first() == globalDomain[0].first()) {
            domain[0] = Index(globalDomain[0].first() + 1, domain[0].last());
        }

        if (domain[1].first() == globalDomain[1].first()) {
            domain[1] = Index(globalDomain[1].first() + 1, domain[1].last());
        }
    }
}


inline
void Utils::addGostCellToLocalDomains(const FieldLayout<DIM> & cellCenteredLayout,
                                      std::vector<NDIndex<DIM> > & vectorLocalDomains)
{
    const NDIndex<DIM> & globalDomain = cellCenteredLayout.getDomain();

    for (std::vector<NDIndex<DIM> >::iterator NDIndexIterator = vectorLocalDomains.begin();
         NDIndexIterator != vectorLocalDomains.end();
         ++ NDIndexIterator) {
        NDIndex<DIM> & domain = *NDIndexIterator;

        if (domain[0].first() != globalDomain[0].first()) {
            domain[0] = Index(domain[0].first() - 1, domain[0].last());
        }

        if (domain[0].last() != globalDomain[0].last()) {
            domain[0] = Index(domain[0].first(), domain[0].last() + 1);
        }

        if (domain[1].first() != globalDomain[1].first()) {
            domain[1] = Index(domain[1].first() - 1, domain[1].last());
        }

        if (domain[1].last() != globalDomain[1].last()) {
            domain[1] = Index(domain[1].first(), domain[1].last() + 1);
        }
    }
}

inline
size_t Utils::calcDualCellBorderCrossings(CompanionSorterItem xy_crossings[],
                                          Vector_t oldR,
                                          Vector_t newR,
                                          const Vector_t & spacing,
                                          const Vector_t & origin,
                                          const Utils::samePos & uni)
{
    oldR -= origin;
    newR -= origin;
    int Ii = static_cast<int>(floor(oldR(0) / spacing(0) - 0.5));
    int Ji = static_cast<int>(floor(oldR(1) / spacing(1) - 0.5));
    int If = static_cast<int>(floor(newR(0) / spacing(0) - 0.5));
    int Jf = static_cast<int>(floor(newR(1) / spacing(1) - 0.5));

    int I = Ii;
    int J = Ji;
    Vector_t Dx = newR - oldR;

    size_t next = 0;
    xy_crossings[next++].fraction = 0.0;
    // xy_crossings[next++].subStep = oldR;
    while (I < If) {
        xy_crossings[next].fraction = (spacing(0) * (I + 1.5) - oldR(0)) / Dx(0);
        xy_crossings[next].subStep = Vector_t(spacing(0) * (I + 1.5),
                                              oldR(1) + xy_crossings[next].fraction * Dx(1));
        ++ I;
        ++ next;
    }
    while (I > If) {
        xy_crossings[next].fraction = (spacing(0) * (I + 0.5) - oldR(0)) / Dx(0);
        xy_crossings[next].subStep = Vector_t(spacing(0) * (I + 0.5),
                                              oldR(1) + xy_crossings[next].fraction * Dx(1));
        -- I;
        ++ next;
    }

    while (J < Jf) {
        xy_crossings[next].fraction = (spacing(1) * (J + 1.5) - oldR(1)) / Dx(1);
        xy_crossings[next].subStep = Vector_t(oldR(0) + xy_crossings[next].fraction * Dx(0),
                                              spacing(1) * (J + 1.5));
        ++ J;
        ++ next;
    }
    while (J > Jf) {
        xy_crossings[next].fraction = (spacing(1) * (J + 0.5) - oldR(1)) / Dx(1);
        xy_crossings[next].subStep = Vector_t(oldR(0) + xy_crossings[next].fraction * Dx(0),
                                              spacing(1) * (J + 0.5));
        -- J;
        ++ next;
    }
    xy_crossings[next++].fraction = 1.0;
    // xy_crossings[next++].subStep = newR;

    std::sort(xy_crossings + 1, xy_crossings + next - 1, CompanionSorter);
    auto end = std::unique(xy_crossings, xy_crossings + next, uni);
    size_t num = (end - xy_crossings);

    xy_crossings[0].subStep = oldR;
    xy_crossings[num - 1].subStep = newR;
    return num;
}

inline
size_t Utils::calcCellBorderCrossings(CompanionSorterItem xy_crossings[],
                                      Vector_t oldR,
                                      Vector_t newR,
                                      const Vector_t & spacing,
                                      const Vector_t & origin,
                                      const Utils::samePos & uni)
{
    oldR -= origin;
    newR -= origin;
    const int Ii = static_cast<int>(floor(oldR(0) / spacing(0)));
    const int Ji = static_cast<int>(floor(oldR(1) / spacing(1)));
    const int If = static_cast<int>(floor(newR(0) / spacing(0)));
    const int Jf = static_cast<int>(floor(newR(1) / spacing(1)));

    int I = Ii;
    int J = Ji;
    const Vector_t Dx = newR - oldR;

    size_t next = 0;
    xy_crossings[next++].fraction = 0.0;
    // xy_crossings[next++].subStep = oldR;
    while (I < If) {
        xy_crossings[next].fraction = (spacing(0) * (I + 1) - oldR(0)) / Dx(0);
        xy_crossings[next].subStep = Vector_t(spacing(0) * (I + 1),
                                              oldR(1) + xy_crossings[next].fraction * Dx(1));
        ++ I;
        ++ next;
    }
    while (I > If) {
        xy_crossings[next].fraction = (spacing(0) * I - oldR(0)) / Dx(0);
        xy_crossings[next].subStep = Vector_t(spacing(0) * I,
                                              oldR(1) + xy_crossings[next].fraction * Dx(1));
        -- I;
        ++ next;
    }

    while (J < Jf) {
        xy_crossings[next].fraction = (spacing(1) * (J + 1) - oldR(1)) / Dx(1);
        xy_crossings[next].subStep = Vector_t(oldR(0) + xy_crossings[next].fraction * Dx(0),
                                              spacing(1) * (J + 1));
        ++ J;
        ++ next;
    }
    while (J > Jf) {
        xy_crossings[next].fraction = (spacing(1) * J - oldR(1)) / Dx(1);
        xy_crossings[next].subStep = Vector_t(oldR(0) + xy_crossings[next].fraction * Dx(0),
                                              spacing(1) * J);
        -- J;
        ++ next;
    }
    xy_crossings[next++].fraction = 1.0;
    // xy_crossings[next++].subStep = newR;

    std::sort(xy_crossings + 1, xy_crossings + next - 1, CompanionSorter);
    auto end = std::unique(xy_crossings, xy_crossings + next, uni);
    size_t num = (end - xy_crossings);

    xy_crossings[0].subStep = oldR;
    xy_crossings[num - 1].subStep = newR;

    return num;
}

inline
float Utils::Q_rsqrt( float x )
{
    float xhalf = 0.5f*x;
    float xqhalf = 0.5010936f*x;
    uint32_t i;
    assert(sizeof(x) == sizeof(i));
    std::memcpy(&i, &x, sizeof(i));
    i = 0x5f3759df - (i>>1);
    std::memcpy(&x, &i, sizeof(i));
    x = x*(1.5f - xhalf*x*x);        // result is less than correct result for all x since 1/sqrt(x) is convex
    x = x*(1.5010936f - xqhalf*x*x); // optimized slope (steeper slope pushes result to higher values)
    return x;
}

inline
Utils::KahanAccumulation::KahanAccumulation():
    sum(0.0),
    correction(0.0)
{ }

inline
Utils::KahanAccumulation Utils::KahanSum(KahanAccumulation accumulation, double value)
{
    Utils::KahanAccumulation result;
    long double y = value - accumulation.correction;
    long double t = accumulation.sum + y;
    result.correction = (t - accumulation.sum) - y;
    result.sum = t;
    return result;
}

#endif
