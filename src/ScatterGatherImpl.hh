/***************************************************************************
                         ScatterGatherImpl.hh
                         -------------------
    begin                : Wed Aug 7 2013
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

#ifndef SCATTERGATHERIMPL_HH
#define SCATTERGATHERIMPL_HH

#include "defs.hh"
#include "FieldPatch.hh"
#include "Physics.hh"

#include <vector>

template<class IntOp, class T>
class ScatterGatherImpl {
public:
    // gather edge centered FieldPatch at particle location
    static
    void gather(const FieldPatch<T> & fp,
                ParticleAttrib< Vector_t > & fld,
                const ParticleAttrib< Vector_t > & R,
                const size_t & comp,
                size_t from,
                size_t to)
    { return; }

    // gather cell centered FieldPatch at particle location
    static
    void gather(const FieldPatch<T> & fp,
                ParticleAttrib< double > & fld,
                const ParticleAttrib< Vector_t > & R,
                size_t from,
                size_t to)
    { return; }

    // scatter edge centered ParticleAttrib onto FieldPatch
    static
    void scatter(FieldPatch<T> & fp,
                 ParticleAttrib< Vector_t > & fld,
                 const ParticleAttrib< Vector_t > & R,
                 const size_t & comp)
    { return; }

    // scatter cell centered ParticleAttrib onto FieldPatch
    static
    void scatter(FieldPatch<T> & fp,
                 ParticleAttrib< Vector_t > & fld,
                 const ParticleAttrib< Vector_t > & R)
    { return; }

};

template<>
class ScatterGatherImpl<IntCIC, double> {
public:
    static
    void gather(const FieldPatch<double> & fp,
                ParticleAttrib< Vector_t > & fld,
                const ParticleAttrib< Vector_t > & R,
                const size_t & comp,
                size_t from,
                size_t to)
    {
        const Vector_t & dx = fp.getSpacing();
        const Vector_t oneoverdx = Vector_t(1.0) / dx;
        const double onehalf = 0.5;
        size_t limitedComp = comp % 2;

        Vector_t add(0.0);
        add(limitedComp) = onehalf;
        to = std::min(to, fld.size());
        from = std::min(from, to);
        for (size_t i = from; i < to; ++ i) {
            Vector_t ppos = (R[i] - fp.getOrigin()) * oneoverdx;
            Vektor<size_t, 2U> lgp(std::floor(ppos(0) + add(0)),
                                   std::floor(ppos(1) + add(1)));
            lgp(limitedComp) -= 1;

            Vector_t gpos(lgp(0), lgp(1));
            gpos(limitedComp) += onehalf;

            Vector_t dpos = (ppos - gpos);
            const size_t & I = lgp(0);
            const size_t & J = lgp(1);

            fld[i](limitedComp) = ((1 - dpos(1)) * ((1 - dpos(0)) * fp.getRelative(I,     J) +
                                                    dpos(0)       * fp.getRelative(I + 1, J)) +
                                   dpos(1)       * ((1 - dpos(0)) * fp.getRelative(I,     J + 1) +
                                                    dpos(0)       * fp.getRelative(I + 1, J + 1)));
        }
    }

    static
    void gather(const FieldPatch<double> & fp,
                ParticleAttrib< double > & fld,
                const ParticleAttrib< Vector_t > & R,
                size_t from,
                size_t to)
    {
        static const double onehalf = 0.5;
        const Vector_t & dx = fp.getSpacing();
        const Vector_t oneoverdx = Vector_t(1.0) / dx;
        to = std::min(to, fld.size());
        from = std::min(from, to);
        for (size_t i = from; i < to; ++ i) {
            Vector_t ppos = (R[i] - fp.getOrigin()) * oneoverdx;
            Vektor<size_t, 2U> lgp(std::floor(ppos(0) - onehalf),
                                   std::floor(ppos(1) - onehalf));
            Vector_t gpos(lgp(0) + onehalf,
                          lgp(1) + onehalf);

            Vector_t dpos = (ppos - gpos);
            const size_t & I = lgp(0);
            const size_t & J = lgp(1);

            fld[i] = ((1 - dpos(1)) * ((1 - dpos(0)) * fp.getRelative(I,     J) +
                                       dpos(0)       * fp.getRelative(I + 1, J)) +
                      dpos(1)       * ((1 - dpos(0)) * fp.getRelative(I,     J + 1) +
                                       dpos(0)       * fp.getRelative(I + 1, J + 1)));
        }
    }

    static
    void scatter(FieldPatch<double> & fp,
                 const ParticleAttrib< Vector_t > & fld,
                 const ParticleAttrib< Vector_t > & R,
                 const size_t & comp)
    {
        const Vector_t & dx = fp.getSpacing();
        const Vector_t oneoverdx = Vector_t(1.0) / dx;
        const double onehalf = 0.5;
        size_t limitedComp = comp % 2;

        Vector_t add(0.0);
        add(limitedComp) = onehalf;
        for (size_t i = 0; i < fld.size(); ++ i) {
            Vector_t ppos = (R[i] - fp.getOrigin()) * oneoverdx;
            Vektor<size_t, 2U> lgp(std::floor(ppos(0) + add(0)),
                                   std::floor(ppos(1) + add(1)));
            lgp(limitedComp) -= 1;

            Vector_t gpos(lgp(0), lgp(1));
            gpos(limitedComp) += onehalf;

            Vector_t dpos = (ppos - gpos);
            const size_t & I = lgp(0);
            const size_t & J = lgp(1);

            fp.getRelative(I,     J)     += (1 - dpos(0)) * (1 - dpos(1)) * fld[i](limitedComp);
            fp.getRelative(I + 1, J)     += dpos(0)       * (1 - dpos(1)) * fld[i](limitedComp);
            fp.getRelative(I,     J + 1) += (1 - dpos(0)) * dpos(1)       * fld[i](limitedComp);
            fp.getRelative(I + 1, J + 1) += dpos(0)       * dpos(1)       * fld[i](limitedComp);
        }
    }

    static
    void scatter(FieldPatch<double> & fp,
                 const ParticleAttrib< double > & fld,
                 const ParticleAttrib< Vector_t > & R)
    {
        static const double onehalf = 0.5;
        const Vector_t & dx = fp.getSpacing();
        const Vector_t oneoverdx = Vector_t(1.0) / dx;
        for (size_t i = 0; i < fld.size(); ++ i) {
            Vector_t ppos = (R[i] - fp.getOrigin()) * oneoverdx;
            Vektor<size_t, 2U> lgp(std::floor(ppos(0) - onehalf),
                                   std::floor(ppos(1) - onehalf));
            Vector_t gpos(lgp(0) + onehalf,
                          lgp(1) + onehalf);

            Vector_t dpos = (ppos - gpos);
            const size_t & I = lgp(0);
            const size_t & J = lgp(1);

            fp.getRelative(I,     J) +=     (1 - dpos(0)) * (1 - dpos(1)) * fld[i];
            fp.getRelative(I + 1, J) +=     dpos(0)       * (1 - dpos(1)) * fld[i];
            fp.getRelative(I,     J + 1) += (1 - dpos(0)) * dpos(1)       * fld[i];
            fp.getRelative(I + 1, J + 1) += dpos(0)       * dpos(1)       * fld[i];

        }
    }

};

template<>
class ScatterGatherImpl<IntS1, double> {
public:
    static
    void gather(const FieldPatch<double> & fp,
                ParticleAttrib< Vector_t > & fld,
                const ParticleAttrib< Vector_t > & R,
                const size_t & comp,
                size_t from,
                size_t to)
    {
        static const double onehalf = 0.5;
        static const double threefourth = 0.75;
        const Vector_t & dx = fp.getSpacing();
        const Vector_t oneoverdx = Vector_t(1.0) / dx;
        const size_t limitedComp = comp % 2;

        Vektor<double, 2U> add(onehalf);
        add(limitedComp) = 0.0;

        to = std::min(to, fld.size());
        from = std::min(from, to);
        for (size_t i = from; i < to; ++ i) {
            const Vector_t ppos = (R[i] - fp.getOrigin()) * oneoverdx;
            const Vektor<size_t, 2U> lgp(std::floor(ppos(0) + add(0)),
                                         std::floor(ppos(1) + add(1)));

            Vector_t gpos(lgp(0), lgp(1));
            gpos(limitedComp) += onehalf;

            const Vector_t dpos = (ppos - gpos);

            // const double ax = onehalf * (onehalf - dpos(0))*(onehalf - dpos(0));
            const double bx = threefourth - dpos(0)*dpos(0);
            const double cx = onehalf * (onehalf + dpos(0))*(onehalf + dpos(0));
            const double ax = 1.0 - (bx + cx);

            // const double ay = onehalf * (onehalf - dpos(1))*(onehalf - dpos(1));
            const double by = threefourth - dpos(1)*dpos(1);
            const double cy = onehalf * (onehalf + dpos(1))*(onehalf + dpos(1));
            const double ay = 1.0 - (by + cy);

            const size_t & I = lgp(0);
            const size_t & J = lgp(1);

            fld[i](limitedComp) = (ay * (ax * fp.getRelative(I - 1, J - 1) +
                                         bx * fp.getRelative(I    , J - 1) +
                                         cx * fp.getRelative(I + 1, J - 1)) +
                                   by * (ax * fp.getRelative(I - 1, J    ) +
                                         bx * fp.getRelative(I    , J    ) +
                                         cx * fp.getRelative(I + 1, J    )) +
                                   cy * (ax * fp.getRelative(I - 1, J + 1) +
                                         bx * fp.getRelative(I    , J + 1) +
                                         cx * fp.getRelative(I + 1, J + 1)));
        }
    }

    static
    void gather(const FieldPatch<double> & fp,
                ParticleAttrib< double > & fld,
                const ParticleAttrib< Vector_t > & R,
                size_t from,
                size_t to)
    {
        static const double onehalf = 0.5;
        static const double threefourth = 0.75;
        const Vector_t & dx = fp.getSpacing();
        const Vector_t oneoverdx = Vector_t(1.0) / dx;

        to = std::min(to, fld.size());
        from = std::min(from, to);
        for (size_t i = from; i < to; ++ i) {
            const Vector_t ppos = (R[i] - fp.getOrigin()) * oneoverdx;
            const Vektor<size_t, 2U> lgp(std::floor(ppos(0)),
                                         std::floor(ppos(1)));
            const Vector_t gpos(lgp(0) + onehalf,
                                lgp(1) + onehalf);

            const Vector_t dpos = (ppos - gpos);

            // const double ax = onehalf * (onehalf - dpos(0))*(onehalf - dpos(0));
            const double bx = threefourth - dpos(0)*dpos(0);
            const double cx = onehalf * (onehalf + dpos(0))*(onehalf + dpos(0));
            const double ax = 1.0 - (bx + cx);

            // const double ay = onehalf * (onehalf - dpos(1))*(onehalf - dpos(1));
            const double by = threefourth - dpos(1)*dpos(1);
            const double cy = onehalf * (onehalf + dpos(1))*(onehalf + dpos(1));
            const double ay = 1.0 - (by + cy);

            const size_t & I = lgp(0);
            const size_t & J = lgp(1);

            fld[i] = (ay * (ax * fp.getRelative(I - 1, J - 1) +
                            bx * fp.getRelative(I    , J - 1) +
                            cx * fp.getRelative(I + 1, J - 1)) +
                      by * (ax * fp.getRelative(I - 1, J    ) +
                            bx * fp.getRelative(I    , J    ) +
                            cx * fp.getRelative(I + 1, J    )) +
                      cy * (ax * fp.getRelative(I - 1, J + 1) +
                            bx * fp.getRelative(I    , J + 1) +
                            cx * fp.getRelative(I + 1, J + 1)));
        }
    }

    static
    void scatter(FieldPatch<double> & fp,
                 const ParticleAttrib< Vector_t > & fld,
                 const ParticleAttrib< Vector_t > & R,
                 const size_t & comp)
    {
        static const double onehalf = 0.5;
        static const double threefourth = 0.75;
        const Vector_t & dx = fp.getSpacing();
        const Vector_t oneoverdx = Vector_t(1.0) / dx;
        const size_t limitedComp = comp % 2;

        Vektor<double, 2U> add(onehalf);
        add(limitedComp) = 0.0;

        for (size_t i = 0; i < fld.size(); ++ i) {
            const Vector_t ppos = (R[i] - fp.getOrigin()) * oneoverdx;
            const Vektor<size_t, 2U> lgp(std::floor(ppos(0) + add(0)),
                                         std::floor(ppos(1) + add(1)));

            Vector_t gpos(lgp(0), lgp(1));
            gpos(limitedComp) += onehalf;

            const Vector_t dpos = (ppos - gpos);

            // const double ax = onehalf * (onehalf - dpos(0))*(onehalf - dpos(0));
            const double bx = threefourth - dpos(0)*dpos(0);
            const double cx = onehalf * (onehalf + dpos(0))*(onehalf + dpos(0));
            const double ax = 1.0 - (bx + cx);

            // const double ay = onehalf * (onehalf - dpos(1))*(onehalf - dpos(1));
            const double by = threefourth - dpos(1)*dpos(1);
            const double cy = onehalf * (onehalf + dpos(1))*(onehalf + dpos(1));
            const double ay = 1.0 - (by + cy);

            const size_t & I = lgp(0);
            const size_t & J = lgp(1);

            fp.getRelative(I - 1, J - 1) += ax * ay * fld[i](limitedComp);
            fp.getRelative(I    , J - 1) += bx * ay * fld[i](limitedComp);
            fp.getRelative(I + 1, J - 1) += cx * ay * fld[i](limitedComp);
            fp.getRelative(I - 1, J    ) += ax * by * fld[i](limitedComp);
            fp.getRelative(I    , J    ) += bx * by * fld[i](limitedComp);
            fp.getRelative(I + 1, J    ) += cx * by * fld[i](limitedComp);
            fp.getRelative(I - 1, J + 1) += ax * cy * fld[i](limitedComp);
            fp.getRelative(I    , J + 1) += bx * cy * fld[i](limitedComp);
            fp.getRelative(I + 1, J + 1) += cx * cy * fld[i](limitedComp);
        }
    }

    static
    void scatter(FieldPatch<double> & fp,
                 const ParticleAttrib< double > & fld,
                 const ParticleAttrib< Vector_t > & R)
    {
        static const double onehalf = 0.5;
        static const double threefourth = 0.75;
        const Vector_t & dx = fp.getSpacing();
        const Vector_t oneoverdx = Vector_t(1.0) / dx;

        for (size_t i = 0; i < fld.size(); ++ i) {
            const Vector_t ppos = (R[i] - fp.getOrigin()) * oneoverdx;
            const Vektor<size_t, 2U> lgp(std::floor(ppos(0)),
                                         std::floor(ppos(1)));
            const Vector_t gpos(lgp(0) + onehalf,
                                lgp(1) + onehalf);

            const Vector_t dpos = (ppos - gpos);

            const double ax = onehalf * (onehalf - dpos(0))*(onehalf - dpos(0));
            const double bx = threefourth - dpos(0)*dpos(0);
            const double cx = onehalf * (onehalf + dpos(0))*(onehalf + dpos(0));
            const double ay = onehalf * (onehalf - dpos(1))*(onehalf - dpos(1));
            const double by = threefourth - dpos(1)*dpos(1);
            const double cy = onehalf * (onehalf + dpos(1))*(onehalf + dpos(1));

            const size_t & I = lgp(0);
            const size_t & J = lgp(1);

            fp.getRelative(I - 1, J - 1) += ax * ay * fld[i];
            fp.getRelative(I    , J - 1) += bx * ay * fld[i];
            fp.getRelative(I + 1, J - 1) += cx * ay * fld[i];
            fp.getRelative(I - 1, J    ) += ax * by * fld[i];
            fp.getRelative(I    , J    ) += bx * by * fld[i];
            fp.getRelative(I + 1, J    ) += cx * by * fld[i];
            fp.getRelative(I - 1, J + 1) += ax * cy * fld[i];
            fp.getRelative(I    , J + 1) += bx * cy * fld[i];
            fp.getRelative(I + 1, J + 1) += cx * cy * fld[i];
        }
    }


};

template<>
class ScatterGatherImpl<IntS2, double> {
private:
    // static
    // double powthree(double x) {
    //     return x*x*x;
    // }

public:
    static
    void gather(const FieldPatch<double> & fp,
                ParticleAttrib< Vector_t > & fld,
                const ParticleAttrib< Vector_t > & R,
                const size_t & comp,
                size_t from,
                size_t to)
    {
        static const double onehalf = 0.5;
        static const double twothirds = 2.0 / 3.0;
        static const double onesixth = 1.0 / 6.0;
        const Vector_t & dx = fp.getSpacing();
        const Vector_t oneoverdx = Vector_t(1.0) / dx;
        const size_t limitedComp = comp % 2;

        Vector_t add(0.0);
        add(limitedComp) = onehalf;

        to = std::min(to, fld.size());
        from = std::min(from, to);
        for (size_t i = from; i < to; ++ i) {
            const Vector_t ppos = (R[i] - fp.getOrigin()) * oneoverdx;
            Vektor<size_t, 2U> lgp(std::floor(ppos(0) + add(0)),
                                   std::floor(ppos(1) + add(1)));
            lgp(limitedComp) -= 1;

            Vector_t gpos(lgp(0), lgp(1));
            gpos(limitedComp) += onehalf;

            const Vector_t dpos = (ppos - gpos);
            const Vector_t dpossqr = dpos * dpos;
            // const double ax = onesixth * powthree(1-dpos(0));
            const double bx = twothirds - dpossqr(0) * (1 -  onehalf * dpos(0));
            const double cx = twothirds - onehalf * (1 - dpos(0))*(1 - dpossqr(0));
            const double ex = onesixth * dpossqr(0) * dpos(0);
            const double ax = 1.0 - (bx + cx + ex);

            // const double ay = onesixth * powthree(1-dpos(1));
            const double by = twothirds - dpossqr(1) * (1 - onehalf * dpos(1));
            const double cy = twothirds - onehalf * (1 - dpos(1))*(1 - dpossqr(1));
            const double ey = onesixth * dpossqr(1) * dpos(1);
            const double ay = 1.0 - (by + cy + ey);

            const size_t & I = lgp(0);
            const size_t & J = lgp(1);

            fld[i](limitedComp) = (ay * (ax * fp.getRelative(I - 1, J - 1) +
                                         bx * fp.getRelative(I    , J - 1) +
                                         cx * fp.getRelative(I + 1, J - 1) +
                                         ex * fp.getRelative(I + 2, J - 1)) +
                                   by * (ax * fp.getRelative(I - 1, J    ) +
                                         bx * fp.getRelative(I    , J    ) +
                                         cx * fp.getRelative(I + 1, J    ) +
                                         ex * fp.getRelative(I + 2, J    )) +
                                   cy * (ax * fp.getRelative(I - 1, J + 1) +
                                         bx * fp.getRelative(I    , J + 1) +
                                         cx * fp.getRelative(I + 1, J + 1) +
                                         ex * fp.getRelative(I + 2, J + 1)) +
                                   ey * (ax * fp.getRelative(I - 1, J + 2) +
                                         bx * fp.getRelative(I    , J + 2) +
                                         cx * fp.getRelative(I + 1, J + 2) +
                                         ex * fp.getRelative(I + 1, J + 2)));
        }
    }

    static
    void gather(const FieldPatch<double> & fp,
                ParticleAttrib< double > & fld,
                const ParticleAttrib< Vector_t > & R,
                size_t from,
                size_t to)
    {
        static const double onehalf = 0.5;
        static const double twothirds = 2.0 / 3.0;
        static const double onesixth = 1.0 / 6.0;
        const Vector_t & dx = fp.getSpacing();
        const Vector_t oneoverdx = Vector_t(1.0) / dx;

        to = std::min(to, fld.size());
        from = std::min(from, to);
        for (size_t i = from; i < to; ++ i) {
            const Vector_t ppos = (R[i] - fp.getOrigin()) * oneoverdx;
            const Vektor<size_t, 2U> lgp(std::floor(ppos(0) - onehalf),
                                         std::floor(ppos(1) - onehalf));
            const Vector_t gpos(lgp(0) + onehalf,
                                lgp(1) + onehalf);

            const Vector_t dpos = (ppos - gpos);
            const Vector_t dpossqr = dpos*dpos;
            // const double ax = onesixth * powthree(1-dpos(0));
            const double bx = twothirds - dpossqr(0) * (1 -  onehalf * dpos(0));
            const double cx = twothirds - onehalf * (1 - dpos(0))*(1 - dpossqr(0));
            const double ex = onesixth * dpossqr(0) * dpos(0);
            const double ax = 1.0 - (bx + cx + ex);

            // const double ay = onesixth * powthree(1-dpos(1));
            const double by = twothirds - dpossqr(1) * (1 - onehalf * dpos(1));
            const double cy = twothirds - onehalf * (1 - dpos(1))*(1 - dpossqr(1));
            const double ey = onesixth * dpossqr(1) * dpos(1);
            const double ay = 1.0 - (by + cy + ey);

            const size_t & I = lgp(0);
            const size_t & J = lgp(1);

            fld[i] = (ay * (ax * fp.getRelative(I - 1, J - 1) +
                            bx * fp.getRelative(I    , J - 1) +
                            cx * fp.getRelative(I + 1, J - 1) +
                            ex * fp.getRelative(I + 2, J - 1)) +
                      by * (ax * fp.getRelative(I - 1, J    ) +
                            bx * fp.getRelative(I    , J    ) +
                            cx * fp.getRelative(I + 1, J    ) +
                            ex * fp.getRelative(I + 2, J    )) +
                      cy * (ax * fp.getRelative(I - 1, J + 1) +
                            bx * fp.getRelative(I    , J + 1) +
                            cx * fp.getRelative(I + 1, J + 1) +
                            ex * fp.getRelative(I + 2, J + 1)) +
                      ey * (ax * fp.getRelative(I - 1, J + 2) +
                            bx * fp.getRelative(I    , J + 2) +
                            cx * fp.getRelative(I + 1, J + 2) +
                            ex * fp.getRelative(I + 1, J + 2)));
        }
    }

    static
    void scatter(FieldPatch<double> & fp,
                 const ParticleAttrib< Vector_t > & fld,
                 const ParticleAttrib< Vector_t > & R,
                 const size_t & comp)
    {
        static const double onehalf = 0.5;
        static const double twothirds = 2.0 / 3.0;
        static const double onesixth = 1.0 / 6.0;
        const Vector_t & dx = fp.getSpacing();
        const Vector_t oneoverdx = Vector_t(1.0) / dx;
        const size_t limitedComp = comp % 2;

        Vector_t add(0.0);
        add(limitedComp) = onehalf;
        for (size_t i = 0; i < fld.size(); ++ i) {
            const Vector_t ppos = (R[i] - fp.getOrigin()) * oneoverdx;
            Vektor<size_t, 2U> lgp(std::floor(ppos(0) + add(0)),
                                   std::floor(ppos(1) + add(1)));
            lgp(limitedComp) -= 1;

            Vector_t gpos(lgp(0), lgp(1));
            gpos(limitedComp) += onehalf;

            const Vector_t dpos = (ppos - gpos);
            const Vector_t dpossqr = dpos * dpos;
            // const double ax = onesixth * powthree(1-dpos(0));
            const double bx = twothirds - dpossqr(0) * (1 -  onehalf * dpos(0));
            const double cx = twothirds - onehalf * (1 - dpos(0))*(1 - dpossqr(0));
            const double ex = onesixth * dpossqr(0) * dpos(0);
            const double ax = 1.0 - (bx + cx + ex);

            // const double ay = onesixth * powthree(1-dpos(1));
            const double by = twothirds - dpossqr(1) * (1 - onehalf * dpos(1));
            const double cy = twothirds - onehalf * (1 - dpos(1))*(1 - dpossqr(1));
            const double ey = onesixth * dpossqr(1) * dpos(1);
            const double ay = 1.0 - (by + cy + ey);

            const size_t & I = lgp(0);
            const size_t & J = lgp(1);

            fp.getRelative(I - 1, J - 1) += ax * ay * fld[i](limitedComp);
            fp.getRelative(I    , J - 1) += bx * ay * fld[i](limitedComp);
            fp.getRelative(I + 1, J - 1) += cx * ay * fld[i](limitedComp);
            fp.getRelative(I + 2, J - 1) += ex * ay * fld[i](limitedComp);
            fp.getRelative(I - 1, J    ) += ax * by * fld[i](limitedComp);
            fp.getRelative(I    , J    ) += bx * by * fld[i](limitedComp);
            fp.getRelative(I + 1, J    ) += cx * by * fld[i](limitedComp);
            fp.getRelative(I + 2, J    ) += ex * by * fld[i](limitedComp);
            fp.getRelative(I - 1, J + 1) += ax * cy * fld[i](limitedComp);
            fp.getRelative(I    , J + 1) += bx * cy * fld[i](limitedComp);
            fp.getRelative(I + 1, J + 1) += cx * cy * fld[i](limitedComp);
            fp.getRelative(I + 2, J + 1) += ex * cy * fld[i](limitedComp);
            fp.getRelative(I - 1, J + 2) += ax * ey * fld[i](limitedComp);
            fp.getRelative(I    , J + 2) += bx * ey * fld[i](limitedComp);
            fp.getRelative(I + 1, J + 2) += cx * ey * fld[i](limitedComp);
            fp.getRelative(I + 2, J + 2) += ex * ey * fld[i](limitedComp);
        }
    }

    static
    void scatter(FieldPatch<double> & fp,
                 const ParticleAttrib< double > & fld,
                 const ParticleAttrib< Vector_t > & R)
    {
        static const double onehalf = 0.5;
        static const double twothirds = 2.0 / 3.0;
        static const double onesixth = 1.0 / 6.0;
        const Vector_t & dx = fp.getSpacing();
        const Vector_t oneoverdx = Vector_t(1.0) / dx;

        for (size_t i = 0; i < fld.size(); ++ i) {
            const Vector_t ppos = (R[i] - fp.getOrigin()) * oneoverdx;
            const Vektor<size_t, 2U> lgp(std::floor(ppos(0) - onehalf),
                                         std::floor(ppos(1) - onehalf));
            const Vector_t gpos(lgp(0) + onehalf,
                                lgp(1) + onehalf);

            const Vector_t dpos = (ppos - gpos);
            const Vector_t dpossqr = dpos*dpos;
            // const double ax = onesixth * powthree(1-dpos(0));
            const double bx = twothirds - dpossqr(0) * (1 -  onehalf * dpos(0));
            const double cx = twothirds - onehalf * (1 - dpos(0))*(1 - dpossqr(0));
            const double ex = onesixth * dpossqr(0) * dpos(0);
            const double ax = 1.0 - (bx + cx + ex);

            // const double ay = onesixth * powthree(1-dpos(1));
            const double by = twothirds - dpossqr(1) * (1 - onehalf * dpos(1));
            const double cy = twothirds - onehalf * (1 - dpos(1))*(1 - dpossqr(1));
            const double ey = onesixth * dpossqr(1) * dpos(1);
            const double ay = 1.0 - (by + cy + ey);

            const size_t & I = lgp(0);
            const size_t & J = lgp(1);

            fp.getRelative(I - 1, J - 1) += ax * ay * fld[i];
            fp.getRelative(I    , J - 1) += bx * ay * fld[i];
            fp.getRelative(I + 1, J - 1) += cx * ay * fld[i];
            fp.getRelative(I + 2, J - 1) += ex * ay * fld[i];
            fp.getRelative(I - 1, J    ) += ax * by * fld[i];
            fp.getRelative(I    , J    ) += bx * by * fld[i];
            fp.getRelative(I + 1, J    ) += cx * by * fld[i];
            fp.getRelative(I + 2, J    ) += ex * by * fld[i];
            fp.getRelative(I - 1, J + 1) += ax * cy * fld[i];
            fp.getRelative(I    , J + 1) += bx * cy * fld[i];
            fp.getRelative(I + 1, J + 1) += cx * cy * fld[i];
            fp.getRelative(I + 2, J + 1) += ex * cy * fld[i];
            fp.getRelative(I - 1, J + 2) += ax * ey * fld[i];
            fp.getRelative(I    , J + 2) += bx * ey * fld[i];
            fp.getRelative(I + 1, J + 2) += cx * ey * fld[i];
            fp.getRelative(I + 2, J + 2) += ex * ey * fld[i];
        }
    }


};

#endif
