/***************************************************************************
                           SaveEMFields.cpp
                         -------------------
    begin                : Tue Jan 14 2014
    copyright            : (C) 2014 by Christof Kraus
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
#include "SaveEMFields.hh"
#include "../PartBunch.hh"

extern std::ofstream dbg;
#define DBGOUT dbg << "SaveEMFields.cpp: " << __LINE__ << "\t"

SaveEMFieldsCommand::SaveEMFieldsCommand():
    _iteration(0)
{
    _saveTimer = IpplTimings::getTimer("line save");
}

SaveEMFieldsCommand::~SaveEMFieldsCommand()
{

}

void SaveEMFieldsCommand::execute(const VField_Edge_t & EFD,
                                  const VField_Cell_t & HFD,
                                  const VField_Edge_t & JFD,
                                  const PartBunchState & state,
                                  Inform & msg)
{
    IpplTimings::startTimer(_saveTimer);
    NDIndex<DIM> globalBounds;
    NDIndex<DIM> localFDomain = EFD.getLayout().getLocalNDIndex();
    NDIndex<DIM> localLDomain;
    unsigned int numSamples;
    long startIdx;
    std::vector<double> lineField;
    Vector_t R = state.get_rmean();
    Vector_t P = state.get_pmean();
    Vector_t h(EFD.get_mesh().get_meshSpacing(0),
               EFD.get_mesh().get_meshSpacing(1));
    Vector_t origin = EFD.get_mesh().get_origin();
    Vector_t l, u, lc, uc;
    bool inside = true;

    NDIndex<DIM> minGC = EFD.get_mesh().getCellContaining(state.get_rmin());
    NDIndex<DIM> maxGC = EFD.get_mesh().getCellContaining(state.get_rmax());
    for (unsigned int d = 0; d < DIM; ++ d) {
        globalBounds[d] = Index(minGC[d].first(), maxGC[d].first());
    }

    // l(0) = globalBounds[0].first();
    // u(0) = globalBounds[0].last() + 1;
    l(0) = std::floor(1.5 * globalBounds[0].first() - 0.5 * globalBounds[0].last());
    u(0) = 2 * globalBounds[0].last() - globalBounds[0].first();

    l(1) = (R(1) + (origin(0) + l(0) * h(0) - R(0)) / P(0) * P(1) - origin(1)) / h(1);
    u(1) = (R(1) + (origin(0) + u(0) * h(0) - R(0)) / P(0) * P(1) - origin(1)) / h(1);

    startIdx = static_cast<long>(l(0));
    numSamples = static_cast<unsigned int>(std::floor(u(0) - l(0) + 0.5)) + 1;
    lineField.resize(4 * numSamples);
    std::fill(lineField.begin(), lineField.end(), 0.0);

    lc = l;
    uc = u;
    getPartInside(lc, uc, EFD);

    localLDomain[0] = Index(static_cast<int>(std::floor(lc(0) + 0.5)),
                            static_cast<int>(std::floor(uc(0) + 0.5)));
    localLDomain[1] = Index(static_cast<int>(std::floor(lc(1))),
                            static_cast<int>(std::floor(uc(1))));

    for (unsigned int d = 0; d < DIM; ++ d) {
        int lower = std::max(localFDomain[d].first(), localLDomain[d].first());
        int upper = std::min(localFDomain[d].last(), localLDomain[d].last());

        if (lower > upper) {
            inside = false;
            break;
        }

        localLDomain[d] = Index(lower, upper);
    }

    if (inside) {
        ParticleAttrib<Vector_t> Pos;
        ParticleAttrib<Vector_t> Ef;
        ParticleAttrib<Vector_t> Hf;
        ParticleAttrib<Vector_t> Jf;
        Pos.create(localLDomain[0].length());
        Ef.create(localLDomain[0].length());
        Hf.create(localLDomain[0].length());
        Jf.create(localLDomain[0].length());

        unsigned int ii = 0;
        for (int i = localLDomain[0].first(); i <= localLDomain[0].last(); ++ i, ++ ii) {
            Pos[ii](0) = origin(0) + i * h(0);
            Pos[ii](1) = R(1) + (Pos[ii](0) - R(0)) / P(0) * P(1);
        }

        Ef.gather(EFD, Pos, IntOp_t());
        Hf.gather(HFD, Pos, IntOp_t());
        Jf.gather(JFD, Pos, IntOp_t());

        double lenP = sqrt(dot(P,P));
        double cosP = P(0) / lenP;
        double sinP = P(1) / lenP;
        ii = localLDomain[0].first() - startIdx ;
        unsigned int iii = 0;
        for (int i = localLDomain[0].first(); i <= localLDomain[0].last(); ++ i, ++ ii, ++ iii) {
            lineField[4 * ii    ] =  cosP * Ef[iii](0) + sinP * Ef[iii](1);
            lineField[4 * ii + 1] = -sinP * Ef[iii](0) + cosP * Ef[iii](1);
            lineField[4 * ii + 2] = Hf[iii](0) + Hf[iii](1);
            lineField[4 * ii + 3] =  cosP * Jf[iii](0) + sinP * Jf[iii](1);
        }
    }

    if (Ippl::myNode() == 0) {
        std::vector<double> cumLineField(4 * numSamples, 0.0);
        DBGOUT << "using collective communication" << std::endl;
        MPI_Reduce(&(lineField[0]), &(cumLineField[0]), 4 * numSamples, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        std::stringstream fn;
        fn << "Data/FieldOnLine_"
           << std::setw(4) << std::setfill('0') << _iteration << ".dat";

        std::ofstream out(fn.str().c_str());
        long ii = startIdx;
        for (unsigned int i = 0; i < numSamples; ++ i, ++ ii) {
            Vector_t curPos;
            curPos(0) = origin(0) + ii * h(0) - R(0);
            curPos(1) = curPos(0) / P(0) * P(1);

            out << sqrt(dot(curPos, curPos)) * (curPos(0) > 0? 1: -1) << "\t"
                << cumLineField[4 * i    ] << "\t"
                << cumLineField[4 * i + 1] << "\t"
                << cumLineField[4 * i + 2] << "\t"
                << cumLineField[4 * i + 3] << "\t"
                << curPos(0) + R(0) << "\t"
                << curPos(1) + R(1) << "\n";
        }
        out.close();
    } else {
        DBGOUT << "using collective communication" << std::endl;
        MPI_Reduce(&(lineField[0]), &(lineField[0]), 4 * numSamples, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    ++ _iteration;
    IpplTimings::stopTimer(_saveTimer);
}

void SaveEMFieldsCommand::getPartInside(Vector_t & begin,
                                        Vector_t & end,
                                        const VField_Edge_t & EFD)
{
    NDIndex<DIM> localFDomain = EFD.getLayout().getLocalNDIndex();
    int lower = static_cast<int>(std::floor(begin(0) + 0.5));
    int upper = static_cast<int>(std::floor(end(0) + 0.5));
    int minL = upper, maxL = upper;
    Vector_t projLength = end - begin;

    for (int i = lower; i <= upper; ++ i) {
        double y = begin(1) + (i - lower) / projLength(0) * projLength(1);
        int yi = static_cast<int>(std::floor(y));
        NDIndex<DIM> elemL(Index(i,i), Index(yi,yi));
        if (localFDomain.contains(elemL)) {
            if (i < minL) minL = i;
            if (i > minL) maxL = i;
        }
    }

    Vector_t tmp(minL,
                 begin(1) + (minL - lower) / projLength(0) * projLength(1));
    end = Vector_t(maxL,
                   begin(1) + (maxL - lower) / projLength(0) * projLength(1));
    begin = tmp;

    return;
}
