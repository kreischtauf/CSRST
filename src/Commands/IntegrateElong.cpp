/***************************************************************************
                          IntegrateElong.cpp
                         -------------------
    begin                : Thu Jan 30 2014
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
#include "IntegrateElong.hh"
#include "../PartBunch.hh"

extern std::ofstream dbg;
#define DBGOUT dbg << "IntegrateElong.cpp: " << __LINE__ << "\t"

IntegrateElongCommand::IntegrateElongCommand(const VField_Edge_t & mmEFD,
                                             const double & dt):
    _dt(dt)
{
    NDIndex<DIM> lDom = mmEFD.getLayout().getDomain();
    _numSamples = lDom[0].length() - 4;
    _integral.resize(_numSamples, Utils::KahanAccumulation());
    _saveTimer = IpplTimings::getTimer("line integrate");
}

IntegrateElongCommand::~IntegrateElongCommand()
{

}

void IntegrateElongCommand::execute(const VField_Edge_t & mmEFD,
                                    const PartBunchState & state,
                                    Inform & msg)
{
    IpplTimings::startTimer(_saveTimer);
    NDIndex<DIM> globalFDomain = mmEFD.getLayout().getDomain();
    NDIndex<DIM> localFDomain = mmEFD.getLayout().getLocalNDIndex();
    Vector_t R = state.get_rmean();
    Vector_t P = state.get_pmean();
    Vector_t h(mmEFD.get_mesh().get_meshSpacing(0),
               mmEFD.get_mesh().get_meshSpacing(1));
    Vector_t origin = mmEFD.get_mesh().get_origin();
    Vector_t l;

    double lenP = sqrt(dot(P,P));
    double cosP = P(0) / lenP;
    double sinP = P(1) / lenP;
    Tenzor<double, DIM> M(cosP, sinP, -sinP, cosP);

    std::vector<Vector_t> samplePos;
    int firstIdx = 0, lastIdx = 0;
    for (int i = globalFDomain[0].first() + 2; i <= globalFDomain[0].last() - 2; ++ i) {
        l = Vector_t(origin(0) + i * h(0) - R(0), 0.0);
        l = dot(M, l) + R;
        NDIndex<DIM> elem = mmEFD.get_mesh().getCellContaining(l);
        if (localFDomain.contains(elem)) {
            samplePos.push_back(l);
            lastIdx = i;
        }
    }
    firstIdx = lastIdx - samplePos.size() + 1;

    if (samplePos.size() > 0) {
        ParticleAttrib<Vector_t> Pos;
        ParticleAttrib<Vector_t> Ef;
        Pos.create(samplePos.size());
        Ef.create(samplePos.size());

        for (size_t i = 0; i < samplePos.size(); ++ i)
            Pos[i] = samplePos[i];

        Ef.gather(mmEFD, Pos, IntOp_t());

        size_t ii = firstIdx - (globalFDomain[0].first() + 2) ;
        for (size_t i = 0; i < samplePos.size(); ++ i, ++ ii) {
            _integral[ii] = Utils::KahanSum(_integral[ii], cosP * Ef[i](0) + sinP * Ef[i](1));
        }
    }


    IpplTimings::stopTimer(_saveTimer);
}

void IntegrateElongCommand::write(std::string fname,
                                  const VField_Edge_t & mmEFD,
                                  const PartBunchState & state) const
{
    std::vector<double> lineField(_numSamples, 0.0);
    for (size_t i = 0; i < _numSamples; ++ i) {
        lineField[i] = _integral[i].sum;
    }

    if (Ippl::myNode() == 0) {
        Vector_t R = state.get_rmean();
        NDIndex<DIM> globalFDomain = mmEFD.getLayout().getDomain();
        Vector_t h(mmEFD.get_mesh().get_meshSpacing(0),
                   mmEFD.get_mesh().get_meshSpacing(1));
        Vector_t origin = mmEFD.get_mesh().get_origin();

        std::vector<double> cumLineField(_numSamples, 0.0);

        DBGOUT << "using collective communication" << std::endl;
        MPI_Reduce(&(lineField[0]), &(cumLineField[0]), _numSamples, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        std::ofstream out(fname.c_str());
        int ii = globalFDomain[0].first() + 2;
        for (unsigned int i = 0; i < _numSamples; ++ i, ++ ii) {
            double curPos = origin(0) + ii * h(0) - R(0);

            out << curPos << "\t" <<  -_dt * cumLineField[i] << "\n";
        }
        out.close();

    } else {
        DBGOUT << "using collective communication" << std::endl;
        MPI_Reduce(&(lineField[0]), &(lineField[0]), _numSamples, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
}
