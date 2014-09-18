/***************************************************************************
                              utils.cpp
                         -------------------
    begin                : Thu Jan 12 2012
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

#include "utils.hh"
#include <unistd.h>
#include <boost/make_shared.hpp>

extern std::ofstream dbg;
#define DBGOUT dbg << __FILE__ << ": " << __LINE__ << "\t"

FieldLayout_Cell_t * Utils::getCellFromVertCenteredLayout(UniformCartesian<DIM> & mesh,
                                                          FieldLayout<DIM> & centeredLayout)
{
    std::vector<NDIndex<DIM> > vectorLocalDomains;
    std::vector<int> vectorNodeIndices;

    getLocalDomains(centeredLayout, vectorLocalDomains);
    decreaseLocalDomainsLast(centeredLayout, vectorLocalDomains);
    getNodeIndices(centeredLayout, vectorNodeIndices);

    return new FieldLayout_Cell_t(mesh,
                                  &vectorLocalDomains[0],
                                  &vectorLocalDomains[0] + Ippl::getNodes(),
                                  &vectorNodeIndices[0],
                                  &vectorNodeIndices[0] + Ippl::getNodes());
}

FieldLayout_Cell_t * Utils::getCellFromEdgeCenteredLayout(UniformCartesian<DIM> & mesh,
                                                          FieldLayout<DIM> & centeredLayout)
{
    std::vector<NDIndex<DIM> > vectorLocalDomains;
    std::vector<int> vectorNodeIndices;

    getLocalDomains(centeredLayout, vectorLocalDomains);
    decreaseLocalDomainsLast(centeredLayout, vectorLocalDomains);
    getNodeIndices(centeredLayout, vectorNodeIndices);

    return new FieldLayout_Cell_t(mesh,
                                  &vectorLocalDomains[0],
                                  &vectorLocalDomains[0] + Ippl::getNodes(),
                                  &vectorNodeIndices[0],
                                  &vectorNodeIndices[0] + Ippl::getNodes());
}

FieldLayout_Vert_t * Utils::getVertFromCellCenteredLayout(UniformCartesian<DIM> & mesh,
                                                          FieldLayout<DIM> & centeredLayout)
{
    std::vector<NDIndex<DIM> > vectorLocalDomains;
    std::vector<int> vectorNodeIndices;

    getLocalDomains(centeredLayout, vectorLocalDomains);
    increaseLocalDomainsLast(centeredLayout, vectorLocalDomains);
    getNodeIndices(centeredLayout, vectorNodeIndices);

    return new FieldLayout_Vert_t(mesh,
                                  &vectorLocalDomains[0],
                                  &vectorLocalDomains[0] + Ippl::getNodes(),
                                  &vectorNodeIndices[0],
                                  &vectorNodeIndices[0] + Ippl::getNodes());
}

FieldLayout_Vert_t * Utils::getVertFromEdgeCenteredLayout(UniformCartesian<DIM> & mesh,
                                                          FieldLayout<DIM> & centeredLayout)
{
    std::vector<NDIndex<DIM> > vectorLocalDomains;
    std::vector<int> vectorNodeIndices;

    getLocalDomains(centeredLayout, vectorLocalDomains);
    getNodeIndices(centeredLayout, vectorNodeIndices);

    return new FieldLayout_Vert_t(mesh,
                                  &vectorLocalDomains[0],
                                  &vectorLocalDomains[0] + Ippl::getNodes(),
                                  &vectorNodeIndices[0],
                                  &vectorNodeIndices[0] + Ippl::getNodes());
}

FieldLayout_Edge_t * Utils::getEdgeFromVertCenteredLayout(UniformCartesian<DIM> & mesh,
                                                          FieldLayout<DIM> & centeredLayout)
{
    std::vector<NDIndex<DIM> > vectorLocalDomains;
    std::vector<int> vectorNodeIndices;

    getLocalDomains(centeredLayout, vectorLocalDomains);
    getNodeIndices(centeredLayout, vectorNodeIndices);

    return new FieldLayout_Edge_t(mesh,
                              &vectorLocalDomains[0],
                              &vectorLocalDomains[0] + Ippl::getNodes(),
                              &vectorNodeIndices[0],
                              &vectorNodeIndices[0] + Ippl::getNodes());
}

FieldLayout_Edge_t * Utils::getEdgeFromCellCenteredLayout(UniformCartesian<DIM> & mesh,
                                                          FieldLayout<DIM> & centeredLayout)
{
    std::vector<NDIndex<DIM> > vectorLocalDomains;
    std::vector<int> vectorNodeIndices;

    getLocalDomains(centeredLayout, vectorLocalDomains);
    increaseLocalDomainsLast(centeredLayout, vectorLocalDomains);
    getNodeIndices(centeredLayout, vectorNodeIndices);

    return new FieldLayout_Edge_t(mesh,
                                  &vectorLocalDomains[0],
                                  &vectorLocalDomains[0] + Ippl::getNodes(),
                                  &vectorNodeIndices[0],
                                  &vectorNodeIndices[0] + Ippl::getNodes());
}

SField_Vert_ptr Utils::getScalarVertField(const VField_Vert_t & vectorField)
{
    Mesh_t & mesh = vectorField.get_mesh();
    const GuardCellSizes<DIM> & gcs = vectorField.getGuardCellSizes();

    return boost::make_shared<SField_Vert_t>(mesh, vectorField.getLayout(), gcs);
}

SField_Vert_ptr Utils::getScalarVertField(const VField_Cell_t & vectorField)
{
    Mesh_t & mesh = vectorField.get_mesh();
    FieldLayout_Vert_t * VertCenteredLayout = getVertFromCellCenteredLayout(mesh, vectorField.getLayout());
    const GuardCellSizes<DIM> & gcs = vectorField.getGuardCellSizes();

    return boost::shared_ptr<SField_Vert_t>(new SField_Vert_t(mesh, *VertCenteredLayout, gcs),
                                            CompanionDeleter<SField_Vert_t, FieldLayout_Vert_t>(VertCenteredLayout));
}

SField_Vert_ptr Utils::getScalarVertField(const VField_Edge_t & vectorField)
{
    Mesh_t & mesh = vectorField.get_mesh();
    FieldLayout_Vert_t * VertCenteredLayout = getVertFromEdgeCenteredLayout(mesh, vectorField.getLayout());
    const GuardCellSizes<DIM> & gcs = vectorField.getGuardCellSizes();

    return boost::shared_ptr<SField_Vert_t>(new SField_Vert_t(mesh, *VertCenteredLayout, gcs),
                                            CompanionDeleter<SField_Vert_t, FieldLayout_Vert_t>(VertCenteredLayout));
}

SField_Cell_ptr Utils::getScalarCellField(const VField_Vert_t & vectorField)
{
    Mesh_t & mesh = vectorField.get_mesh();
    FieldLayout_Cell_t * CellCenteredLayout = getCellFromVertCenteredLayout(mesh, vectorField.getLayout());
    const GuardCellSizes<DIM> & gcs = vectorField.getGuardCellSizes();

    return boost::shared_ptr<SField_Cell_t>(new SField_Cell_t(mesh, *CellCenteredLayout, gcs),
                                            CompanionDeleter<SField_Cell_t, FieldLayout_Cell_t>(CellCenteredLayout));
}

SField_Cell_ptr Utils::getScalarCellField(const VField_Cell_t & vectorField)
{
    Mesh_t & mesh = vectorField.get_mesh();
    const GuardCellSizes<DIM> & gcs = vectorField.getGuardCellSizes();

    return boost::make_shared<SField_Cell_t>(mesh, vectorField.getLayout(), gcs);
}

SField_Cell_ptr Utils::getScalarCellField(const VField_Edge_t & vectorField)
{
    Mesh_t & mesh = vectorField.get_mesh();
    FieldLayout_Cell_t * CellCenteredLayout = getCellFromEdgeCenteredLayout(mesh, vectorField.getLayout());
    const GuardCellSizes<DIM> & gcs = vectorField.getGuardCellSizes();

    return boost::shared_ptr<SField_Cell_t>(new SField_Cell_t(mesh, vectorField.getLayout(), gcs),
                                            CompanionDeleter<SField_Cell_t, FieldLayout_Cell_t>(CellCenteredLayout));
}

int Utils::ipow(int base, unsigned int exp)
{
    int result = 1;
    while (exp)
        {
            if (exp & 1)
                result *= base;
            exp >>= 1;
            base *= base;
        }

    return result;
}

void Utils::process_mem_usage(double& vm_usage, double& resident_set)
{
    using std::ios_base;
    using std::ifstream;
    using std::string;

    vm_usage     = 0.0;
    resident_set = 0.0;

    // 'file' stat seems to give the most reliable results
    //
    ifstream stat_stream("/proc/self/stat",ios_base::in);

    // dummy vars for leading entries in stat that we don't care about
    //
    string pid, comm, state, ppid, pgrp, session, tty_nr;
    string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    string utime, stime, cutime, cstime, priority, nice;
    string O, itrealvalue, starttime;

    // the two fields we want
    //
    unsigned long vsize;
    long rss;

    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
                >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
                >> utime >> stime >> cutime >> cstime >> priority >> nice
                >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

    stat_stream.close();

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage     = vsize / 1024.0;
    resident_set = rss * page_size_kb;
}

#ifdef TIMINGSDEBUGOUTPUT
Timings::TimerRef Timings::_current = -1U;
boost::timer::cpu_timer Timings::_mainTimer;
std::map<Timings::TimerRef, std::string> Timings::_ref2Name;
#endif
