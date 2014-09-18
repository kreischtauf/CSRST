/***************************************************************************
                          PlotBinaryVTK.cpp
                         -------------------
    begin                : Tue Jun 23 2009
    copyright            : (C) 2009 by Christof Kraus
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

#include "PlotBinaryVTK.hh"
#include "BinaryVTKFile.hh"
#include "../utils.hh"

#include <boost/regex.hpp>

#include <fstream>
#include <sstream>

extern std::ofstream dbg;
#define DBGOUT dbg << "PlotBinaryVTK.cpp: " << __LINE__ << "\t"

PlotBinaryVTKCommand::PlotBinaryVTKCommand(const VField_Edge_t & EFD,
                                           const VField_Cell_t & HFD,
                                           const VField_Edge_t & JFD,
                                           const std::string & baseName):
    _dx(EFD.get_mesh().get_meshSpacing(0)),
    _dy(EFD.get_mesh().get_meshSpacing(1)),
    _EFD(EFD),
    _HFD(HFD),
    _JFD(JFD),
    _baseName(baseName),
    _iteration(0)
{
    _outputTimer = Timings::getTimer("binary_vtk_output");

    FieldLayout<DIM> & FL = EFD.getLayout();
    Utils::getLocalDomains(FL, _lDoms);
    Utils::addGostCellToLocalDomains(FL, _lDoms);
    Utils::decreaseLocalDomainsLast(FL, _lDoms);

    parseCollectionFile();
}

void PlotBinaryVTKCommand::execute(const double & t,
                                   const PartBunch & bunch,
                                   Inform & msg)
{
    NDIndex<DIM> elem;
    const NDIndex<DIM> & localDomain = _lDoms[Ippl::myNode()];
    std::stringstream fname;

    _timesteps.push_back(t);

    Timings::startTimer(_outputTimer);

    fname << "Data/" << _baseName << "_";
    fname << std::setw(4) << std::setfill('0') << _iteration;

    BinaryVtkFile newStep;
    newStep.addVectorField(_EFD, "E-Field");
    newStep.addVectorField(_JFD, "J-Field");

    boost::shared_ptr<SField_Cell_t> HField = Utils::getScalarCellField(_HFD);

    for (int j = localDomain[1].first(); j <= localDomain[1].last(); ++ j) {
        elem[1] = Index(j,j);
        for (int i = localDomain[0].first() ; i <= localDomain[0].last(); ++ i) {
            elem[0] = Index(i,i);
            Vector_t hfd = _HFD.localElement(elem);
            HField->localElement(elem) = hfd(0) + hfd(1);
        }
    }
    newStep.addScalarField(*HField, "H-Field");

    newStep.writeFile(fname.str());

    Timings::stopTimer(_outputTimer);

    ++_iteration;
}

void PlotBinaryVTKCommand::parseCollectionFile()
{
    std::string fileName("Data/collection.pvd");
    std::ifstream is(fileName.c_str());
    // is.exceptions(std::ifstream::failbit | std::ifstream::badbit);

    std::string fullLine;

    std::string::const_iterator start;
    std::string::const_iterator end;
    boost::regex expression("<DataSet timestep=\"(.*?)\"");//\"(\\d*\\.\\d*)\"");
    boost::match_results<std::string::const_iterator> what;
    boost::match_flag_type flags = boost::match_default;

    while (is.good() && is.peek() != EOF) {
        getline(is, fullLine);
        start = fullLine.begin();
        end = fullLine.end();
        if (regex_search(start, end, what, expression, flags)) {
            std::stringstream ss;
            double time;

            ss << what[1];
            ss >> time;
            _timesteps.push_back(time);
            _iteration ++;
        }
    }
    is.close();
}
