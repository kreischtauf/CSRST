/***************************************************************************
                           PlotASCIIVTK.hh
                         -------------------
    begin                : Mon Oct 4 2010
    copyright            : (C) 2010 by Christof Kraus
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

#ifndef PLOTASCIIVTKCOMMAND_HH
#define PLOTASCIIVTKCOMMAND_HH

#include <vector>
#include "../defs.hh"
#include "../PartBunch.hh"

class PlotASCIIVTKCommand {
public:
    PlotASCIIVTKCommand(const VField_Edge_t & EFD,
                        const VField_Cell_t & HFD,
                        const VField_Edge_t & JFD,
                        const std::string & baseName = "c");

    void closeFile();

    virtual void execute(const double & t,
                         const PartBunch & bunch,
                         Inform & msg);

private:
    const VField_Edge_t & _EFD;
    const VField_Cell_t & _HFD;
    const VField_Edge_t & _JFD;

    std::string _baseName;

    const double _dx;
    const double _dy;

    unsigned int _iteration;

    std::vector<NDIndex<DIM> > _lDoms;
    std::vector<double> _timesteps;
    Timings::TimerRef _outputTimer;
};

inline
void PlotASCIIVTKCommand::closeFile()
{
    std::ofstream vtkout;
    short EndianTest_s;
    char *EndianTest = reinterpret_cast<char*>(&EndianTest_s);
    EndianTest[0] = 1;
    EndianTest[1] = 0;

    std::stringstream fname;
    fname << "Data/" << _baseName << "_collection.pvd";
    vtkout.open(fname.str().c_str());
    vtkout << "<?xml version=\"1.0\"?>\n"
           << indent_l0 << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"" << (EndianTest_s == 1? "LittleEndian": "BigEndian") << "\">\n"
           << indent_l1 << "<Collection>\n";
    for (unsigned int i = 0; i < _iteration; ++ i) {
        vtkout << indent_l2 << "<DataSet timestep=\"" << _timesteps[i] << "\" group=\"\" part=\"0\" file=\"c_" << std::setw(4) << std::setfill('0') << i << ".pvtr\"/>\n";
    }
    vtkout << indent_l1 << "</Collection>\n"
           << indent_l0 << "</VTKFile>\n"
           << std::endl;
    vtkout.close();

}

#endif
