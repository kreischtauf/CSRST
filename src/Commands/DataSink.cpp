/***************************************************************************
                             DataSink.cpp
                         -------------------
    begin                : Fri Jan 10 2014
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

#include "DataSink.hh"
#include "../PartBunch.hh"
#include "../csrst_config.h"

DataSinkCommand::DataSinkCommand(const VField_Edge_t & EFD):
    _cec(EFD)
{
    if (Ippl::myNode() > 0) return;

    _sink.open("Data/maxwadis.stat", std::ios::out);
    _sink.precision(15);
    _sink.setf(std::ios::scientific, std::ios::floatfield);
    writeSDDSHeader();
}

DataSinkCommand::~DataSinkCommand()
{
    if (Ippl::myNode() == 0)
        _sink.close();
}

void DataSinkCommand::writeSDDSHeader()
{
    _sink << "SDDS1" << std::endl;
    _sink << "&description text=\"Statistics data \", contents=\"stat parameters\" &end" << std::endl;

    _sink << "&parameter name=processors, type=long, "
          << "description=\"Number of Processors\" &end" << std::endl;

    _sink << "&parameter name=revision, type=string, "
          << "description=\"svn revision of CSRST\" &end" << std::endl;;

    _sink << "&column name=t,      type=double, units=s, "
          << "description=\"1 Time\" &end" << std::endl;

    _sink << "&column name=energy, type=double, units=MeV, "
          << "description=\"2 Mean Energy\" &end" << std::endl;

    _sink << "&column name=rms_x,  type=double, units=m , "
          << "description=\"3 RMS Beamsize in x  \" &end" << std::endl;
    _sink << "&column name=rms_y,  type=double, units=m , "
          << "description=\"4 RMS Beamsize in y  \" &end" << std::endl;

    _sink << "&column name=rms_px, type=double, units=1 , "
          << "description=\"5 RMS Momenta in x  \" &end" << std::endl;
    _sink << "&column name=rms_py, type=double, units=1 , "
          << "description=\"6 RMS Momenta in y  \" &end" << std::endl;

    _sink << "&column name=emit_x, type=double, units=m , "
          << "description=\"7 Normalized Emittance x  \" &end" << std::endl;
    _sink << "&column name=emit_y, type=double, units=m , "
          << "description=\"8 Normalized Emittance y  \" &end" << std::endl;

    _sink << "&column name=mean_x, type=double, units=m , "
          << "description=\"9 Mean Beam Position in x  \" &end" << std::endl;
    _sink << "&column name=mean_y, type=double, units=m , "
          << "description=\"10 Mean Beam Position in y  \" &end" << std::endl;

    _sink << "&column name=mean_px, type=double, units=1 , "
          << "description=\"11 Mean Beam Momentum in x  \" &end" << std::endl;
    _sink << "&column name=mean_py, type=double, units=1 , "
          << "description=\"12 Mean Beam Momentum in y  \" &end" << std::endl;

    _sink << "&column name=Dy,     type=double, units=m , "
          << "description=\"13 Dispersion in y  \" &end" << std::endl;
    _sink << "&column name=DDy,    type=double, units=1 , "
          << "description=\"14 Derivative of dispersion in y  \" &end" << std::endl;

    _sink << "&column name=dE,     type=double, units=MeV , "
          << "description=\"15 energy spread of the beam  \" &end" << std::endl;

    _sink << "&column name=Efld,     type=double, units=nJ , "
          << "description=\"16 total field energy   \" &end" << std::endl;

    _sink << "&data mode=ascii &end" << std::endl;

    _sink << "Cores used " << Ippl::getNodes() << std::endl;
    _sink << "CSRST " << PACKAGE_VERSION << " (" << VERSION_PATCH << ")" << std::endl;
}

void DataSinkCommand::execute(const PartBunch & bunch,
                              const double & t,
                              const VField_Edge_t & EFD,
                              const VField_Cell_t & HFD,
                              const VField_Edge_t & mmEFD,
                              const VField_Cell_t & mmHFD,
                              const Vektor<int, DIM> & start,
                              PartBunchState & bunchState,
                              Inform & msg)
{
    /// Set width of write fields in output files.
    const unsigned int pwi = 25;

    /// Calculate beam statistics and gather load balance statistics.
    bunch.calcBeamParameters(bunchState);

    double energy = _cec.execute(EFD, HFD, mmEFD, mmHFD, start, t, msg);

    if (Ippl::myNode() > 0) return;

    _sink << std::setw(pwi) << t                                                       // 1
          << std::setw(pwi) << bunchState.get_meanEnergy()                             // 2

        // rms r
          << std::setw(pwi) << bunchState.get_rrms()(0)                                // 3
          << std::setw(pwi) << bunchState.get_rrms()(1)                                // 4

        // rms p
          << std::setw(pwi) << bunchState.get_prms()(0)                                // 5
          << std::setw(pwi) << bunchState.get_prms()(1)                                // 6

        // normalized emittance
          << std::setw(pwi) << bunchState.get_norm_emit()(0)                           // 7
          << std::setw(pwi) << bunchState.get_norm_emit()(1)                           // 8

        // mean r
          << std::setw(pwi) << bunchState.get_rmean()(0)                               // 9
          << std::setw(pwi) << bunchState.get_rmean()(1)                               // 10

        // mean p
          << std::setw(pwi) << bunchState.get_pmean()(0)                               // 11
          << std::setw(pwi) << bunchState.get_pmean()(1)                               // 12

          << std::setw(pwi) << bunchState.get_Dy()                                     // 13
          << std::setw(pwi) << bunchState.get_DDy()                                    // 14

        // dE energy spread
          << std::setw(pwi) << bunchState.get_dE()                                     // 15

          << std::setw(pwi) << energy * 1e9                                            // 16
          << std::endl;
}
