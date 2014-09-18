/***************************************************************************
                            CalcPwrSpec.hh
                         -------------------
    begin                : Wed Sep 7 2011
    copyright            : (C) 2011 by Christof Kraus
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

#ifndef CALCPWRSPEC_HH
#define CALCPWRSPEC_HH

#include "Ippl.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <set>

#include "../Physics.hh"
#include "../defs.hh"

#define GUARDCELL 1


/***************************************************************************
 * PwrSpec
 *
 *
 ***************************************************************************/

template <class T, unsigned int Dim>
class CalcPwrSpec
{
public:

    typedef std::vector<T>                                  Spectrum_t;
    typedef std::list<Spectrum_t>                           SpecList_t;
    typedef Field<dcomplex, Dim, Mesh_t, Center_t>          CxField_t;
    typedef Field<T, Dim, Mesh_t, Center_t>                 RxField_t;
    typedef FFT<CCTransform, Dim, T>                        FFT_t;

    // constructor and destructor
    CalcPwrSpec(VField_t & field, const int & comp = -1);

    ~CalcPwrSpec() {
        if (fft_m)
            delete fft_m;
        fft_m = 0;
    }

    void execute();

    void save(const std::string & fn);

    void saveASCII(const std::string & fn);

    void save2D(const std::string & fn);

    void save2DASCII(const std::string & fn);

private:

    void savePPM(const std::string & fn,
                 const unsigned int & res_x,
                 const unsigned int & res_y,
                 double * data);

    void saveASCII(const std::string & fn,
                   const unsigned int & res_x,
                   const unsigned int & res_y,
                   double * data);

    void writePixel(std::ostream & out,
                    const double & value);

    void writeValue(std::ostream & out,
                    const double & value);

    /// fortrans nint function
    inline T nint(const T x)
    {
        return ceil(x + 0.5) - (fmod(x*0.5 + 0.25, 1.0) != 0);
    }

    FFT_t *fft_m;

    VField_t & field_;

    const int component_;

    // mesh and layout objects for rho_m
    Mesh_t *mesh_m;
    FieldLayout_t *layout_m;

    /// global domain for the various fields
    NDIndex<Dim> gDomain_m;
    /// local domain for the various fields
    NDIndex<Dim> lDomain_m;

    /// global domain for the enlarged fields
    NDIndex<Dim> gDomainL_m;

    BConds<T,Dim,Mesh_t,Center_t> bc_m;
    BConds<T,Dim,Mesh_t,Center_t> zerobc_m;

    e_dim_tag dcomp_m[Dim];
    Vektor<int, Dim> nr_m;

    /// Fourier transformed density field
    CxField_t rho_m;

    /// power spectra kmax
    int kmax_m;

    /// 1D power spectra
    T *spectra1D_m;

    /// Nk power spectra
    int *Nk_m;

    // time series of power spectra
    SpecList_t spectra_;
};

template <class T, unsigned int Dim>
CalcPwrSpec<T, Dim>::CalcPwrSpec(VField_t & field, const int & comp):
    field_(field),
    component_(comp),
    layout_m(&field.getLayout()),
    mesh_m(&field.get_mesh())
{

    Inform msg ("FFT doInit");
    gDomain_m = layout_m->getDomain();       // global domain
    lDomain_m = layout_m->getLocalNDIndex(); // local domain

#if Dim == 3
    gDomainL_m = NDIndex<Dim>(Index(gDomain_m[0].length()+1),
                              Index(gDomain_m[1].length()+1),
                              Index(gDomain_m[2].length()+1));
#else
    gDomainL_m = NDIndex<Dim>(Index(gDomain_m[0].length()+1),
                              Index(gDomain_m[1].length()+1));
#endif
    msg << "GDomain " << gDomain_m << " GDomainL " << gDomainL_m << endl;


    for (int i=0; i < 2*Dim; ++i) {
        bc_m[i] = new ParallelPeriodicFace<T,Dim,Mesh_t,Center_t>(i);
        zerobc_m[i] = new ZeroFace<T,Dim,Mesh_t,Center_t>(i);
    }

    for(int d=0; d<Dim; d++) {
        dcomp_m[d] = layout_m->getRequestedDistribution(d);
        nr_m[d] = gDomain_m[d].length();
    }

    // create additional objects for the use in the FFT's
    rho_m.initialize(*mesh_m, *layout_m, GuardCellSizes<Dim>(1));

    // create the FFT object
    bool compressTemps = true;
    fft_m = new FFT_t(layout_m->getDomain(), compressTemps);
    fft_m->setDirectionName(+1, "forward");

    // power spectra calc
    kmax_m = nr_m[Dim - 1] / 2;
    Nk_m = (int *) malloc (kmax_m * sizeof(int));

    if (Ippl::myNode() > 0)
        spectra_.push_back(std::vector<T>(kmax_m + 1, T()));
}


template <class T, unsigned int Dim>
void CalcPwrSpec<T, Dim>::execute() {
    Inform msg("CalcPwrSpec::execute ");


    if (Ippl::myNode() == 0)
        spectra_.push_back(std::vector<T>(kmax_m + 1, T()));

    unsigned long volume = nr_m[0] * nr_m[1];
    Index I = lDomain_m[0], J = lDomain_m[1];

#if Dim == 3
    Index K = lDomain_m[2];

    switch (component_) {

    case 0:
        assign(rho_m[I][J][K], field_[I][J][K](0));
        break;
    case 1:
        assign(rho_m[I][J][K], field_[I][J][K](1));
        break;
    case 2:
        assign(rho_m[I][J][K], field_[I][J][K](2));
        break;
    default:
        assign(rho_m[I][J][K],
               sqrt(field_[I][J][K](0) * field_[I][J][K](0) +
                    field_[I][J][K](1) * field_[I][J][K](1)));
    }
    volume *= nr_m[2];
#else
    switch (component_) {

    case 0:
        assign(rho_m[I][J], field_[I][J](0));
        break;
    case 1:
        assign(rho_m[I][J], field_[I][J](1));
        break;
    default:
        assign(rho_m[I][J],
               sqrt(field_[I][J](0) * field_[I][J](0) +
                    field_[I][J](1) * field_[I][J](1)));
    }
#endif

    unsigned int klen;

    T rho_0 = real(sum(rho_m)) / volume;
    rho_m = (rho_m - rho_0) / rho_0;

    fft_m->transform("forward", rho_m);

    rho_m = rho_m * conj(rho_m);

    for (int i = 0; i < kmax_m + 1; ++ i) {
        Nk_m[i]=0;
        (spectra_.back())[i] = 0.0;
    }

    msg << "Sum psp=real" << sum(real(rho_m)) << " kmax= " << kmax_m << endl;

    NDIndex<Dim> loop;

    // This computes the 1-D power spectrum of a 3-D field (rho_m)
    // by binning values
#if Dim == 3
    for (int i=lDomain_m[0].first(); i<=lDomain_m[0].last() && i <= kmax_m; i++) {
        loop[0]=Index(i,i);
        for (int j=lDomain_m[1].first(); j<=lDomain_m[1].last() && j <= kmax_m; j++) {
            loop[1]=Index(j,j);
            for (int k=lDomain_m[2].first(); k<=lDomain_m[2].last() && k <= kmax_m; k++) {

                loop[2]=Index(k,k);
                int ii = i;
                int jj = j;
                int kk = k;
                //FIXME: this seems to have no change
                //if(i >= (gDomain_m[0].max()+1)/2)
                if(i >= gDomain_m[0].max()/2)
                    ii -= nr_m[0];
                if(j >= gDomain_m[1].max()/2)
                    jj -= nr_m[1];
                if(k >= gDomain_m[2].max()/2)
                    kk -= nr_m[2];

                klen=(int)nint(sqrt(ii*ii + jj*jj + kk*kk));
                if (klen <= kmax_m) {
                    (spectra_.back())[klen] += real(rho_m.localElement(loop));
                    Nk_m[klen]++;
                }
            }
        }
    }
#else
    for (int i=lDomain_m[0].first(); i<=lDomain_m[0].last() && i <= kmax_m; i++) {
        loop[0]=Index(i,i);
        for (int k=lDomain_m[1].first(); k<=lDomain_m[1].last() && k <= kmax_m; k++) {
            loop[1]=Index(k,k);

            int ii = i;
            int kk = k;
            //FIXME: this seems to have no change
            //if(i >= (gDomain_m[0].max()+1)/2)
            if(i >= gDomain_m[0].max()/2)
                ii -= nr_m[0];
            if(k >= gDomain_m[1].max()/2)
                kk -= nr_m[1];

            klen=(int)nint(sqrt(ii*ii + kk*kk));
            if (klen <= kmax_m) {
                (spectra_.back())[klen] += real(rho_m.localElement(loop));
                Nk_m[klen]++;
            }
        }
    }
#endif
    INFOMSG("Loops done" << endl);

    reduce(&(Nk_m[0]), &(Nk_m[0]) + kmax_m, &(Nk_m[0]), OpAddAssign());
    reduce(&((spectra_.back())[0]), &((spectra_.back())[0]) + kmax_m, &((spectra_.back())[0]) ,OpAddAssign());


    T scale = 1.0; // std::pow((T)(simData_m.rL/simData_m.ng_comp),(T)3.0);
    //     int sumNk = 0;
    for (int i = 0; i < kmax_m; ++ i) {
        //         sumNk += Nk_m[i];
        (spectra_.back())[i] /= 1.0 * (Nk_m[i] > 0? Nk_m[i]: 1);
        (spectra_.back())[i] *= scale;
    }
}


template <class T, unsigned int Dim>
void CalcPwrSpec<T, Dim>::save(const std::string & fn) {
    Inform msg("CalcPwrSpec::save ");

    if (Ippl::myNode() == 0) {
        const int res_x = spectra_.size();
        const int res_y = (spectra_.front()).size();
        double * values = new double[res_x * res_y];

        int i = 0;
        for (typename SpecList_t::iterator it = spectra_.begin();
             it != spectra_.end();
             ++ it, ++ i) {
            int j = res_y - 1;
            for (typename Spectrum_t::iterator it2 = (*it).begin();
                 it2 != (*it).end();
                 ++ it2, -- j) {
                values[j * res_x + i] = *it2;
            }
        }

        savePPM(fn, res_x, res_y, values);

        delete[] values;
    }
}


template <class T, unsigned int Dim>
void CalcPwrSpec<T, Dim>::saveASCII(const std::string & fn) {
    Inform msg("CalcPwrSpec::save ");

    if (Ippl::myNode() == 0) {
        const int res_x = spectra_.size();
        const int res_y = (spectra_.front()).size();
        double * values = new double[res_x * res_y];

        int i = 0;
        for (typename SpecList_t::iterator it = spectra_.begin();
             it != spectra_.end();
             ++ it, ++ i) {
            int j = res_y - 1;
            for (typename Spectrum_t::iterator it2 = (*it).begin();
                 it2 != (*it).end();
                 ++ it2, -- j) {
                values[j * res_x + i] = *it2;
            }
        }

        saveASCII(fn, res_x, res_y, values);

        delete[] values;
    }
}


template <class T, unsigned int Dim>
void CalcPwrSpec<T, Dim>::save2D(const std::string & fn) {
    Inform msg("CalcPwrSpec::save2D ");

    const int res_x = nr_m[0] / 2;
    const int res_y = nr_m[1] / 2;
    double * values = new double[res_x * res_y];

    for (int i = 0; i < res_x * res_y; ++ i)
        values[i] = 0.0;

    int XLowLimit = max(lDomain_m[0].first(), 0);
    int XHighLimit = min((lDomain_m[0].last() + 1) / 2, res_x);
    int YLowLimit = max(lDomain_m[1].first(), 0);
    int YHighLimit = min((lDomain_m[1].last() + 1) / 2, res_y);

    NDIndex<Dim> loop;
    int k = 0;
    for (int j = YHighLimit - 1; j >= YLowLimit; -- j) {
//     for (int j = YLowLimit; j < YHighLimit; ++ j) {
        loop[1] = Index(j,j);
//         k = j * res_x + XLowLimit;
        for (int i = XLowLimit; i < XHighLimit; ++ i) {
            loop[0] = Index(i,i);
            values[k ++] = real(rho_m.localElement(loop));
        }
    }

    reduce(values, values + res_x * res_y, values, OpAddAssign());

    if (Ippl::myNode() == 0)
        savePPM(fn, res_x, res_y, values);

    delete[] values;
}


template <class T, unsigned int Dim>
void CalcPwrSpec<T, Dim>::save2DASCII(const std::string & fn) {
    Inform msg("CalcPwrSpec::save2DASCII ");

    const int res_x = nr_m[0] / 2;
    const int res_y = nr_m[1] / 2;
    double * values = new double[res_x * res_y];

    for (int i = 0; i < res_x * res_y; ++ i)
        values[i] = 0.0;

    int XLowLimit = max(lDomain_m[0].first(), 0);
    int XHighLimit = min((lDomain_m[0].last() + 1) / 2, res_x);
    int YLowLimit = max(lDomain_m[1].first(), 0);
    int YHighLimit = min((lDomain_m[1].last() + 1) / 2, res_y);

    NDIndex<Dim> loop;
    int k = 0;
    for (int j = YHighLimit - 1; j >= YLowLimit; -- j) {
        loop[1] = Index(j,j);
        k = j * res_x + XLowLimit;
        for (int i = XLowLimit; i < XHighLimit; ++ i) {
            loop[0] = Index(i,i);
            values[k ++] = real(rho_m.localElement(loop));
        }
    }

    reduce(values, values + res_x * res_y, values, OpAddAssign());

    if (Ippl::myNode() == 0)
        saveASCII(fn, res_x, res_y, values);

    delete[] values;
}


template <class T, unsigned int Dim>
void CalcPwrSpec<T, Dim>::savePPM(const std::string & fn,
                                  const unsigned int & res_x,
                                  const unsigned int & res_y,
                                  double * data) {
    Inform msg("CalcPwrSpec::savePPM ");

    int k = 0;
    while (data[k] != data[k])
        k++;
    double maxValue = data[k];
    double minValue = maxValue;
    double sign = 1.;

    for (unsigned int i = 0; i < res_x * res_y; ++ i) {

        if (data[i] == data[i]) {
            if (data[i] > maxValue) {
                maxValue = data[i];
            }
            if (data[i] < minValue) {
                minValue = data[i];
            }
        }
    }

    if (std::abs(minValue) > std::abs(maxValue)) {
        sign = -1.;
        double temp = maxValue;
        maxValue = -minValue;
        minValue = -temp;
    }

    msg << "min: " << 0.0 << ", max: " << maxValue - minValue<< endl;

    // Die Datei wird binaer geschrieben
    std::ofstream outfile(fn.c_str(), std::ios::binary | std::ios::out);

    // Definiert den Header der PPM-Datei
    outfile << "P6" << endl;
    outfile << res_x + 8 << endl;
    outfile << res_y << endl;
    outfile << 255 << endl;


    for (unsigned int i = 0; i < res_x * res_y; ++ i) {
        data[i] = (data[i] == data[i] ? sign * data[i] - minValue: 0.0);
        data[i] = log(data[i] + 1) / log(maxValue - minValue + 1);
    }

    // Fuer jeden Pixel werden die Farben geschrieben
    k = 0;
    for (int j = 0; j < res_y; ++ j) {
        for (int i = 0; i < res_x; ++ i) {
           double & relValue = data[k++];
           writePixel(outfile, relValue);
        }
        for (int i = 0; i < (8-2)*3; ++ i) {
            outfile << (unsigned char)(215);
        }
        double relValue = log((res_y - j) * 1.0 / res_y + 1) / log(2);// - log(maxValue - minValue + 1);
        writePixel(outfile, relValue);
        writePixel(outfile, relValue);
    }
}

template <class T, unsigned int Dim>
void CalcPwrSpec<T, Dim>::saveASCII(const std::string & fn,
                                    const unsigned int & res_x,
                                    const unsigned int & res_y,
                                    double * data) {
    Inform msg("CalcPwrSpec::saveASCII ");


//     double maxValue = data[0];
//     double minValue = maxValue;

//     for (unsigned int i = 0; i < res_x * res_y; ++ i) {

//         if (data[i] == data[i]) {
//             if (data[i] > maxValue) {
//                 maxValue = data[i];
//             }
//             if (data[i] < minValue) {
//                 minValue = data[i];
//             }
//         }
//     }

//     msg << "min: " << 0.0 << ", max: " << maxValue - minValue<< endl;

    std::ofstream outfile(fn.c_str());


    for (unsigned int i = 0; i < res_x * res_y; ++ i) {
        data[i] = (data[i] == data[i] ? data[i] : 0.0);
//         data[i] = (data[i] == data[i] ? data[i] - minValue: 0.0);
//         data[i] = data[i] / (maxValue - minValue);
    }

    int k = 0;
    for (int j = 0; j < res_y; ++ j) {
        for (int i = 0; i < res_x; ++ i) {
           double relValue = data[k++];
           writeValue(outfile, relValue);
        }
        outfile << endl;
    }
}

template <class T, unsigned int Dim>
void CalcPwrSpec<T, Dim>::writePixel(std::ostream & out, const double & value) {
    switch ((int)(3 * value)) {
    case 0:
        out << (unsigned char)(value * 765)  // R
            << (unsigned char)(0)            // G
            << (unsigned char)(0);           // B
        break;
    case 1:
        out << (unsigned char)(255)
            << (unsigned char)((value - 0.3333) * 765)
            << (unsigned char)(0);
        break;
    default:
        out << (unsigned char)(255)
            << (unsigned char)(255)
            << (unsigned char)((value - 0.6666) * 765);
    }
}

template <class T, unsigned int Dim>
void CalcPwrSpec<T, Dim>::writeValue(std::ostream & out, const double & value) {
    out << std::setw(16) << std::setprecision(8) << value;
}

#endif
