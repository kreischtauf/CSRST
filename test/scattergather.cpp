#include <initializer_list>
#include "defs.hh"
#include "FieldPatch.hh"
#include "Utility/RNGBitReverse.h"

template<class PL>
class ChargedParticles : public ParticleBase<PL> {
public:
  ParticleAttrib<double>     Q;       // charge-to-mass ratio
  ParticleAttrib<double>     B;        // particle velocity
  typename PL::ParticlePos_t E;        // electric field at particle position
  ChargedParticles(PL* pl) : ParticleBase<PL>(pl) {
    // register the particle attributes
    this->addAttribute(Q);
    this->addAttribute(B);
    this->addAttribute(E);
  }
};

struct Vcheck {
    unsigned int i;
    unsigned int j;
    Vector_t value;
    Vcheck(const unsigned int & ii,
           const unsigned int & ji,
           const Vector_t & valuei):
        i(ii),
        j(ji),
        value(valuei)
    { }
};

struct Scheck {
    unsigned int i;
    unsigned int j;
    double value;
    Scheck(const unsigned int & ii,
           const unsigned int & ji,
           const double & valuei):
        i(ii),
        j(ji),
        value(valuei)
    { }
};

const int nx=20, ny=20;              // size of domain is nx X ny
const unsigned int totalP = 1;      // number of particles to create
const int nt = 0;                    // total number of timesteps

const double pi = acos(-1.0);
const double qmmax = 1.0;       // maximum value for particle q/m
const double dt = 1.0;          // size of timestep

int main(int argc, char *argv[]){
    TAU_PROFILE("main()", "int (int, char **)", TAU_DEFAULT);

    Ippl ippl(argc, argv);
    Inform testmsg(argv[0]);
    testmsg << "Particle test EDGE: Begin." << endl;

    bool errors = false;

    // create layout objects
    Index I1(nx+1), J1(ny+1);
    Mesh_t mymesh(I1,J1);
    // FieldLayout_Vert_t FL_vert(mymesh);
    // FieldLayout_Edge_t FL_edge(mymesh);
    // FieldLayout_Cell_t FL_cell(mymesh);


    FieldPatch<double> ex(NDIndex<2U>(I1, J1));
    FieldPatch<double> ey(NDIndex<2U>(I1, J1));
    FieldPatch<double> bz(NDIndex<2U>(I1, J1));

    ex.setOrigin(Vector_t(0.0));
    ey.setOrigin(Vector_t(0.0));
    bz.setOrigin(Vector_t(0.0));

    ex.setSpacing(Vector_t(1.0));
    ey.setSpacing(Vector_t(1.0));
    bz.setSpacing(Vector_t(1.0));

    // create an empty ChargedParticles object, setting it to use periodic BC's
    playout_t* PL = new playout_t();
    ChargedParticles<playout_t> P(PL);
    P.create(totalP);
    P.R[0] = Vector_t(3.3,3.2);

    // CHECK CIC ////////////////////////////////////////////////////////////
    std::vector<Vcheck> EFDcheck = {Vcheck(2,3,Vector_t(0.16,0)),
                                    Vcheck(2,4,Vector_t(0.04,0)),
                                    Vcheck(3,2,Vector_t(0,0.21)),
                                    Vcheck(3,3,Vector_t(0.64,0.49)),
                                    Vcheck(3,4,Vector_t(0.16,0)),
                                    Vcheck(4,2,Vector_t(0,0.09)),
                                    Vcheck(4,3,Vector_t(0,0.21))};
    std::vector<Scheck> BFDcheck = {Scheck(2,2,0.06),
                                    Scheck(2,3,0.14),
                                    Scheck(3,2,0.24),
                                    Scheck(3,3,0.56)};

    P.E[0] = Vector_t(1.0,1.0);
    P.B[0] = 1.0;

    ex.reset();
    ey.reset();
    bz.reset();

    ex.scatter(P.E, P.R, 0, IntCIC());
    ey.scatter(P.E, P.R, 1, IntCIC());
    bz.scatter(P.B, P.R, IntCIC());

    testmsg << "-- Testing IntCIC scatter on edge centered fields" << endl;
    for (const Vcheck c: EFDcheck) {
        if (std::abs(ex(c.i,c.j) - c.value(0)) > 1e-8 ||
            std::abs(ey(c.i,c.j) - c.value(1)) > 1e-8) {
            ERRORMSG("CIC scatter not working for edge centered fields" << endl);
            errors = true;
        }
    }

    testmsg << "-- Testing IntCIC scatter on cell centered fields" << endl;
    for (const Scheck c: BFDcheck) {
        if (std::abs(bz(c.i, c.j) - c.value) > 1e-8) {
            ERRORMSG("CIC scatter not working for cell centered fields" << endl);
            errors = true;
        }
    }

    P.E[0] = Vector_t(0.0);
    P.B[0] = 0.0;

    ex.reset();
    ex(2,3) = 1.0;
    ex(3,3) = 1.0;
    ex(2,4) = 1.0;
    ex(3,4) = 1.0;

    ey.reset();
    ey(3,2) = 1.0;
    ey(4,2) = 1.0;
    ey(3,3) = 1.0;
    ey(4,3) = 1.0;

    bz.reset();
    bz(2,2) = 1.0;
    bz(3,2) = 1.0;
    bz(2,3) = 1.0;
    bz(3,3) = 1.0;

    ex.gather(P.E, P.R, 0, IntCIC());
    ey.gather(P.E, P.R, 1, IntCIC());
    bz.gather(P.B, P.R, IntCIC());

    testmsg << "-- Testing IntCIC gather on edge centered fields" << endl;
    if (std::abs(P.E[0](0) - 1.0) > 1e-8 ||
        std::abs(P.E[0](1) - 1.0) > 1e-8) {
        ERRORMSG("CIC gather not working for edge centered fields" << endl);
        errors = true;
    }

    testmsg << "-- Testing IntCIC gather on cell centered fields" << endl;
    if (std::abs(P.B[0] - 1.0) > 1e-8) {
        ERRORMSG("CIC gather not working for cell centered fields" << endl);
        errors = true;
    }

    // CHECK S1 /////////////////////////////////////////////////////////////
    EFDcheck.assign({Vcheck(2,2,Vector_t(0.011025,0.0064)),
                     Vcheck(2,3,Vector_t(0.173950,0.0132)),
                     Vcheck(2,4,Vector_t(0.060025,0.0004)),
                     Vcheck(3,2,Vector_t(0.031950,0.2112)),
                     Vcheck(3,3,Vector_t(0.504100,0.4356)),
                     Vcheck(3,4,Vector_t(0.173950,0.0132)),
                     Vcheck(4,2,Vector_t(0.002025,0.1024)),
                     Vcheck(4,3,Vector_t(0.031950,0.2112)),
                     Vcheck(4,4,Vector_t(0.011025,0.0064))});

    BFDcheck.assign({Scheck(2,2,0.0784),
                     Scheck(2,3,0.1617),
                     Scheck(2,4,0.0049),
                     Scheck(3,2,0.2272),
                     Scheck(3,3,0.4686),
                     Scheck(3,4,0.0142),
                     Scheck(4,2,0.0144),
                     Scheck(4,3,0.0297),
                     Scheck(4,4,0.0009)});

    P.E[0] = Vector_t(1.0);
    P.B[0] = 1.0;

    ex.reset();
    ey.reset();
    bz.reset();

    ex.scatter(P.E, P.R, 0, IntS1());
    ey.scatter(P.E, P.R, 1, IntS1());
    bz.scatter(P.B, P.R, IntS1());

    testmsg << "-- Testing IntS1 scatter on edge centered fields" << endl;
    for (const Vcheck c: EFDcheck) {
        if (std::abs(ex(c.i,c.j) - c.value(0)) > 1e-8 ||
            std::abs(ey(c.i,c.j) - c.value(1)) > 1e-8) {
            ERRORMSG("S1 scatter interpolation not working for edge based fields" << endl);
            errors = true;
        }
    }


    testmsg << "-- Testing IntS1 scatter on cell centered fields" << endl;
    for (const Scheck c: BFDcheck) {
        if (std::abs(bz(c.i, c.j) - c.value) > 1e-8) {
            ERRORMSG("S1 scatter interpolation not working for cell centered fields" << endl);
            errors = true;
        }
    }

    P.E[0] = Vector_t(0.0);
    P.B[0] = 0.0;

    ex.reset();
    ex(2,2) = 1.0;
    ex(3,2) = 1.0;
    ex(4,2) = 1.0;
    ex(2,3) = 1.0;
    ex(3,3) = 1.0;
    ex(4,3) = 1.0;
    ex(2,4) = 1.0;
    ex(3,4) = 1.0;
    ex(4,4) = 1.0;

    ey.reset();
    ey(2,2) = 1.0;
    ey(3,2) = 1.0;
    ey(4,2) = 1.0;
    ey(2,3) = 1.0;
    ey(3,3) = 1.0;
    ey(4,3) = 1.0;
    ey(2,4) = 1.0;
    ey(3,4) = 1.0;
    ey(4,4) = 1.0;

    bz.reset();
    bz(2,2) = 1.0;
    bz(3,2) = 1.0;
    bz(4,2) = 1.0;
    bz(2,3) = 1.0;
    bz(3,3) = 1.0;
    bz(4,3) = 1.0;
    bz(2,4) = 1.0;
    bz(3,4) = 1.0;
    bz(4,4) = 1.0;

    ex.gather(P.E, P.R, 0, IntS1());
    ey.gather(P.E, P.R, 1, IntS1());
    bz.gather(P.B, P.R, IntS1());

    testmsg << "-- Testing IntS1 gather on edge centered fields" << endl;
    if (std::abs(P.E[0](0) - 1.0) > 1e-8 ||
        std::abs(P.E[0](1) - 1.0) > 1e-8) {
        ERRORMSG("S1 gather not working for edge centered fields" << endl);
        errors = true;
    }

    testmsg << "-- Testing IntS1 gather on cell centered fields" << endl;
    if (std::abs(P.B[0] - 1.0) > 1e-8) {
        ERRORMSG("S1 gather not working for cell centered fields" << endl);
        errors = true;
    }

    // CHECK S2 /////////////////////////////////////////////////////////////
    EFDcheck.assign({Vcheck(1,2,Vector_t(0.000113777,0.000000000)),
                     Vcheck(1,3,Vector_t(0.000840888,0.000000000)),
                     Vcheck(1,4,Vector_t(0.000376888,0.000000000)),
                     Vcheck(1,5,Vector_t(0.000001777,0.000000000)),
                     Vcheck(2,1,Vector_t(0.000000000,0.000257250)),
                     Vcheck(2,2,Vector_t(0.024120888,0.019903527)),
                     Vcheck(2,3,Vector_t(0.178268444,0.033737861)),
                     Vcheck(2,4,Vector_t(0.079900444,0.003268027)),
                     Vcheck(2,5,Vector_t(0.000376888,0.000000000)),
                     Vcheck(3,1,Vector_t(0.000000000,0.002655750)),
                     Vcheck(3,2,Vector_t(0.053816888,0.205476361)),
                     Vcheck(3,3,Vector_t(0.397740444,0.348296694)),
                     Vcheck(3,4,Vector_t(0.178268444,0.033737861)),
                     Vcheck(3,5,Vector_t(0.000840888,0.000000000)),
                     Vcheck(4,1,Vector_t(0.000000000,0.001566750)),
                     Vcheck(4,2,Vector_t(0.007281777,0.121220027)),
                     Vcheck(4,3,Vector_t(0.053816888,0.205476361)),
                     Vcheck(4,4,Vector_t(0.024120888,0.019903527)),
                     Vcheck(4,5,Vector_t(0.000113777,0.000000000)),
                     Vcheck(5,1,Vector_t(0.000000000,0.000020250)),
                     Vcheck(5,2,Vector_t(0.000000000,0.001566750)),
                     Vcheck(5,3,Vector_t(0.000000000,0.002655750)),
                     Vcheck(5,4,Vector_t(0.000000000,0.000257250))});

    BFDcheck.assign({Scheck(1,1,0.000006000),
                     Scheck(1,2,0.000464222),
                     Scheck(1,3,0.000786888),
                     Scheck(1,4,0.000076222),
                     Scheck(2,1,0.001272000),
                     Scheck(2,2,0.098415111),
                     Scheck(2,3,0.166820444),
                     Scheck(2,4,0.016159111),
                     Scheck(3,1,0.002838000),
                     Scheck(3,2,0.219577111),
                     Scheck(3,3,0.372198444),
                     Scheck(3,4,0.036053111),
                     Scheck(4,1,0.000384000),
                     Scheck(4,2,0.029710222),
                     Scheck(4,3,0.050360888),
                     Scheck(4,4,0.004878222)});

    P.E[0] = Vector_t(1.0);
    P.B[0] = 1.0;

    ex.reset();
    ey.reset();
    bz.reset();

    ex.scatter(P.E, P.R, 0, IntS2());
    ey.scatter(P.E, P.R, 1, IntS2());
    bz.scatter(P.B, P.R, IntS2());

    testmsg << "-- Testing IntS2 scatter on edge centered fields" << endl;
    for (const Vcheck c: EFDcheck) {
        if (std::abs(ex(c.i,c.j) - c.value(0)) > 1e-8 ||
            std::abs(ey(c.i,c.j) - c.value(1)) > 1e-8) {
            std::cerr << c.i << ", " << c.j << "\t" << ex(c.i,c.j) << " / " << c.value(0) << "\n\t\t"
                      << ey(c.i,c.j) << " / " << c.value(1) << std::endl;
            ERRORMSG("S2 scatter interpolation not working for edge based fields" << endl);
            errors = true;
        }
    }


    testmsg << "-- Testing IntS2 scatter on cell centered fields" << endl;
    for (const Scheck c: BFDcheck) {
        if (std::abs(bz(c.i, c.j) - c.value) > 1e-8) {
            ERRORMSG("S2 scatter interpolation not working for cell centered fields" << endl);
            errors = true;
        }
    }

    P.E[0] = Vector_t(0.0);
    P.B[0] = 0.0;

    ex.reset();
    ex(1,2) = 1.0;
    ex(2,2) = 1.0;
    ex(3,2) = 1.0;
    ex(4,2) = 1.0;
    ex(1,3) = 1.0;
    ex(2,3) = 1.0;
    ex(3,3) = 1.0;
    ex(4,3) = 1.0;
    ex(1,4) = 1.0;
    ex(2,4) = 1.0;
    ex(3,4) = 1.0;
    ex(4,4) = 1.0;
    ex(1,5) = 1.0;
    ex(2,5) = 1.0;
    ex(3,5) = 1.0;
    ex(4,5) = 1.0;

    ey.reset();
    ey(2,1) = 1.0;
    ey(3,1) = 1.0;
    ey(4,1) = 1.0;
    ey(5,1) = 1.0;
    ey(2,2) = 1.0;
    ey(3,2) = 1.0;
    ey(4,2) = 1.0;
    ey(5,2) = 1.0;
    ey(2,3) = 1.0;
    ey(3,3) = 1.0;
    ey(4,3) = 1.0;
    ey(5,3) = 1.0;
    ey(2,4) = 1.0;
    ey(3,4) = 1.0;
    ey(4,4) = 1.0;
    ey(5,4) = 1.0;

    bz.reset();
    bz(1,1) = 1.0;
    bz(2,1) = 1.0;
    bz(3,1) = 1.0;
    bz(4,1) = 1.0;
    bz(1,2) = 1.0;
    bz(2,2) = 1.0;
    bz(3,2) = 1.0;
    bz(4,2) = 1.0;
    bz(1,3) = 1.0;
    bz(2,3) = 1.0;
    bz(3,3) = 1.0;
    bz(4,3) = 1.0;
    bz(1,4) = 1.0;
    bz(2,4) = 1.0;
    bz(3,4) = 1.0;
    bz(4,4) = 1.0;

    ex.gather(P.E, P.R, 0, IntS2());
    ey.gather(P.E, P.R, 1, IntS2());
    bz.gather(P.B, P.R, IntS2());

    testmsg << "-- Testing IntS2 gather on edge centered fields" << endl;
    if (std::abs(P.E[0](0) - 1.0) > 1e-8 ||
        std::abs(P.E[0](1) - 1.0) > 1e-8) {
        ERRORMSG("S2 gather not working for edge centered fields" << endl);
        errors = true;
    }

    testmsg << "-- Testing IntCIC gather on cell centered fields" << endl;
    if (std::abs(P.B[0] - 1.0) > 1e-8) {
        ERRORMSG("S2 gather not working for cell centered fields" << endl);
        errors = true;
    }

    if (errors) {
        ERRORMSG("ERRORS OCCURED" << endl);
        return 1;
    }

    testmsg << "Everything seems to work fine" << endl;
    return 0;
}

/***************************************************************************
 * $RCSfile: edge.cpp,v $   $Author: chkraus $
 * $Revision: 1.1.1.1 $   $Date: 2013/08/23 07:40:38 $
 * IPPL_VERSION_ID: $Id: edge.cpp,v 1.1.1.1 2013/08/23 07:40:38 chkraus Exp $
 ***************************************************************************/
