#include <initializer_list>
#include "defs.hh"
#include "utils.hh"
#include "Bend.hh"
#include "Distribution.hh"
#include "PartBunch.hh"
#include "FieldPatch.hh"
#include "Utility/RNGBitReverse.h"

void getCurrent(double J[],
                const double & tau,
                const Vector_t & a,
                const Vector_t & b,
                const Vector_t & qoverdxdt,
                const Vector_t & delta);

void getCurrentDensityImpl(VField_Edge_t& JFD,
                           const Vector_t & oldR,
                           const Vector_t & R,
                           const double & q,
                           const double & dt);

const int nx=20, ny=20;              // size of domain is nx X ny
const unsigned int totalP = 1;      // number of particles to create
const int nt = 0;                    // total number of timesteps

const double pi = acos(-1.0);
const double qmmax = 1.0;       // maximum value for particle q/m
const double dt = 1.0;          // size of timestep

std::ofstream dbg;

int main(int argc, char *argv[]){
    TAU_PROFILE("main()", "int (int, char **)", TAU_DEFAULT);

    dbg.open("calcCurrent_dbg.out");
    Ippl ippl(argc, argv);
    Inform testmsg(argv[0]);
    testmsg << "Particle test calcCurrent: Begin." << endl;

    bool errors = false;
    double dx[] = {0.04, 0.05};
    double dt = 1.0;

    // create layout objects
    Index I1(nx+1), J1(ny+1);
    Mesh_t mymesh(I1,J1,dx,Vector_t(0.0));

    FieldLayout_Vert_t FL_vert(mymesh);
    FieldLayout_Edge_t FL_edge(mymesh);

    VField_Edge_t J(mymesh, FL_edge, GuardCellSizes<DIM>(0U));
    SField_Vert_t rho(mymesh, FL_vert, GuardCellSizes<DIM>(0U));
    SField_Vert_t rho2(mymesh, FL_vert, GuardCellSizes<DIM>(0U));
    SField_Vert_t oldrho(mymesh, FL_vert, GuardCellSizes<DIM>(0U));
    rho = 0.0;
    rho2 = 0.0;
    oldrho = 0.0;
    J = 0.0;

    Bend domain(nx*dx[0]/2, ny*dx[1], 0.0, 0.0, 0.0);

    std::vector<Vector_t> oldR(1,Vector_t(0.164,0.36));
    std::vector<Vector_t> R(1,Vector_t(0.128,0.315));
    std::vector<Vector_t> P(1,0.0);
    std::vector<double> Q(1,-1.0);
    Distribution dist(R, oldR, P, Q);
    PartBunch bunch(FL_vert, FL_edge, mymesh, dist, domain, 1.0);

    bunch.Q.scatter(oldrho, bunch.oldR, IntOp_t());
    oldrho /= (dx[0] * dx[1]);
    bunch.Q[0] = 1.0;
    bunch.Q.scatter(rho, bunch.R, IntOp_t());
    rho /= (dx[0] * dx[1]);
    ShapeFunction<IntS1> sh;
    sh.getChargeDensityDiff(rho2, bunch);

    getCurrentDensityImpl(J, bunch.oldR[0], bunch.R[0], bunch.Q[0], dt);

    NDIndex<DIM> elem,elemmx,elemmy;
    for (int j = 1; j < ny + 1; ++ j) {
        elem[1] = elemmx[1] = Index(j,j);
        elemmy[1] = Index(j-1,j-1);
        for (int i = 1; i < nx + 1; ++ i) {
            elem[0] = elemmy[0] = Index(i,i);
            elemmx[0] = Index(i-1,i-1);
            const double locRho = rho2.localElement(elem);// + oldrho.localElement(elem);
            const Vector_t locJ = J.localElement(elem);
            if (std::abs(locRho) > SPACE_EPS || std::abs(locJ(0)) > SPACE_EPS || std::abs(locJ(1)) > SPACE_EPS) {
                double balance = (locRho / dt +
                                  (J.localElement(elem)(0) - J.localElement(elemmx)(0)) / dx[0] +
                                  (J.localElement(elem)(1) - J.localElement(elemmy)(1)) / dx[1]);

                if (std::abs(balance) > SPACE_EPS) {
                    errors = true;

                    std::cout << "node; <" << i << ", " << j << ">\n"
                              << "remainder: " << balance << ",\t "
                              << "ratio: " << std::abs(balance/locRho) * dt << "\033[0m\n\n";

                    std::cout << "\033[31;1m"
                              << "local rho: " << locRho / dt << "    " << rho.localElement(elem) + oldrho.localElement(elem)
                              << "    ratio: " << (rho.localElement(elem) + oldrho.localElement(elem)) / locRho
                              << "\033[0m\n";
                    std::cout << std::setw(30) << std::setprecision(8) << J.localElement(elem)(1) / dx[1] << "\n";
                    std::cout << std::setw(14) << std::setprecision(8) << J.localElement(elemmx)(0) / dx[0];
                    std::cout << std::setw(30) << std::setprecision(8) << J.localElement(elem)(0) / dx[0] << "\n";
                    std::cout << std::setw(30) << std::setprecision(8) << J.localElement(elemmy)(1) / dx[1] << "\n\n";
                }
            }
        }
    }

    if (errors) {
        testmsg << "\033[31;1mErrors occured!\033[0m" << endl;
    } else {
        testmsg << "Everything seems to work fine" << endl;
    }
    return 0;
}

void getCurrentDensityImpl(VField_Edge_t& JFD,
                           const Vector_t & oldR,
                           const Vector_t & R,
                           const double & q,
                           const double & dt)
{
    const Vector_t dx(JFD.get_mesh().get_meshSpacing(0),
                      JFD.get_mesh().get_meshSpacing(1));
    const Vector_t origin(JFD.get_mesh().get_origin());
    const Vector_t oneoverdx = Vector_t(1.0) / dx;
    const Vector_t oneoverdxdt = oneoverdx / dt;
    Utils::samePos uni(sqrt(dot(dx,dx)));

    NDIndex<DIM> lDom = JFD.getLayout().getLocalNDIndex();
    FieldPatch<double> jx(lDom);
    FieldPatch<double> jy(lDom);
    dbg << jx.size() << "\t" << jy.size() << std::endl;
    Utils::CompanionSorterItem xy_crossings[10];
    size_t numXs = Utils::calcCellBorderCrossings(xy_crossings, oldR, R, dx, origin, uni);

    const Vector_t Dx = R - oldR;
    const Vector_t delta = Dx * oneoverdx;
    const double coveredDistance = sqrt(dot(Dx, Dx));
    for (size_t j = 0; j < numXs - 1; ++ j) {
        const Vector_t xi = xy_crossings[j].subStep * oneoverdx;
        const Vector_t xf = xy_crossings[j+1].subStep * oneoverdx;
        const Vector_t segment = xy_crossings[j+1].subStep - xy_crossings[j].subStep;
        const double tau = sqrt(dot(segment, segment)) / coveredDistance;
        const Vector_t center = (xi + xf) / 2;
        const int I = static_cast<int>(floor(center(0) + 0.5));
        const int J = static_cast<int>(floor(center(1) + 0.5));
        const Vector_t a(I + 0.5 - xi(0), J + 0.5 - xi(1));
        const Vector_t b(I - xi(0), J - xi(1));

        double localCurrent[12];
        getCurrent(localCurrent, tau, a, b, q*oneoverdxdt, delta);

        jx(I-1,J-1) += localCurrent[0];
        jx(I  ,J-1) += localCurrent[1];
        jx(I-1,J  ) += localCurrent[2];
        jx(I  ,J  ) += localCurrent[3];
        jx(I-1,J+1) += localCurrent[4];
        jx(I  ,J+1) += localCurrent[5];

        jy(I-1,J-1) += localCurrent[6];
        jy(I  ,J-1) += localCurrent[7];
        jy(I+1,J-1) += localCurrent[8];
        jy(I-1,J  ) += localCurrent[9];
        jy(I  ,J  ) += localCurrent[10];
        jy(I+1,J  ) += localCurrent[11];


        dbg << I << ", " << J << "\t" << jx(I,J) << ", " << jy(I,J) << std::endl;

    }
    NDIndex<DIM> elem;
    for (int j = lDom[1].first(); j <= lDom[1].last(); ++ j) {
        elem[1] = Index(j,j);
        for (int i = lDom[0].first(); i <= lDom[0].last(); ++ i) {
            elem[0] = Index(i,i);
            JFD.localElement(elem) = Vector_t(jx(i,j), jy(i,j));
        }
    }

    JFD.accumGuardCells();
}

void getCurrent(double J[],
                const double & tau,
                const Vector_t & a,
                const Vector_t & b,
                const Vector_t & qoverdxdt,
                const Vector_t & delta)
{
    const Vector_t current_density(qoverdxdt(1) * delta(0),
                                   qoverdxdt(0) * delta(1));

    const static double half = 0.5;
    const static double onethird = 1.0 / 3.0;
    const static double onefourth = half * half;
    const double taupow2 = tau*tau * half;
    const double taupow3 = taupow2*tau * 2 * onethird;

    const double Ax1Im1 = a(0);
    const double Bx1Im1 = -delta(0) * tau;
    const double Ax1I = 1.0 - a(0);
    const double Bx1I = delta(0) * tau;

    const double Ay2Jm1 = half* (half+ b(1))*(half+ b(1)) * tau;
    const double By2Jm1 = -(half+ b(1)) * delta(1) * taupow2;
    const double Cy2Jm1 = half* delta(1)*delta(1) * taupow3;
    const double Ay2J = (0.75 - b(1)*b(1)) * tau;
    const double By2J = 2 * b(1) * delta(1) * taupow2;
    const double Cy2J = -delta(1)*delta(1) * taupow3;
    const double Ay2Jp1 = half* (half- b(1))*(half- b(1)) * tau;
    const double By2Jp1 = (half- b(1)) * delta(1) * taupow2;
    const double Cy2Jp1 = half* delta(1)*delta(1) * taupow3;

    const double Ax2Im1 = half* (half+ b(0))*(half+ b(0)) * tau;
    const double Bx2Im1 = -(half+ b(0)) * delta(0) * taupow2;
    const double Cx2Im1 = half* delta(0)*delta(0) * taupow3;
    const double Ax2I = (0.75 - b(0)*b(0)) * tau;
    const double Bx2I = 2 * b(0) * delta(0) * taupow2;
    const double Cx2I = -delta(0)*delta(0) * taupow3;
    const double Ax2Ip1 = half* (half- b(0))*(half- b(0)) * tau;
    const double Bx2Ip1 = (half- b(0)) * delta(0) * taupow2;
    const double Cx2Ip1 = half* delta(0)*delta(0) * taupow3;

    const double Ay1Jm1 = a(1);
    const double By1Jm1 = -delta(1) * tau;
    const double Ay1J = 1.0 - a(1);
    const double By1J = delta(1) * tau;

    // Jx: I-1, J-1
    J[0] = current_density(0) * ((Ax1Im1 + Bx1Im1) * (Ay2Jm1 + By2Jm1 + Cy2Jm1)
                                 - Bx1Im1 * (half * Ay2Jm1 + onethird * By2Jm1 + onefourth * Cy2Jm1));

    // Jx: I, J-1
    J[1] = current_density(0) * ((Ax1I + Bx1I) * (Ay2Jm1 + By2Jm1 + Cy2Jm1)
                                 - Bx1I * (half * Ay2Jm1 + onethird * By2Jm1 + onefourth * Cy2Jm1));

    // Jx: I-1, J
    J[2] = current_density(0) * ((Ax1Im1 + Bx1Im1) * (Ay2J + By2J + Cy2J)
                                 - Bx1Im1 * (half * Ay2J + onethird * By2J + onefourth * Cy2J));

    // Jx: I, J
    J[3] = current_density(0) * ((Ax1I + Bx1I) * (Ay2J + By2J + Cy2J)
                                 - Bx1I * (half * Ay2J + onethird * By2J + onefourth * Cy2J));

    // Jx: I-1, J+1
    J[4] = current_density(0) * ((Ax1Im1 + Bx1Im1) * (Ay2Jp1 + By2Jp1 + Cy2Jp1)
                                 - Bx1Im1 * (half * Ay2Jp1 + onethird * By2Jp1 + onefourth * Cy2Jp1));

    // Jx: I, J+1
    J[5] = current_density(0) * ((Ax1I + Bx1I) * (Ay2Jp1 + By2Jp1 + Cy2Jp1)
                                 - Bx1I * (half * Ay2Jp1 + onethird * By2Jp1 + onefourth * Cy2Jp1));

    // Jy: I-1, J-1
    J[6] = current_density(1) * ((Ay1Jm1 + By1Jm1) * (Ax2Im1 + Bx2Im1 + Cx2Im1)
                                 - By1Jm1 * (half * Ax2Im1 + onethird * Bx2Im1 + onefourth * Cx2Im1));

    // Jy: I, J-1
    J[7] = current_density(1) * ((Ay1Jm1 + By1Jm1) * (Ax2I + Bx2I + Cx2I)
                                 - By1Jm1 * (half * Ax2I + onethird * Bx2I + onefourth * Cx2I));

    // Jy: I+1, J-1
    J[8] = current_density(1) * ((Ay1Jm1 + By1Jm1) * (Ax2Ip1 + Bx2Ip1 + Cx2Ip1)
                                 - By1Jm1 * (half * Ax2Ip1 + onethird * Bx2Ip1 + onefourth * Cx2Ip1));

    // Jy: I-1, J
    J[9] = current_density(1) * ((Ay1J + By1J) * (Ax2Im1 + Bx2Im1 + Cx2Im1)
                                 - By1J * (half * Ax2Im1 + onethird * Bx2Im1 + onefourth * Cx2Im1));

    // Jy: I, J
    J[10] = current_density(1) * ((Ay1J + By1J) * (Ax2I + Bx2I + Cx2I)
                                  - By1J * (half * Ax2I + onethird * Bx2I + onefourth * Cx2I));

    // Jy: I+1, J
    J[11] = current_density(1) * ((Ay1J + By1J) * (Ax2Ip1 + Bx2Ip1 + Cx2Ip1)
                                  - By1J * (half * Ax2Ip1 + onethird * Bx2Ip1 + onefourth * Cx2Ip1));
}
/***************************************************************************
 * $RCSfile: edge.cpp,v $   $Author: chkraus $
 * $Revision: 1.1.1.1 $   $Date: 2013/08/23 07:40:38 $
 * IPPL_VERSION_ID: $Id: edge.cpp,v 1.1.1.1 2013/08/23 07:40:38 chkraus Exp $
 ***************************************************************************/
