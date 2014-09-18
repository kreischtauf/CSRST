#include <initializer_list>
#include "defs.hh"
#include "utils.hh"
#include "FieldPatch.hh"
#include "Utility/RNGBitReverse.h"

void getCurrent(double J[],
                const double & Q,
                const Vector_t & tau_i,
                const Vector_t & dxdt,
                const Vector_t & Dx);

void getCurrentDensityImpl(VField_Edge_t& JFD,
                           const Vector_t & oldR,
                           const Vector_t & R,
                           const double & q,
                           const double & dt);

template<class PL>
class ChargedParticles : public ParticleBase<PL> {
public:
    ParticleAttrib<double>     Q;       // charge-to-mass ratio
    ChargedParticles(PL* pl) : ParticleBase<PL>(pl) {
        // register the particle attributes
        this->addAttribute(Q);
  }
};

const int nx=20, ny=20;              // size of domain is nx X ny
const unsigned int totalP = 1;      // number of particles to create
const int nt = 0;                    // total number of timesteps

const double pi = acos(-1.0);
const double qmmax = 1.0;       // maximum value for particle q/m
const double dt = 1.0;          // size of timestep

std::ofstream dbg;

int main(int argc, char *argv[]){
    TAU_PROFILE("main()", "int (int, char **)", TAU_DEFAULT);

    Ippl ippl(argc, argv);
    Inform testmsg(argv[0]);
    testmsg << "Particle test calcCurrent: Begin." << endl;

    bool errors = false;
    double dx[] = {1.0, 0.5};
    double dt = 1.0;

    // create layout objects
    Index I1(nx+1), J1(ny+1);
    Mesh_t mymesh(I1,J1,dx,Vector_t(0.0));

    FieldLayout_Vert_t FL_vert(mymesh);
    FieldLayout_Edge_t FL_edge(mymesh);

    VField_Edge_t J(mymesh, FL_edge, GuardCellSizes<DIM>(0U));
    SField_Vert_t rho(mymesh, FL_vert, GuardCellSizes<DIM>(0U));
    rho = 0.0;
    J = 0.0;

    // create an empty ChargedParticles object, setting it to use periodic BC's
    playout_t* PL = new playout_t();
    ChargedParticles<playout_t> P(PL);
    P.create(2);
    P.R[0] = Vector_t(3.3,3.2);
    P.R[1] = Vector_t(4.2,3.3);
    P.Q[0] = -1.0;
    P.Q[1] = 1.0;
    P.Q.scatter(rho, P.R, IntOp_t());
    rho /= (dx[0] * dx[1]);

    getCurrentDensityImpl(J, P.R[0], P.R[1], P.Q[1], dt);

    NDIndex<DIM> elem,elemmx,elemmy;
    for (int j = 1; j < ny + 1; ++ j) {
        elem[1] = elemmx[1] = Index(j,j);
        elemmy[1] = Index(j-1,j-1);
        for (int i = 1; i < nx + 1; ++ i) {
            elem[0] = elemmy[0] = Index(i,i);
            elemmx[0] = Index(i-1,i-1);
            const double locRho = rho.localElement(elem);
            const Vector_t locJ = J.localElement(elem);
            if (std::abs(locRho) > SPACE_EPS || std::abs(locJ(0)) > SPACE_EPS || std::abs(locJ(1)) > SPACE_EPS) {
                double balance = (locRho / dt +
                                  (J.localElement(elem)(0) - J.localElement(elemmx)(0)) / dx[0] +
                                  (J.localElement(elem)(1) - J.localElement(elemmy)(1)) / dx[1]);

                if (std::abs(balance) > SPACE_EPS) {
                    errors = true;

                    std::cout << "node; <" << i << ", " << j << ">\n"
                              << "remainder: " << balance << ",\t "
                              << "\033[0mratio: " << std::abs(balance/locRho) * dt << "\033[0m\n\n";

                    std::cout << "local rho: " << locRho / dt << "\n";
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
    NDIndex<DIM> lDom = JFD.getLayout().getLocalNDIndex();
    const Vector_t dx(JFD.get_mesh().get_meshSpacing(0),
                      JFD.get_mesh().get_meshSpacing(1));
    const Vector_t origin(JFD.get_mesh().get_origin());
    Utils::samePos uni(sqrt(dot(dx,dx)));

    Utils::CompanionSorterItem xy_crossings[10];
    size_t numXs = Utils::calcCellBorderCrossings(xy_crossings, oldR, R, dx, origin, uni);

    Vector_t Dx = R - oldR;
    Vector_t rDx = Dx / dx;
    double coveredDistance = sqrt(dot(Dx, Dx));
    for (size_t i = 0; i < numXs - 1; ++ i) {
        const Vector_t & xi = xy_crossings[i].subStep / dx;
        const Vector_t & xf = xy_crossings[i+1].subStep / dx;

        const Vector_t center = (xi + xf) / 2;
        int I = static_cast<long>(floor(center(0)));
        int J = static_cast<long>(floor(center(1)));

        NDIndex<DIM> elem(Index(I,I), Index(J,J));

        const Vector_t segment = xy_crossings[i+1].subStep - xy_crossings[i].subStep;
        const double length_segment = sqrt(dot(segment, segment));
        const double relative_length_segment = length_segment / coveredDistance;

        Vector_t tau_i(I + 1 - xi(0), J + 1 - xi(1));

        const Vector_t start = xy_crossings[i].subStep - (oldR - origin);
        const Vector_t end = xy_crossings[i+1].subStep - (oldR - origin);
        // const Vector_t tau(sqrt(dot(start,start)) / coveredDistance,
        //                    sqrt(dot(end,end)) / coveredDistance);
        const double tau = (sqrt(dot(end,end)) - sqrt(dot(start,start))) / coveredDistance;
        const Vector_t a(I + 0.5 - xi(0), J + 0.5 - xi(1));
        const Vector_t b(I - xi(0), J - xi(1));
        const Vector_t delta = rDx;

        const double taupow2 = tau*tau / 2;
        const double taupow3 = taupow2*tau * 2 / 3;
        const double taupow4 = taupow2*taupow2;
        const double taupow5 = taupow3*taupow2 * 6 / 5;
        const double taupow6 = taupow3*taupow3 * 9 / 6;

        const double Ax2Im1 = (0.5 + a(0))*(0.5 + a(0)) / 2;
        const double Bx2Im1 = -(0.5 + a(0)) * delta(0);
        const double Cx2Im1 = delta(0)*delta(0) / 2;
        const double Ax2I = 0.75 - a(0)*a(0);
        const double Bx2I = 2 * a(0) * delta(0);
        const double Cx2I = -delta(0)*delta(0);
        const double Ax2Ip1 = (0.5 - a(0))*(0.5 - a(0)) / 2;
        const double Bx2Ip1 = (0.5 - a(0)) * delta(0);
        const double Cx2Ip1 = delta(0)*delta(0) / 2;

        const double Ay3Jm1 = (1 + b(1))*(1 + b(1))*(1 + b(1)) / 6;
        const double By3Jm1 = -(1 + b(1))*(1 + b(1)) * delta(1) / 2;
        const double Cy3Jm1 = (1 + b(1)) * delta(1)*delta(1) / 2;
        const double Dy3Jm1 = -delta(1)*delta(1)*delta(1) / 6;
        const double Ay3J = (4 - 3 * b(1)*b(1) * (2 + b(1))) / 6;
        const double By3J = b(1) * delta(1) * (3 * b(1) + 4) / 2;
        const double Cy3J = -delta(1)*delta(1) * (2 + 3 * b(1)) / 2;
        const double Dy3J = delta(1)*delta(1)*delta(1) / 2;
        const double Ay3Jp1 = (4 - 3 * (b(1) + 1)*(b(1) + 1) * (1 - b(1))) / 6;
        const double By3Jp1 = (1 - 2 * b(1) - 3 * b(1)*b(1)) * delta(1) / 2;
        const double Cy3Jp1 = (1 + 3 * b(1)) * delta(1)*delta(1) / 2;
        const double Dy3Jp1 = -delta(1)*delta(1)*delta(1) / 2;
        const double Ay3Jp2 = -b(1)*b(1)*b(1) / 6;
        const double By3Jp2 = b(1)*b(1) * delta(1) / 2;
        const double Cy3Jp2 = -b(1) * delta(1)*delta(1) / 2;
        const double Dy3Jp2 = delta(1)*delta(1)*delta(1) / 6;

        const double Ax3Im1 = (1 + b(0))*(1 + b(0))*(1 + b(0)) / 6;
        const double Bx3Im1 = -(1 + b(0))*(1 + b(0)) * delta(0) / 2;
        const double Cx3Im1 = (1 + b(0)) * delta(0)*delta(0) / 2;
        const double Dx3Im1 = -delta(0)*delta(0)*delta(0) / 6;
        const double Ax3I = (4 - 3 * b(0)*b(0) * (2 + b(0))) / 6;
        const double Bx3I = b(0) * delta(0) * (3 * b(0) + 4) / 2;
        const double Cx3I = -delta(0)*delta(0) * (2 + 3 * b(0)) / 2;
        const double Dx3I = delta(0)*delta(0)*delta(0) / 2;
        const double Ax3Ip1 = (4 - 3 * (b(0) + 1)*(b(0) + 1) * (1 - b(0))) / 6;
        const double Bx3Ip1 = (1 - 2 * b(0) - 3 * b(0)*b(0)) * delta(0) / 2;
        const double Cx3Ip1 = (1 + 3 * b(0)) * delta(0)*delta(0) / 2;
        const double Dx3Ip1 = -delta(0)*delta(0)*delta(0) / 2;
        const double Ax3Ip2 = -b(0)*b(0)*b(0) / 6;
        const double Bx3Ip2 = b(0)*b(0) * delta(0) / 2;
        const double Cx3Ip2 = -b(0) * delta(0)*delta(0) / 2;
        const double Dx3Ip2 = delta(0)*delta(0)*delta(0) / 6;

        const double Ay2Jm1 = (0.5 + a(1))*(0.5 + a(1)) / 2;
        const double By2Jm1 = -(0.5 + a(1)) * delta(1);
        const double Cy2Jm1 = delta(1)*delta(1) / 2;
        const double Ay2J = 0.75 - a(1)*a(1);
        const double By2J = 2 * a(1) * delta(1);
        const double Cy2J = -delta(1)*delta(1);
        const double Ay2Jp1 = (0.5 - a(1))*(0.5 - a(1)) / 2;
        const double By2Jp1 = (0.5 - a(1)) * delta(1);
        const double Cy2Jp1 = delta(1)*delta(1) / 2;

        double current[24];
        getCurrent(current, q, tau_i, dt * dx, relative_length_segment * rDx);
        double * Jx = current;
        double * Jy = current + 12;

        // I-1, J-1
        double newJxIm1Jm1 = q * rDx(0) / (dx(1) * dt) * (Cx2Im1 * Dy3Jm1                                       * taupow6 +
                                                          (Cx2Im1 * Cy3Jm1 + Bx2Im1 * Dy3Jm1)                   * taupow5 +
                                                          (Cx2Im1 * By3Jm1 + Bx2Im1 * Cy3Jm1 + Ax2Im1 * Dy3Jm1) * taupow4 +
                                                          (Cx2Im1 * Ay3Jm1 + Bx2Im1 * By3Jm1 + Ax2Im1 * Cy3Jm1) * taupow3 +
                                                          (                  Bx2Im1 * Ay3Jm1 + Ax2Im1 * By3Jm1) * taupow2 +
                                                          (                                    Ax2Im1 * Ay3Jm1) * tau);

        // I, J-1
        double newJxIJm1 = q * rDx(0) / (dx(1) * dt) * (Cx2I * Dy3Jm1                                   * taupow6 +
                                                        (Cx2I * Cy3Jm1 + Bx2I * Dy3Jm1)                 * taupow5 +
                                                        (Cx2I * By3Jm1 + Bx2I * Cy3Jm1 + Ax2I * Dy3Jm1) * taupow4 +
                                                        (Cx2I * Ay3Jm1 + Bx2I * By3Jm1 + Ax2I * Cy3Jm1) * taupow3 +
                                                        (                Bx2I * Ay3Jm1 + Ax2I * By3Jm1) * taupow2 +
                                                        (                                Ax2I * Ay3Jm1) * tau);

        // I+1, J-1
        double newJxIp1Jm1 = q * rDx(0) / (dx(1) * dt) * (Cx2Ip1 * Dy3Jm1                                       * taupow6 +
                                                          (Cx2Ip1 * Cy3Jm1 + Bx2Ip1 * Dy3Jm1)                   * taupow5 +
                                                          (Cx2Ip1 * By3Jm1 + Bx2Ip1 * Cy3Jm1 + Ax2Ip1 * Dy3Jm1) * taupow4 +
                                                          (Cx2Ip1 * Ay3Jm1 + Bx2Ip1 * By3Jm1 + Ax2Ip1 * Cy3Jm1) * taupow3 +
                                                          (                  Bx2Ip1 * Ay3Jm1 + Ax2Ip1 * By3Jm1) * taupow2 +
                                                          (                                    Ax2Ip1 * Ay3Jm1) * tau);

        // I-1, J
        double newJxIm1J = q * rDx(0) / (dx(1) * dt) * (Cx2Im1 * Dy3J                                   * taupow6 +
                                                        (Cx2Im1 * Cy3J + Bx2Im1 * Dy3J)                 * taupow5 +
                                                        (Cx2Im1 * By3J + Bx2Im1 * Cy3J + Ax2Im1 * Dy3J) * taupow4 +
                                                        (Cx2Im1 * Ay3J + Bx2Im1 * By3J + Ax2Im1 * Cy3J) * taupow3 +
                                                        (                Bx2Im1 * Ay3J + Ax2Im1 * By3J) * taupow2 +
                                                        (                                Ax2Im1 * Ay3J) * tau);

        // I, J
        double newJxIJ = q * rDx(0) / (dx(1) * dt) * (Cx2I * Dy3J                               * taupow6 +
                                                      (Cx2I * Cy3J + Bx2I * Dy3J)               * taupow5 +
                                                      (Cx2I * By3J + Bx2I * Cy3J + Ax2I * Dy3J) * taupow4 +
                                                      (Cx2I * Ay3J + Bx2I * By3J + Ax2I * Cy3J) * taupow3 +
                                                      (              Bx2I * Ay3J + Ax2I * By3J) * taupow2 +
                                                      (                            Ax2I * Ay3J) * tau);

        // I+1, J
        double newJxIp1J = q * rDx(0) / (dx(1) * dt) * (Cx2Ip1 * Dy3J                                   * taupow6 +
                                                        (Cx2Ip1 * Cy3J + Bx2Ip1 * Dy3J)                 * taupow5 +
                                                        (Cx2Ip1 * By3J + Bx2Ip1 * Cy3J + Ax2Ip1 * Dy3J) * taupow4 +
                                                        (Cx2Ip1 * Ay3J + Bx2Ip1 * By3J + Ax2Ip1 * Cy3J) * taupow3 +
                                                        (                Bx2Ip1 * Ay3J + Ax2Ip1 * By3J) * taupow2 +
                                                        (                                Ax2Ip1 * Ay3J) * tau);

        // I-1, J+1
        double newJxIm1Jp1 = q * rDx(0) / (dx(1) * dt) * (Cx2Im1 * Dy3Jp1                                       * taupow6 +
                                                          (Cx2Im1 * Cy3Jp1 + Bx2Im1 * Dy3Jp1)                   * taupow5 +
                                                          (Cx2Im1 * By3Jp1 + Bx2Im1 * Cy3Jp1 + Ax2Im1 * Dy3Jp1) * taupow4 +
                                                          (Cx2Im1 * Ay3Jp1 + Bx2Im1 * By3Jp1 + Ax2Im1 * Cy3Jp1) * taupow3 +
                                                          (                  Bx2Im1 * Ay3Jp1 + Ax2Im1 * By3Jp1) * taupow2 +
                                                          (                                    Ax2Im1 * Ay3Jp1) * tau);

        // I, J+1
        double newJxIJp1 = q * rDx(0) / (dx(1) * dt) * (Cx2I * Dy3Jp1                                   * taupow6 +
                                                        (Cx2I * Cy3Jp1 + Bx2I * Dy3Jp1)                 * taupow5 +
                                                        (Cx2I * By3Jp1 + Bx2I * Cy3Jp1 + Ax2I * Dy3Jp1) * taupow4 +
                                                        (Cx2I * Ay3Jp1 + Bx2I * By3Jp1 + Ax2I * Cy3Jp1) * taupow3 +
                                                        (                Bx2I * Ay3Jp1 + Ax2I * By3Jp1) * taupow2 +
                                                        (                                Ax2I * Ay3Jp1) * tau);

        // I+1, J+1
        double newJxIp1Jp1 = q * rDx(0) / (dx(1) * dt) * (Cx2Ip1 * Dy3Jp1                                       * taupow6 +
                                                          (Cx2Ip1 * Cy3Jp1 + Bx2Ip1 * Dy3Jp1)                   * taupow5 +
                                                          (Cx2Ip1 * By3Jp1 + Bx2Ip1 * Cy3Jp1 + Ax2Ip1 * Dy3Jp1) * taupow4 +
                                                          (Cx2Ip1 * Ay3Jp1 + Bx2Ip1 * By3Jp1 + Ax2Ip1 * Cy3Jp1) * taupow3 +
                                                          (                  Bx2Ip1 * Ay3Jp1 + Ax2Ip1 * By3Jp1) * taupow2 +
                                                          (                                    Ax2Ip1 * Ay3Jp1) * tau);

        // I-1, J+2
        double newJxIm1Jp2 = q * rDx(0) / (dx(1) * dt) * (Cx2Im1 * Dy3Jp2                                       * taupow6 +
                                                          (Cx2Im1 * Cy3Jp2 + Bx2Im1 * Dy3Jp2)                   * taupow5 +
                                                          (Cx2Im1 * By3Jp2 + Bx2Im1 * Cy3Jp2 + Ax2Im1 * Dy3Jp2) * taupow4 +
                                                          (Cx2Im1 * Ay3Jp2 + Bx2Im1 * By3Jp2 + Ax2Im1 * Cy3Jp2) * taupow3 +
                                                          (                  Bx2Im1 * Ay3Jp2 + Ax2Im1 * By3Jp2) * taupow2 +
                                                          (                                    Ax2Im1 * Ay3Jp2) * tau);

        // I, J+2
        double newJxIJp2 = q * rDx(0) / (dx(1) * dt) * (Cx2I * Dy3Jp2                                   * taupow6 +
                                                        (Cx2I * Cy3Jp2 + Bx2I * Dy3Jp2)                 * taupow5 +
                                                        (Cx2I * By3Jp2 + Bx2I * Cy3Jp2 + Ax2I * Dy3Jp2) * taupow4 +
                                                        (Cx2I * Ay3Jp2 + Bx2I * By3Jp2 + Ax2I * Cy3Jp2) * taupow3 +
                                                        (                Bx2I * Ay3Jp2 + Ax2I * By3Jp2) * taupow2 +
                                                        (                                Ax2I * Ay3Jp2) * tau);

        // I+1, J+2
        double newJxIp1Jp2 = q * rDx(0) / (dx(1) * dt) * (Cx2Ip1 * Dy3Jp2                                       * taupow6 +
                                                          (Cx2Ip1 * Cy3Jp2 + Bx2Ip1 * Dy3Jp2)                   * taupow5 +
                                                          (Cx2Ip1 * By3Jp2 + Bx2Ip1 * Cy3Jp2 + Ax2Ip1 * Dy3Jp2) * taupow4 +
                                                          (Cx2Ip1 * Ay3Jp2 + Bx2Ip1 * By3Jp2 + Ax2Ip1 * Cy3Jp2) * taupow3 +
                                                          (                  Bx2Ip1 * Ay3Jp2 + Ax2Ip1 * By3Jp2) * taupow2 +
                                                          (                                    Ax2Ip1 * Ay3Jp2) * tau);

        // I-1, J-1
        double newJyIm1Jm1 = q * rDx(1) / (dx(0) * dt) * (Cy2Jm1 * Dx3Im1                                       * taupow6 +
                                                          (Cy2Jm1 * Cx3Im1 + By2Jm1 * Dx3Im1)                   * taupow5 +
                                                          (Cy2Jm1 * Bx3Im1 + By2Jm1 * Cx3Im1 + Ay2Jm1 * Dx3Im1) * taupow4 +
                                                          (Cy2Jm1 * Ax3Im1 + By2Jm1 * Bx3Im1 + Ay2Jm1 * Cx3Im1) * taupow3 +
                                                          (                  By2Jm1 * Ax3Im1 + Ay2Jm1 * Bx3Im1) * taupow2 +
                                                          (                                    Ay2Jm1 * Ax3Im1) * tau);

        // I-1, J
        double newJyIm1J = q * rDx(1) / (dx(0) * dt) * (Cy2J * Dx3Im1                                   * taupow6 +
                                                        (Cy2J * Cx3Im1 + By2J * Dx3Im1)                 * taupow5 +
                                                        (Cy2J * Bx3Im1 + By2J * Cx3Im1 + Ay2J * Dx3Im1) * taupow4 +
                                                        (Cy2J * Ax3Im1 + By2J * Bx3Im1 + Ay2J * Cx3Im1) * taupow3 +
                                                        (                By2J * Ax3Im1 + Ay2J * Bx3Im1) * taupow2 +
                                                        (                                Ay2J * Ax3Im1) * tau);

        // I-1, J+1
        double newJyIm1Jp1 = q * rDx(1) / (dx(0) * dt) * (Cy2Jp1 * Dx3Im1                                       * taupow6 +
                                                          (Cy2Jp1 * Cx3Im1 + By2Jp1 * Dx3Im1)                   * taupow5 +
                                                          (Cy2Jp1 * Bx3Im1 + By2Jp1 * Cx3Im1 + Ay2Jp1 * Dx3Im1) * taupow4 +
                                                          (Cy2Jp1 * Ax3Im1 + By2Jp1 * Bx3Im1 + Ay2Jp1 * Cx3Im1) * taupow3 +
                                                          (                  By2Jp1 * Ax3Im1 + Ay2Jp1 * Bx3Im1) * taupow2 +
                                                          (                                    Ay2Jp1 * Ax3Im1) * tau);

        // I, J-1
        double newJyIJm1 = q * rDx(1) / (dx(0) * dt) * (Cy2Jm1 * Dx3I                                   * taupow6 +
                                                        (Cy2Jm1 * Cx3I + By2Jm1 * Dx3I)                 * taupow5 +
                                                        (Cy2Jm1 * Bx3I + By2Jm1 * Cx3I + Ay2Jm1 * Dx3I) * taupow4 +
                                                        (Cy2Jm1 * Ax3I + By2Jm1 * Bx3I + Ay2Jm1 * Cx3I) * taupow3 +
                                                        (                By2Jm1 * Ax3I + Ay2Jm1 * Bx3I) * taupow2 +
                                                        (                                Ay2Jm1 * Ax3I) * tau);

        // I, J
        double newJyIJ = q * rDx(1) / (dx(0) * dt) * (Cy2J * Dx3I                               * taupow6 +
                                                      (Cy2J * Cx3I + By2J * Dx3I)               * taupow5 +
                                                      (Cy2J * Bx3I + By2J * Cx3I + Ay2J * Dx3I) * taupow4 +
                                                      (Cy2J * Ax3I + By2J * Bx3I + Ay2J * Cx3I) * taupow3 +
                                                      (              By2J * Ax3I + Ay2J * Bx3I) * taupow2 +
                                                      (                            Ay2J * Ax3I) * tau);

        // I, J+1
        double newJyIJp1 = q * rDx(1) / (dx(0) * dt) * (Cy2Jp1 * Dx3I                                   * taupow6 +
                                                        (Cy2Jp1 * Cx3I + By2Jp1 * Dx3I)                 * taupow5 +
                                                        (Cy2Jp1 * Bx3I + By2Jp1 * Cx3I + Ay2Jp1 * Dx3I) * taupow4 +
                                                        (Cy2Jp1 * Ax3I + By2Jp1 * Bx3I + Ay2Jp1 * Cx3I) * taupow3 +
                                                        (                By2Jp1 * Ax3I + Ay2Jp1 * Bx3I) * taupow2 +
                                                        (                                Ay2Jp1 * Ax3I) * tau);

        // I+1, J-1
        double newJyIp1Jm1 = q * rDx(1) / (dx(0) * dt) * (Cy2Jm1 * Dx3Ip1                                       * taupow6 +
                                                          (Cy2Jm1 * Cx3Ip1 + By2Jm1 * Dx3Ip1)                   * taupow5 +
                                                          (Cy2Jm1 * Bx3Ip1 + By2Jm1 * Cx3Ip1 + Ay2Jm1 * Dx3Ip1) * taupow4 +
                                                          (Cy2Jm1 * Ax3Ip1 + By2Jm1 * Bx3Ip1 + Ay2Jm1 * Cx3Ip1) * taupow3 +
                                                          (                  By2Jm1 * Ax3Ip1 + Ay2Jm1 * Bx3Ip1) * taupow2 +
                                                          (                                    Ay2Jm1 * Ax3Ip1) * tau);

        // I+1, J
        double newJyIp1J = q * rDx(1) / (dx(0) * dt) * (Cy2J * Dx3Ip1                                   * taupow6 +
                                                        (Cy2J * Cx3Ip1 + By2J * Dx3Ip1)                 * taupow5 +
                                                        (Cy2J * Bx3Ip1 + By2J * Cx3Ip1 + Ay2J * Dx3Ip1) * taupow4 +
                                                        (Cy2J * Ax3Ip1 + By2J * Bx3Ip1 + Ay2J * Cx3Ip1) * taupow3 +
                                                        (                By2J * Ax3Ip1 + Ay2J * Bx3Ip1) * taupow2 +
                                                        (                                Ay2J * Ax3Ip1) * tau);

        // I+1, J+1
        double newJyIp1Jp1 = q * rDx(1) / (dx(0) * dt) * (Cy2Jp1 * Dx3Ip1                                       * taupow6 +
                                                          (Cy2Jp1 * Cx3Ip1 + By2Jp1 * Dx3Ip1)                   * taupow5 +
                                                          (Cy2Jp1 * Bx3Ip1 + By2Jp1 * Cx3Ip1 + Ay2Jp1 * Dx3Ip1) * taupow4 +
                                                          (Cy2Jp1 * Ax3Ip1 + By2Jp1 * Bx3Ip1 + Ay2Jp1 * Cx3Ip1) * taupow3 +
                                                          (                  By2Jp1 * Ax3Ip1 + Ay2Jp1 * Bx3Ip1) * taupow2 +
                                                          (                                    Ay2Jp1 * Ax3Ip1) * tau);

        // I+2, J-1
        double newJyIp2Jm1 = q * rDx(1) / (dx(0) * dt) * (Cy2Jm1 * Dx3Ip2                                       * taupow6 +
                                                          (Cy2Jm1 * Cx3Ip2 + By2Jm1 * Dx3Ip2)                   * taupow5 +
                                                          (Cy2Jm1 * Bx3Ip2 + By2Jm1 * Cx3Ip2 + Ay2Jm1 * Dx3Ip2) * taupow4 +
                                                          (Cy2Jm1 * Ax3Ip2 + By2Jm1 * Bx3Ip2 + Ay2Jm1 * Cx3Ip2) * taupow3 +
                                                          (                  By2Jm1 * Ax3Ip2 + Ay2Jm1 * Bx3Ip2) * taupow2 +
                                                          (                                    Ay2Jm1 * Ax3Ip2) * tau);

        // I+2, J
        double newJyIp2J = q * rDx(1) / (dx(0) * dt) * (Cy2J * Dx3Ip2                                   * taupow6 +
                                                        (Cy2J * Cx3Ip2 + By2J * Dx3Ip2)                 * taupow5 +
                                                        (Cy2J * Bx3Ip2 + By2J * Cx3Ip2 + Ay2J * Dx3Ip2) * taupow4 +
                                                        (Cy2J * Ax3Ip2 + By2J * Bx3Ip2 + Ay2J * Cx3Ip2) * taupow3 +
                                                        (                By2J * Ax3Ip2 + Ay2J * Bx3Ip2) * taupow2 +
                                                        (                                Ay2J * Ax3Ip2) * tau);

        // I+2, J+1
        double newJyIp2Jp1 = q * rDx(1) / (dx(0) * dt) * (Cy2Jp1 * Dx3Ip2                                       * taupow6 +
                                                          (Cy2Jp1 * Cx3Ip2 + By2Jp1 * Dx3Ip2)                   * taupow5 +
                                                          (Cy2Jp1 * Bx3Ip2 + By2Jp1 * Cx3Ip2 + Ay2Jp1 * Dx3Ip2) * taupow4 +
                                                          (Cy2Jp1 * Ax3Ip2 + By2Jp1 * Bx3Ip2 + Ay2Jp1 * Cx3Ip2) * taupow3 +
                                                          (                  By2Jp1 * Ax3Ip2 + Ay2Jp1 * Bx3Ip2) * taupow2 +
                                                          (                                    Ay2Jp1 * Ax3Ip2) * tau);

        std::cout << Jx[0]  << "\t" << newJxIm1Jm1 << "\t" << Jx[0]  / newJxIm1Jm1 << "\n"
                  << Jx[1]  << "\t" << newJxIJm1   << "\t" << Jx[1]  / newJxIJm1   << "\n"
                  << Jx[2]  << "\t" << newJxIp1Jm1 << "\t" << Jx[2]  / newJxIp1Jm1 << "\n"
                  << Jx[3]  << "\t" << newJxIm1J   << "\t" << Jx[3]  / newJxIm1J   << "\n"
                  << Jx[4]  << "\t" << newJxIJ     << "\t" << Jx[4]  / newJxIJ     << "\n"
                  << Jx[5]  << "\t" << newJxIp1J   << "\t" << Jx[5]  / newJxIp1J   << "\n"
                  << Jx[6]  << "\t" << newJxIm1Jp1 << "\t" << Jx[6]  / newJxIm1Jp1 << "\n"
                  << Jx[7]  << "\t" << newJxIJp1   << "\t" << Jx[7]  / newJxIJp1   << "\n"
                  << Jx[8]  << "\t" << newJxIp1Jp1 << "\t" << Jx[8]  / newJxIp1Jp1 << "\n"
                  << Jx[9]  << "\t" << newJxIm1Jp2 << "\t" << Jx[9]  / newJxIm1Jp2 << "\n"
                  << Jx[10] << "\t" << newJxIJp2   << "\t" << Jx[10] / newJxIJp2   << "\n"
                  << Jx[11] << "\t" << newJxIp1Jp2 << "\t" << Jx[11] / newJxIp1Jp2 << "\n"
                  << Jy[0]  << "\t" << newJyIm1Jm1 << "\t" << Jy[0]  / newJyIm1Jm1 << "\n"
                  << Jy[1]  << "\t" << newJyIm1J   << "\t" << Jy[1]  / newJyIm1J   << "\n"
                  << Jy[2]  << "\t" << newJyIm1Jp1 << "\t" << Jy[2]  / newJyIm1Jp1 << "\n"
                  << Jy[3]  << "\t" << newJyIJm1   << "\t" << Jy[3]  / newJyIJm1   << "\n"
                  << Jy[4]  << "\t" << newJyIJ     << "\t" << Jy[4]  / newJyIJ     << "\n"
                  << Jy[5]  << "\t" << newJyIJp1   << "\t" << Jy[5]  / newJyIJp1   << "\n"
                  << Jy[6]  << "\t" << newJyIp1Jm1 << "\t" << Jy[6]  / newJyIp1Jm1 << "\n"
                  << Jy[7]  << "\t" << newJyIp1J   << "\t" << Jy[7]  / newJyIp1J   << "\n"
                  << Jy[8]  << "\t" << newJyIp1Jp1 << "\t" << Jy[8]  / newJyIp1Jp1 << "\n"
                  << Jy[9]  << "\t" << newJyIp2Jm1 << "\t" << Jy[9]  / newJyIp2Jm1 << "\n"
                  << Jy[10] << "\t" << newJyIp2J   << "\t" << Jy[10] / newJyIp2J   << "\n"
                  << Jy[11] << "\t" << newJyIp2Jp1 << "\t" << Jy[11] / newJyIp2Jp1 << "\n"
                  << std::endl;

        if (I > lDom[0].first() + 1 && I < lDom[0].last() - 1&&
            J > lDom[1].first() + 1 && J < lDom[1].last() - 1) {

            elem[0] = Index(I - 1, I - 1);
            elem[1] = Index(J - 1, J - 1);
            JFD.localElement(elem) += Vector_t(Jx[0], Jy[0]);

            elem[0] = Index(I, I);
            JFD.localElement(elem) += Vector_t(Jx[1], Jy[3]);

            elem[0] = Index(I + 1, I + 1);
            JFD.localElement(elem) += Vector_t(Jx[2], Jy[6]);

            elem[0] = Index(I + 2, I + 2);
            JFD.localElement(elem)(1) += Jy[9];

            elem[1] = Index(J, J);
            elem[0] = Index(I - 1, I - 1);
            JFD.localElement(elem) += Vector_t(Jx[3], Jy[1]);

            elem[0] = Index(I, I);
            JFD.localElement(elem) += Vector_t(Jx[4], Jy[4]);

            elem[0] = Index(I + 1, I + 1);
            JFD.localElement(elem) += Vector_t(Jx[5], Jy[7]);

            elem[0] = Index(I + 2, I + 2);
            JFD.localElement(elem)(1) += Jy[10];

            elem[1] = Index(J + 1, J + 1);
            elem[0] = Index(I - 1, I - 1);
            JFD.localElement(elem) += Vector_t(Jx[6], Jy[2]);

            elem[0] = Index(I, I);
            JFD.localElement(elem) += Vector_t(Jx[7], Jy[5]);

            elem[0] = Index(I + 1, I + 1);
            JFD.localElement(elem) += Vector_t(Jx[8], Jy[8]);

            elem[0] = Index(I + 2, I + 2);
            JFD.localElement(elem)(1) += Jy[11];

            elem[1] = Index(J + 2, J + 2);
            elem[0] = Index(I - 1, I - 1);
            JFD.localElement(elem)(0) += Jx[9];

            elem[0] = Index(I, I);
            JFD.localElement(elem)(0) += Jx[10];

            elem[0] = Index(I + 1, I + 1);
            JFD.localElement(elem)(0) += Jx[11];

        } else {
            int k = 0;
            for (int j = std::max(J - 1, lDom[1].first());
                 j <= std::min(J + 2, lDom[1].last());
                 ++ j) {
                elem[1] = Index(j, j);
                for (int i = std::max(I - 1, lDom[0].first());
                     i <= std::min(I + 1, lDom[0].last());
                     ++ i) {
                    elem[0] = Index(i, i);

                    k = (j - J) * 3 + i - I + 4;
                    JFD.localElement(elem)(0) += Jx[k];
                }
            }
            for (int i = std::max(I - 1, lDom[0].first());
                 i <= std::min(I + 2, lDom[0].last());
                 ++ i) {
                elem[0] = Index(i, i);
                for (int j = std::max(J - 1, lDom[1].first());
                     j <= std::min(J + 1, lDom[1].last());
                     ++ j) {
                    elem[1] = Index(j, j);

                    k = (i - I) * 3 + j - J + 4;
                    JFD.localElement(elem)(1) += Jy[k];
                }
            }
        }
    }
}


void getCurrent(double J[],
                const double & Q,
                const Vector_t & tau_i,
                const Vector_t & dxdt,
                const Vector_t & Dx)
{
    const Vector_t current_density(Q * Dx(0) / dxdt(1),
                                   Q * Dx(1) / dxdt(0));
    double as = tau_i(0)*tau_i(0), at = as * tau_i(0), bs = tau_i(1)*tau_i(1), bt = bs * tau_i(1);
    double c = tau_i(0) - 1, cs = c*c, ct = cs * c;
    double e = tau_i(1) - 1, es = e*e, et = es * e;
    double dxs = Dx(0)*Dx(0), dxt = dxs * Dx(0), dys = Dx(1)*Dx(1), dyt = dys * Dx(1);

    J[0] = current_density(0) *
        ((as*bt)/12. - (tau_i(0)*bt*Dx(0))/12. + (bt*dxs)/36. +
         (-(as*bs)/8. + (tau_i(0)*bs*Dx(0))/6. - (bs*dxs)/16.)*Dx(1) +
         ((as*tau_i(1))/12. - (tau_i(0)*tau_i(1)*Dx(0))/8. + (tau_i(1)*dxs)/20.)*dys +
         (-as/48. + (tau_i(0)*Dx(0))/30. - dxs/72.)*dyt);

    J[1] = current_density(0) *
        (-(bt*(-3 + 6*as - 6*tau_i(0)*(1 + Dx(0)) + Dx(0)*(3 + 2*Dx(0))))/36. +
         (bs*(-3 + 6*as + Dx(0)*(4 + 3*Dx(0)) - 2*tau_i(0)*(3 + 4*Dx(0)))*Dx(1))/24. +
         (tau_i(1)*(10 - 20*as + 10*tau_i(0)*(2 + 3*Dx(0)) - 3*Dx(0)*(5 + 4*Dx(0)))*dys)/120. +
         ((-15 + 30*as + 4*Dx(0)*(6 + 5*Dx(0)) - 6*tau_i(0)*(5 + 8*Dx(0)))*dyt)/720.);

    J[2] = current_density(0) *
        ((bt*(3 + 3*as - 3*tau_i(0)*(2 + Dx(0)) + Dx(0)*(3 + Dx(0))))/36. -
         (bs*(6*cs - 8*c*Dx(0) + 3*dxs)*Dx(1))/48. +
         (tau_i(1)*(10*cs - 15*c*Dx(0) + 6*dxs)*dys)/120. +
         ((-15*cs + 24*c*Dx(0) - 10*dxs)*dyt)/720.);

    J[3] = current_density(0) *
        (-((-1 + 3*tau_i(1)*(-1 + e*tau_i(1)))*(3*as - 3*tau_i(0)*Dx(0) + dxs))/36. +
         (e*(1 + 3*tau_i(1))*(6*as - 8*tau_i(0)*Dx(0) + 3*dxs)*Dx(1))/48. -
         ((-1 + 3*tau_i(1))*(10*as - 15*tau_i(0)*Dx(0) + 6*dxs)*dys)/120. +
         ((15*as - 24*tau_i(0)*Dx(0) + 10*dxs)*dyt)/240.);

    J[4] = current_density(0) *
        (((-1 + 3*tau_i(1)*(-1 + e*tau_i(1)))*(-3 + 6*as - 6*tau_i(0)*(1 + Dx(0)) + Dx(0)*(3 + 2*Dx(0))))/36. -
         (e*(1 + 3*tau_i(1))*(-3 + 6*c*tau_i(0) + 4*Dx(0) - 8*tau_i(0)*Dx(0) + 3*dxs)*Dx(1))/24. +
         ((-1 + 3*tau_i(1))*(-10 + 20*as - 10*tau_i(0)*(2 + 3*Dx(0)) + 3*Dx(0)*(5 + 4*Dx(0)))*dys)/120. +
         ((15 - 30*as - 4*Dx(0)*(6 + 5*Dx(0)) + 6*tau_i(0)*(5 + 8*Dx(0)))*dyt)/240.);

    J[5] = current_density(0) *
        (-((-1 + 3*tau_i(1)*(-1 + e*tau_i(1)))*(3 + 3*as - 3*tau_i(0)*(2 + Dx(0)) + Dx(0)*(3 + Dx(0))))/36. +
         (e*(1 + 3*tau_i(1))*(6*cs - 8*c*Dx(0) + 3*dxs)*Dx(1))/48. -
         ((-1 + 3*tau_i(1))*(10*cs - 15*c*Dx(0) + 6*dxs)*dys)/120. +
         ((15*cs - 24*c*Dx(0) + 10*dxs)*dyt)/240.);

    J[6] = current_density(0) *
        (((4 + 3*(-2 + tau_i(1))*bs)*(3*as - 3*tau_i(0)*Dx(0) + dxs))/36. -
         (tau_i(1)*(-4 + 3*tau_i(1))*(6*as - 8*tau_i(0)*Dx(0) + 3*dxs)*Dx(1))/48. +
         ((-2 + 3*tau_i(1))*(10*as - 15*tau_i(0)*Dx(0) + 6*dxs)*dys)/120. -
         ((15*as - 24*tau_i(0)*Dx(0) + 10*dxs)*dyt)/240.);

    J[7] = current_density(0) *
        (-((4 + 3*(-2 + tau_i(1))*bs)*(-3 + 6*as - 6*tau_i(0)*(1 + Dx(0)) + Dx(0)*(3 + 2*Dx(0))))/36. +
         (tau_i(1)*(-4 + 3*tau_i(1))*(-3 + 6*c*tau_i(0) + 4*Dx(0) - 8*tau_i(0)*Dx(0) + 3*dxs)*Dx(1))/24. -
         ((-2 + 3*tau_i(1))*(-10 + 20*as - 10*tau_i(0)*(2 + 3*Dx(0)) + 3*Dx(0)*(5 + 4*Dx(0)))*dys)/120. +
         ((-15 + 30*as + 4*Dx(0)*(6 + 5*Dx(0)) - 6*tau_i(0)*(5 + 8*Dx(0)))*dyt)/240.);

    J[8] = current_density(0) *
        (((4 + 3*(-2 + tau_i(1))*bs)*(3 + 3*as - 3*tau_i(0)*(2 + Dx(0)) + Dx(0)*(3 + Dx(0))))/36. -
         (tau_i(1)*(-4 + 3*tau_i(1))*(6*cs - 8*c*Dx(0) + 3*dxs)*Dx(1))/48. +
         ((-2 + 3*tau_i(1))*(10*cs - 15*c*Dx(0) + 6*dxs)*dys)/120. -
         ((15*cs - 24*c*Dx(0) + 10*dxs)*dyt)/240.);

    J[9] = current_density(0) *
        (-(et*(3*as - 3*tau_i(0)*Dx(0) + dxs))/36. +
         (es*(6*as - 8*tau_i(0)*Dx(0) + 3*dxs)*Dx(1))/48. -
         (e*(10*as - 15*tau_i(0)*Dx(0) + 6*dxs)*dys)/120. +
         ((15*as - 24*tau_i(0)*Dx(0) + 10*dxs)*dyt)/720.);

    J[10] = current_density(0) *
        ((et*(-3 + 6*as - 6*tau_i(0)*(1 + Dx(0)) + Dx(0)*(3 + 2*Dx(0))))/36. -
         (es*(-3 + 6*c*tau_i(0) + 4*Dx(0) - 8*tau_i(0)*Dx(0) + 3*dxs)*Dx(1))/24. +
         (e*(-10 + 20*as - 10*tau_i(0)*(2 + 3*Dx(0)) + 3*Dx(0)*(5 + 4*Dx(0)))*dys)/120. +
         ((15 - 30*as - 4*Dx(0)*(6 + 5*Dx(0)) + 6*tau_i(0)*(5 + 8*Dx(0)))*dyt)/720.);

    J[11] = current_density(0) *
        (-(et*(3 + 3*as - 3*tau_i(0)*(2 + Dx(0)) + Dx(0)*(3 + Dx(0))))/36. +
         (es*(6*cs - 8*c*Dx(0) + 3*dxs)*Dx(1))/48. -
         (e*(10*cs - 15*c*Dx(0) + 6*dxs)*dys)/120. +
         ((15*cs - 24*c*Dx(0) + 10*dxs)*dyt)/720.);

    J[12] = current_density(1) *
        ((at*bs)/12. - (at*tau_i(1)*Dx(1))/12. + (at*dys)/36. +
         dxt*(-bs/48. + (tau_i(1)*Dx(1))/30. - dys/72.) +
         dxs*((tau_i(0)*bs)/12. - (tau_i(0)*tau_i(1)*Dx(1))/8. + (tau_i(0)*dys)/20.) + Dx(0)*
         (-(as*bs)/8. + (as*tau_i(1)*Dx(1))/6. -
          (as*dys)/16.));

    J[13] = current_density(1) *
        (-(at*(-3 + 6*bs - 6*tau_i(1)*(1 + Dx(1)) + Dx(1)*(3 + 2*Dx(1))))/36. +
         (as*Dx(0)*(-3 + 6*bs + Dx(1)*(4 + 3*Dx(1)) - 2*tau_i(1)*(3 + 4*Dx(1))))/24. +
         (tau_i(0)*dxs*(10 - 20*bs + 10*tau_i(1)*(2 + 3*Dx(1)) - 3*Dx(1)*(5 + 4*Dx(1))))/120. +
         (dxt*(-15 + 30*bs + 4*Dx(1)*(6 + 5*Dx(1)) - 6*tau_i(1)*(5 + 8*Dx(1))))/720.);

    J[14] = current_density(1) *
        ((dxt*(-15*es + 24*e*Dx(1) - 10*dys))/720. -
         (as*Dx(0)*(6*es - 8*e*Dx(1) + 3*dys))/48. +
         (tau_i(0)*dxs*(10*es - 15*e*Dx(1) + 6*dys))/120. +
         (at*(3 + 3*bs - 3*tau_i(1)*(2 + Dx(1)) + Dx(1)*(3 + Dx(1))))/36.);

    J[15] = current_density(1) *
        (-((-1 + 3*tau_i(0)*(-1 + c*tau_i(0)))*(3*bs - 3*tau_i(1)*Dx(1) + dys))/36. +
         (c*(1 + 3*tau_i(0))*Dx(0)*(6*bs - 8*tau_i(1)*Dx(1) + 3*dys))/48. -
         ((-1 + 3*tau_i(0))*dxs*(10*bs - 15*tau_i(1)*Dx(1) + 6*dys))/120. +
         (dxt*(15*bs - 24*tau_i(1)*Dx(1) + 10*dys))/240.);

    J[16] = current_density(1) *
        (-(c*(1 + 3*tau_i(0))*Dx(0)*(-3 + 6*e*tau_i(1) + 4*Dx(1) - 8*tau_i(1)*Dx(1) + 3*dys))/24. +
         ((-1 + 3*tau_i(0)*(-1 + c*tau_i(0)))*(-3 + 6*bs - 6*tau_i(1)*(1 + Dx(1)) + Dx(1)*(3 + 2*Dx(1))))/36. +
         ((-1 + 3*tau_i(0))*dxs*(-10 + 20*bs - 10*tau_i(1)*(2 + 3*Dx(1)) + 3*Dx(1)*(5 + 4*Dx(1))))/120. +
         (dxt*(15 - 30*bs - 4*Dx(1)*(6 + 5*Dx(1)) + 6*tau_i(1)*(5 + 8*Dx(1))))/ 240.);

    J[17] = current_density(1) *
        ((c*(1 + 3*tau_i(0))*Dx(0)*(6*es - 8*e*Dx(1) + 3*dys))/48. -
         ((-1 + 3*tau_i(0))*dxs*(10*es - 15*e*Dx(1) + 6*dys))/120. +
         (dxt*(15*es - 24*e*Dx(1) + 10*dys))/240. -
         ((-1 + 3*tau_i(0)*(-1 + c*tau_i(0)))*(3 + 3*bs - 3*tau_i(1)*(2 + Dx(1)) + Dx(1)*(3 + Dx(1))))/36.);

    J[18] = current_density(1) *
        (((4 + 3*(-2 + tau_i(0))*as)*(3*bs - 3*tau_i(1)*Dx(1) + dys))/36. -
         (tau_i(0)*(-4 + 3*tau_i(0))*Dx(0)*(6*bs - 8*tau_i(1)*Dx(1) + 3*dys))/48. +
         ((-2 + 3*tau_i(0))*dxs*(10*bs - 15*tau_i(1)*Dx(1) + 6*dys))/120. -
         (dxt*(15*bs - 24*tau_i(1)*Dx(1) + 10*dys))/240.);

    J[19] = current_density(1) *
        ((tau_i(0)*(-4 + 3*tau_i(0))*Dx(0)*(-3 + 6*e*tau_i(1) + 4*Dx(1) - 8*tau_i(1)*Dx(1) + 3*dys))/24. -
         ((4 + 3*(-2 + tau_i(0))*as)*(-3 + 6*bs - 6*tau_i(1)*(1 + Dx(1)) + Dx(1)*(3 + 2*Dx(1))))/36. -
         ((-2 + 3*tau_i(0))*dxs*(-10 + 20*bs - 10*tau_i(1)*(2 + 3*Dx(1)) + 3*Dx(1)*(5 + 4*Dx(1))))/120. +
         (dxt*(-15 + 30*bs + 4*Dx(1)*(6 + 5*Dx(1)) - 6*tau_i(1)*(5 + 8*Dx(1))))/240.);

    J[20] = current_density(1) *
        (-(tau_i(0)*(-4 + 3*tau_i(0))*Dx(0)*(6*es - 8*e*Dx(1) + 3*dys))/48 +
         ((-2 + 3*tau_i(0))*dxs*(10*es - 15*e*Dx(1) + 6*dys))/120. -
         (dxt*(15*es - 24*e*Dx(1) + 10*dys))/240. +
         ((4 + 3*(-2 + tau_i(0))*as)*(3 + 3*bs - 3*tau_i(1)*(2 + Dx(1)) + Dx(1)*(3 + Dx(1))))/36.);

    J[21] = current_density(1) *
        (-(ct*(3*bs - 3*tau_i(1)*Dx(1) + dys))/36. +
         (cs*Dx(0)*(6*bs - 8*tau_i(1)*Dx(1) + 3*dys))/48. -
         (c*dxs*(10*bs - 15*tau_i(1)*Dx(1) + 6*dys))/120. +
         (dxt*(15*bs - 24*tau_i(1)*Dx(1) + 10*dys))/720.);

    J[22] = current_density(1) *
        (-(cs*Dx(0)*(-3 + 6*e*tau_i(1) + 4*Dx(1) - 8*tau_i(1)*Dx(1) + 3*dys))/24. +
         (ct*(-3 + 6*bs - 6*tau_i(1)*(1 + Dx(1)) + Dx(1)*(3 + 2*Dx(1))))/36. +
         (c*dxs*(-10 + 20*bs - 10*tau_i(1)*(2 + 3*Dx(1)) + 3*Dx(1)*(5 + 4*Dx(1))))/120. +
         (dxt*(15 - 30*bs - 4*Dx(1)*(6 + 5*Dx(1)) + 6*tau_i(1)*(5 + 8*Dx(1))))/720.);

    J[23] = current_density(1) *
        ((cs*Dx(0)*(6*es - 8*e*Dx(1) + 3*dys))/48. -
         (c*dxs*(10*es - 15*e*Dx(1) + 6*dys))/120. +
         (dxt*(15*es - 24*e*Dx(1) + 10*dys))/720. -
         (ct*(3 + 3*bs - 3*tau_i(1)*(2 + Dx(1)) + Dx(1)*(3 + Dx(1))))/36.);
}

/***************************************************************************
 * $RCSfile: edge.cpp,v $   $Author: chkraus $
 * $Revision: 1.1.1.1 $   $Date: 2013/08/23 07:40:38 $
 * IPPL_VERSION_ID: $Id: edge.cpp,v 1.1.1.1 2013/08/23 07:40:38 chkraus Exp $
 ***************************************************************************/
