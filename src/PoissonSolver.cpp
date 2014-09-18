/***************************************************************************
                          PoissonSolver.cpp
                         -------------------
    begin                : Tue Jun 28 2011
    copyright            : (C) 2011 by Christof Kraus
    email                : christof.kraus-csrst@my.mail.de
***************************************************************************/

#include "PoissonSolver.hh"
#include "utils.hh"

#include "Epetra_VbrMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Operator.h"
#include "Epetra_Time.h"

#include "EpetraExt_RowMatrixOut.h"

#include "Teuchos_Array.hpp"
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockCGSolMgr.hpp"

#include "Ifpack.h"

#include "Physics.hh"

extern ostream dbg;
#define DBGOUT dbg << "PoissonSolver.cpp: " << __LINE__ << "\t"

struct IndexComp {
    bool operator()(const Index & a,
                    const Index & b) {
        return (a.first() < b.first());
    }
};

typedef double                          ST;
typedef Epetra_Operator                 OP;
typedef Epetra_MultiVector              MV;
typedef Belos::OperatorTraits<ST,MV,OP> OPT;
typedef Belos::MultiVecTraits<ST,MV>    MVT;

PoissonSolver::PoissonSolver(PartBunch & bunch,
                             const Bend & bend,
                             const std::vector<BoundaryCell>& bc,
                             const Mesh_t & mesh,
                             const FieldLayout<DIM> & fl,
                             const BCType & bctX,
                             const BCType & bctY):
    _gamma(0.0),
    _Nx(0),
    _Ny(0),
    _hr(DIM),
    _bctX(bctX),
    _bctY(bctY),
    _bunch(bunch),
    _bend(bend),
    _boundaryCells(bc),
    _mesh(mesh),
    _FL(fl)
{
    Vector_t pmean;
    _bunch.get_pmean(pmean);
    _gamma = sqrt(dot(pmean, pmean) + 1);

    _lDom = _FL.getLocalNDIndex();
    _gDom = _FL.getDomain();

    _hr[0] = _mesh.get_meshSpacing(0) * _gamma;
    _hr[1] = _mesh.get_meshSpacing(1);
#if DIM>2
    _hr[2] = _mesh.get_meshSpacing(2);
#endif

    _Nx = _gDom[0].length();
    _Ny = _gDom[1].length();

    _setupTimer = Timings::getTimer("setup poisson");
    _solveTimer = Timings::getTimer("solve poisson");
}

/// actual computation
void PoissonSolver::computeField(VField_Edge_t & EFD,
                                 VField_Cell_t & HFD,
                                 VField_Edge_t & JFD,
                                 const double & dt,
                                 const double & tol,
                                 const int & maxIterations)
{
    Mesh_t & mesh = EFD.get_mesh();
    FieldLayout<DIM> & FL = EFD.getLayout();
    const NDIndex<DIM> lDom = FL.getLocalNDIndex();

    _bunch.saveR();

    const GCS_t & gcs = EFD.getGuardCellSizes();

    Timings::startTimer(_setupTimer);
    Epetra_Map* Map = getGlobalElements();

    Teuchos::RCP<Epetra_Vector> LHS;
    LHS = Teuchos::rcp(new Epetra_Vector(*Map));

    Teuchos::RCP<Epetra_Vector> RHS;
    RHS = Teuchos::rcp(new Epetra_Vector(*Map));

    Teuchos::RCP<Epetra_CrsMatrix> A;
    A = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *Map, 5));

    StencilGeometry(RHS, A);

    // print("PoissonMatrix.dat", RHS, A);

    Teuchos::ParameterList belosList;
    belosList.set( "Maximum Iterations", maxIterations );  // Maximum number of iterations allowed
    belosList.set( "Convergence Tolerance", tol );
    belosList.set( "Verbosity", (Belos::Errors +
                                 Belos::Warnings +
                                 Belos::TimingDetails +
                                 Belos::FinalSummary +
                                 Belos::StatusTestDetails) );

    Teuchos::ParameterList MLList;
    SetupMLList(MLList);

    Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> MLPrec =
        Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*A, MLList,false));

    MLPrec->ComputePreconditioner();
    Teuchos::RCP<Belos::EpetraPrecOp> prec = Teuchos::rcp(new Belos::EpetraPrecOp(MLPrec));

    Belos::LinearProblem<ST, MV, OP> problem;
    problem.setOperator(A);
    problem.setLHS(LHS);
    problem.setRHS(RHS);
    problem.setLeftPrec(prec);

    if (!problem.isProblemSet()) {
        if (!problem.setProblem()) {
            std::cerr << "\nERROR: Belos::LinearProblem failed to set up correctly!" << std::endl;
        }
    }

    Teuchos::RCP< Belos::SolverManager<ST, MV, OP> > solver;
    solver = Teuchos::rcp( new Belos::BlockCGSolMgr<ST, MV, OP>(Teuchos::rcp(&problem, false), Teuchos::rcp(&belosList, false)));
    Timings::stopTimer(_setupTimer);

    BinaryVtkFile vtkFile;
    SField_t rho(mesh, FL, gcs);

    _bunch.drift_particles(dt / 2);
    _bunch.scatterQ(rho);
    _bunch.drift_particles(-dt / 2);

    fillRHS(rho, RHS);
    LHS->Random();
    problem.setProblem(Teuchos::null, RHS);
    Timings::startTimer(_solveTimer);
    solver->solve();
    Timings::stopTimer(_solveTimer);
    plotPotential(vtkFile, LHS, EFD, 1);
    fillFields(LHS, EFD, HFD, 1);

    _bunch.move_by(Vector_t(-0.5 * _hr[0] / _gamma, 0.0));
    _bunch.scatterQ(rho);

    fillRHS(rho, RHS);
    shiftLHS(rho, LHS, 1.0);
    problem.setProblem(Teuchos::null, RHS);
    belosList.set("Convergence Tolerance", 1e-6);
    solver->setParameters(Teuchos::rcp(&belosList, false));
    Timings::startTimer(_solveTimer);
    solver->solve();
    Timings::stopTimer(_solveTimer);
    plotPotential(vtkFile, LHS, EFD, 2);
    fillFields(LHS, EFD, HFD, 2);

    _bunch.restoreR();

    _bunch.scatterQ(rho);

    vtkFile.addScalarField(rho, "rho");

    JFD = 0.0;
    JFD[lDom[0]][lDom[1]](0) += Physics::c * sqrt(1.0 - 1.0 / (_gamma * _gamma)) * rho[lDom[0]][lDom[1]];

    fillRHS(rho, RHS);
    shiftLHS(rho, LHS, -0.5);
    problem.setProblem(Teuchos::null, RHS);
    Timings::startTimer(_solveTimer);
    solver->solve();
    Timings::stopTimer(_solveTimer);
    plotPotential(vtkFile, LHS, EFD, 0);
    fillFields(LHS, EFD, HFD, 0);

    vtkFile.writeFile("Data/potential");

    delete Map;
}

/// building the stencil
void PoissonSolver::StencilGeometry(Teuchos::RCP<Epetra_Vector> & RHS,
                                    Teuchos::RCP<Epetra_CrsMatrix> & A)
{
    const Epetra_BlockMap & MyMap = RHS->Map();
    int NumMyElements = MyMap.NumMyElements();
    int* MyGlobalElements = MyMap.MyGlobalElements();
    double * rhsvalues = RHS->Values();

    std::vector<double> Values(5);
    std::vector<int> Indices(5);

    const auto & inside = _bend.getInsideMask();

    for (int lid = 0; lid < NumMyElements; ++ lid) {
        size_t NumEntries = 0;

        const size_t & gid = MyGlobalElements[lid];

        cutoffStencil(Indices,
                      Values,
                      rhsvalues[lid],
                      NumEntries,
                      inside,
                      gid);
        A->InsertGlobalValues(gid, NumEntries, &Values[0], &Indices[0]);
    }

    A->FillComplete();
    A->OptimizeStorage();
}

void PoissonSolver::initializeLHS(Teuchos::RCP<Epetra_Vector> & LHS) const
{
    const Epetra_BlockMap & Map = LHS->Map();
    const int * MyGlobalElements = Map.MyGlobalElements();
    const int NumMyElements = Map.NumMyElements();

    ST * values = LHS->Values();

    // size_t length = _bend.getLengthStraightSection() + 1;
    // size_t width = _bend.getWidthStraightSection() + 1;
    Vector_t position;
    _bunch.get_rmean(position);
    _bend.getPositionInCells(position);

    int lastX = (_bend.getGlobalDomain())[0].last();
    int lastY = (_bend.getGlobalDomain())[1].last();

    // NDIndex<DIM> initDom(Index(0,length), Index(lastY - width, lastY));

    NDIndex<DIM> elem;
    for (int lid = 0; lid < NumMyElements; ++ lid) {
        const size_t idx = MyGlobalElements[lid] % _Nx;
        const size_t idy = MyGlobalElements[lid] / _Nx;
        elem[0] = Index(idx, idx);
        elem[1] = Index(idy, idy);

        Vector_t contr;
        if (idx > position(0)) {
            contr(0) = (idx - position(0)) / (lastX - position(0));
        } else {
            contr(0) = (position(0) - idx) / position(0);
        }
        if (idy > position(1)) {
            contr(1) = (idy - position(1)) / (lastY - position(1));
        } else {
            contr(1) = (position(1) - idy) / position(1);
        }

        const double val = sqrt(dot(contr, contr));
        dbg << std::min(10.0, -log(val)) << std::endl;
        values[lid] = -std::max(0.0, std::min(10.0,-log(val)));
    }
}

void PoissonSolver::shiftLHS(SField_t & rho,
                             Teuchos::RCP<Epetra_Vector> & LHS,
                             double tau) const
{
    const Epetra_BlockMap & Map = LHS->Map();
    const int * MyGlobalElements = Map.MyGlobalElements();
    const int NumMyElements = Map.NumMyElements();

    NDIndex<DIM> elem;
    NDIndex<DIM> ldom = rho.getLayout().getLocalNDIndex();
    Index II = ldom[0], JJ = ldom[1];
    ST * values = LHS->Values();

    for (int lid = 0; lid < NumMyElements; ++ lid) {
        const size_t idx = MyGlobalElements[lid] % _Nx;
        const size_t idy = MyGlobalElements[lid] / _Nx;
        elem[0] = Index(idx, idx);
        elem[1] = Index(idy, idy);

        rho.localElement(elem) = values[lid];
    }
    rho.fillGuardCells();

    if (tau > 0) {
        if (tau > 1) tau = 1;
        rho[II][JJ] += tau * rho[II+1][JJ] - tau * rho[II][JJ];
    } else {
        tau = - tau;
        if (tau > 1) tau = 1;
        rho[II][JJ] += tau * rho[II-1][JJ] - tau * rho[II][JJ];
    }

    for (int lid = 0; lid < NumMyElements; ++ lid) {
        const int idx = MyGlobalElements[lid] % _Nx;
        const int idy = MyGlobalElements[lid] / _Nx;
        elem[0] = Index(idx, idx);
        elem[1] = Index(idy, idy);

        values[lid] = rho.localElement(elem);
    }
}

void PoissonSolver::subtract(const SField_t & rho,
                             Teuchos::RCP<Epetra_Vector> & LHS) const
{
    const Epetra_BlockMap & Map = LHS->Map();
    const int * MyGlobalElements = Map.MyGlobalElements();
    const int NumMyElements = Map.NumMyElements();

    NDIndex<DIM> elem;
    NDIndex<DIM> ldom = rho.getLayout().getLocalNDIndex();
    Index II = ldom[0], JJ = ldom[1];
    ST * values = LHS->Values();

    for (int lid = 0; lid < NumMyElements; ++ lid) {
        const size_t idx = MyGlobalElements[lid] % _Nx;
        const size_t idy = MyGlobalElements[lid] / _Nx;
        elem[0] = Index(idx, idx);
        elem[1] = Index(idy, idy);

        values[lid] -= rho.localElement(elem);
    }
}

void PoissonSolver::add(const SField_t & rho,
                        Teuchos::RCP<Epetra_Vector> & LHS) const
{
    const Epetra_BlockMap & Map = LHS->Map();
    const int * MyGlobalElements = Map.MyGlobalElements();
    const int NumMyElements = Map.NumMyElements();

    NDIndex<DIM> elem;
    NDIndex<DIM> ldom = rho.getLayout().getLocalNDIndex();
    Index II = ldom[0], JJ = ldom[1];
    ST * values = LHS->Values();

    for (int lid = 0; lid < NumMyElements; ++ lid) {
        const size_t idx = MyGlobalElements[lid] % _Nx;
        const size_t idy = MyGlobalElements[lid] / _Nx;
        elem[0] = Index(idx, idx);
        elem[1] = Index(idy, idy);

        values[lid] += rho.localElement(elem);
    }
}

void PoissonSolver::initialGuessLHS(Teuchos::RCP<Epetra_Vector> & LHS) const
{
    const Epetra_BlockMap & Map = LHS->Map();
    const int * MyGlobalElements = Map.MyGlobalElements();
    const int NumMyElements = Map.NumMyElements();
    Vector_t spos;
    NDIndex<DIM> elem;
    Vector_t origin = _mesh.get_origin();
    double totalQ = -_bunch.get_qtotal() / (2 * Physics::pi * Physics::epsilon_0);
    totalQ /= 97.63;
    double *values = LHS->Values();

    _bunch.get_rmean(spos);
    spos[0] = (spos[0] - origin(0)) / _mesh.get_meshSpacing(0);
    spos[1] = (spos[1] - origin(1)) / _mesh.get_meshSpacing(1);

    for (int lid = 0; lid < NumMyElements; ++ lid) {
        const double dx = _hr[0] * (spos[0] - MyGlobalElements[lid] % _Nx);
        const double dy = _hr[1] * (spos[1] - MyGlobalElements[lid] / _Nx);
        values[lid] = totalQ / sqrt(dx*dx + dy*dy + _hr[0]*_hr[1]);
    }
}

void PoissonSolver::fillRHS(const SField_t & rho,
                            Teuchos::RCP<Epetra_Vector> & RHS)
{
    const ST couplingConstant = 1 / (Physics::epsilon_0 * _gamma);
    const Epetra_BlockMap & Map = RHS->Map();
    const int * MyGlobalElements = Map.MyGlobalElements();
    const int NumMyElements = Map.NumMyElements();

    NDIndex<DIM> elem;
    ST * values = RHS->Values();
    // double minv = 0.0;
    // double maxv = 0.0;
    for (int lid = 0; lid < NumMyElements; ++ lid) {
        const size_t idx = MyGlobalElements[lid] % _Nx;
        const size_t idy = MyGlobalElements[lid] / _Nx;
        elem[0] = Index(idx, idx);
        elem[1] = Index(idy, idy);

        values[lid] = rho.localElement(elem) * couplingConstant;
        // if (minv > values[lid]) minv = values[lid];
        // if (maxv < values[lid]) maxv = values[lid];
    }
//     dbg << "rhs_{min}: " << minv << " - rhs_{max}: " << maxv << endl;
}

void PoissonSolver::plotPotential(BinaryVtkFile & vtkFile,
                                  const Teuchos::RCP<Epetra_Vector> & LHS,
                                  const VField_Edge_t & EFD,
                                  const int & component)
{
    const Epetra_BlockMap & Map = LHS->Map();
    const int * MyGlobalElements = Map.MyGlobalElements();
    const int NumMyElements = Map.NumMyElements();

    boost::shared_ptr<SField_Vert_t> scalarField = Utils::getScalarVertField(EFD);
    NDIndex<DIM> elem;
    char fieldName[50];
    ST * values = LHS->Values();

    sprintf(fieldName, "potential_%d", component);

    for (int lid = 0; lid < NumMyElements; ++ lid) {
        const size_t idx = MyGlobalElements[lid] % _Nx;
        const size_t idy = MyGlobalElements[lid] / _Nx;
        elem[0] = Index(idx, idx);
        elem[1] = Index(idy, idy);

        scalarField->localElement(elem) = values[lid];
    }

    vtkFile.addScalarField(*scalarField, fieldName);
}

void PoissonSolver::fillFields(const Teuchos::RCP<Epetra_Vector> & LHS,
                               VField_Edge_t & EFD,
                               VField_Cell_t & HFD,
                               const int & components)
{
    const Epetra_BlockMap & Map = LHS->Map();
    const int * MyGlobalElements = Map.MyGlobalElements();
    const int NumMyElements = Map.NumMyElements();

    const Vector_t hr(_mesh.get_meshSpacing(0),
                      _mesh.get_meshSpacing(1));

    Index II = _lDom[0], JJ = _lDom[1];
    NDIndex<DIM> elem, elem2;
    ST * values = LHS->Values();

    switch (components) {
    case 0:
        for (int lid = 0; lid < NumMyElements; ++ lid) {
            const int idx = MyGlobalElements[lid] % _Nx;
            const int idy = MyGlobalElements[lid] / _Nx;
            elem[0] = Index(idx, idx);
            elem[1] = Index(idy, idy);

            EFD.localElement(elem)(0) = values[lid] / _hr[0];
        }

        EFD.fillGuardCells();

        EFD[II][JJ](0) -= EFD[II+1][JJ](0);
        break;
    case 1:
        for (int lid = 0; lid < NumMyElements; ++ lid) {
            const int idx = MyGlobalElements[lid] % _Nx;
            const int idy = MyGlobalElements[lid] / _Nx;
            elem[0] = Index(idx, idx);
            elem[1] = Index(idy, idy);

            EFD.localElement(elem)(1) = _gamma * values[lid] / _hr[1];
        }

        EFD.fillGuardCells();

        EFD[II][JJ](1) -= EFD[II][JJ+1](1);
        break;
    case 2:
        {
            double betaGamma = sqrt(_gamma * _gamma - 1);
            for (int lid = 0; lid < NumMyElements; ++ lid) {
                const int idx = MyGlobalElements[lid] % _Nx;
                const int idy = MyGlobalElements[lid] / _Nx;
                elem[0] = Index(idx, idx);
                elem[1] = Index(idy, idy);
                double hfd =  betaGamma * values[lid] / (2 * Physics::c * Physics::mu_0 * _hr[1]);
                HFD.localElement(elem) = Vector_t(hfd, hfd);
            }
            HFD.fillGuardCells();

            HFD[II][JJ](0) -= HFD[II][JJ+1](0);
            HFD[II][JJ](1) -= HFD[II][JJ+1](1);
            break;
        }
    case -1:
    default:
        ;
    }
}

bool PoissonSolver::isInside(const int & idx, const int & idy)
{
    if (idx < 0 || idx > _gDom[0].last() ||
        idy < 0 || idy > _gDom[1].last()) {
        return false;
    }
    return true;
}

void PoissonSolver::distance_bunch_to_boundary(std::vector<double> & dist,
                                               const double & x,
                                               const double & y)
{
    Vector_t spos;
    _bunch.get_rmean(spos);
    dist[0] = std::abs(x - spos[0]);
    dist[1] = std::abs(y - spos[1]);
}

void PoissonSolver::print(const std::string & filename,
                          Teuchos::RCP<Epetra_Vector> & RHS,
                          Teuchos::RCP<Epetra_CrsMatrix> & A) const
{
    EpetraExt::RowMatrixToMatlabFile(filename.c_str(), *A);
//     if (Ippl::getNodes() == 1) {
//         const Epetra_BlockMap & Map = RHS->Map();
//         const int NumMyElements = Map.NumMyElements();

//         ST * values = RHS->Values();
//         ofstream out("rhs.dat");
//         for (size_t lid = 0; lid < NumMyElements; ++ lid) {
//             out << values[lid] << "\n";
//         }
//         out.close();
//     }
}

void PoissonSolver::linearStencil(std::vector<int> & Indices,
                                  std::vector<double> & Values,
                                  double & rhs,
                                  size_t & NumEntries,
                                  const Field<bool, DIM, Mesh_t, Vert> & isInside,
                                  const size_t & gid)
{
    const double EW = -1 / (_hr[0]*_hr[0]);
    const double NS = -1 / (_hr[1]*_hr[1]);

    const int idx = gid % _Nx;
    const int idy = gid / _Nx;
    const NDIndex<DIM> elem(Index(idx, idx), Index(idy, idy));

    if (isInside.localElement(elem)) {
        NumEntries = 5;

        Values[0] = NS;
        Indices[0] = gid - _Nx;
        Values[1] = EW;
        Indices[1] = gid - 1;
        Values[2] = -2 * (EW + NS);
        Indices[2] = gid;
        Values[3] = EW;
        Indices[3] = gid + 1;
        Values[4] = NS;
        Indices[4] = gid + _Nx;

        return;
    }

    const std::vector<int> & localToBCellNr = _bend.getInverseMapDualGrid();
    const int & localID = localToBCellNr[_bend.getLocalID(elem)];

    if (localID == -1) {
        NumEntries = 1;
        Values[0] = 1;
        Indices[0] = gid;

        NDIndex<DIM> elemmy(elem[0], elem[1] - 1);
        int localIDmy = localToBCellNr[_bend.getLocalID(elemmy)];
        if (localIDmy != -1 && _boundaryCells[localIDmy].lambda_m(1) > SPACE_EPS) {
            const BoundaryCell & bc = _boundaryCells[localIDmy];
            const double & s = bc.lambda_m(1);

            Values[NumEntries] = (1/s - 1);
            Indices[NumEntries] = gid - _Nx;

            dbg << gid << "    "
                << gid - _Nx << "    "
                << Values[NumEntries] << std::endl;

            Values[0] = NumEntries;
            NumEntries += 1;
        }
        NDIndex<DIM> elemmx(elem[0] - 1, elem[1]);
        int localIDmx = localToBCellNr[_bend.getLocalID(elemmx)];
        if (localIDmx != -1 && _boundaryCells[localIDmx].lambda_m(0) > SPACE_EPS) {
            const BoundaryCell & bc = _boundaryCells[localIDmx];
            const double & s = bc.lambda_m(0);

            Values[NumEntries] = (1/s - 1);
            Indices[NumEntries] = gid - 1;

            dbg << gid << "    "
                << gid - 1 << "    "
                << Values[NumEntries] << std::endl;

            Values[0] = NumEntries;
            NumEntries += 1;
        }

        rhs = 0.0;
        return;
    }

    const BoundaryCell & bc = _boundaryCells[localID];
    if (bc.lambda_m[0] > 1 - SPACE_EPS ||
        bc.lambda_m[1] > 1 - SPACE_EPS) {

        NumEntries = 5;

        Values[0] = NS;
        Indices[0] = gid - _Nx;
        Values[1] = EW;
        Indices[1] = gid - 1;
        Values[2] = -2 * (EW + NS);
        Indices[2] = gid;
        Values[3] = EW;
        Indices[3] = gid + 1;
        Values[4] = NS;
        Indices[4] = gid + _Nx;

    } else {
        NumEntries = 1;
        Values[0] = 1;
        Indices[0] = gid;

        if (bc.lambda_m(0) > SPACE_EPS) {
            const double & s = bc.lambda_m(0);

            Values[NumEntries] = (1/s - 1);
            Indices[NumEntries] = gid + 1;

            dbg << gid << "    "
                << gid + 1 << "    "
                << Values[NumEntries] << std::endl;

            Values[0] = NumEntries;
            NumEntries += 1;
        }

        if (bc.lambda_m(1) > SPACE_EPS) {
            const double & s = bc.lambda_m(1);

            Values[NumEntries] = (1/s - 1);
            Indices[NumEntries] = gid + _Nx;

            dbg << gid << "    "
                << gid + _Nx << "    "
                << Values[NumEntries] << std::endl;

            Values[0] = NumEntries;
            NumEntries += 1;
        }

        rhs = 0.0;
    }

    return;
}

void PoissonSolver::cutoffStencil(std::vector<int> & Indices,
                                  std::vector<double> & Values,
                                  double & rhs,
                                  size_t & NumEntries,
                                  const Field<bool, DIM, Mesh_t, Cell> & isInside,
                                  const size_t & gid)
{
    const double EW = -1 / (_hr[0]*_hr[0]);
    const double NS = -1 / (_hr[1]*_hr[1]);

    const int idx = gid % _Nx;
    const int idy = gid / _Nx;
    const NDIndex<DIM> elem(Index(idx, idx), Index(idy, idy));
    const NDIndex<DIM> elemmx(elem[0] - 1, elem[1]);
    const NDIndex<DIM> elemmy(elem[0], elem[1] - 1);
    const NDIndex<DIM> elempx(elem[0] + 1, elem[1]);
    const NDIndex<DIM> elempy(elem[0], elem[1] + 1);

    const std::vector<int> & localToBCellNr = _bend.getInverseMap();
    const int & localID = localToBCellNr[_bend.getLocalID(elem)];
    const int & localIDmx = localToBCellNr[_bend.getLocalID(elemmx)];
    const int & localIDmy = localToBCellNr[_bend.getLocalID(elemmy)];
    const int & localIDpx = localToBCellNr[_bend.getLocalID(elempx)];
    const int & localIDpy = localToBCellNr[_bend.getLocalID(elempy)];

    long maxX = _bend.getLengthStraightSection();

    if (idx < maxX &&
        (isInside.localElement(elem) ||
         (localID > -1 &&
          _boundaryCells[localID].area_m > 0.5))) {

        NumEntries = 0;

        if (isInside.localElement(elemmy) ||
            (localIDmy > -1 &&
             _boundaryCells[localIDmy].area_m > 0.5)) {

            Indices[NumEntries] = gid - _Nx;
            Values[NumEntries] = NS;

            NumEntries ++;
        }

        if (isInside.localElement(elemmx) ||
            (localIDmx > -1 &&
             _boundaryCells[localIDmx].area_m > 0.5)) {

            Indices[NumEntries] = gid - 1;
            Values[NumEntries] = EW;

            NumEntries ++;
        }

        if (isInside.localElement(elempx) ||
            (localIDpx > -1 &&
             _boundaryCells[localIDpx].area_m > 0.5)) {

            Indices[NumEntries] = gid + 1;
            Values[NumEntries] = EW;

            NumEntries ++;
        }

        if (isInside.localElement(elempy) ||
            (localIDpy > -1 &&
             _boundaryCells[localIDpy].area_m > 0.5)) {

            Indices[NumEntries] = gid + _Nx;
            Values[NumEntries] = NS;

            NumEntries ++;
        }

        Indices[NumEntries] = gid;
        Values[NumEntries] = -2 * (EW + NS);

        NumEntries ++;

        return;
    }

    NumEntries = 1;
    Indices[0] = gid;
    Values[0] = 1;
    rhs = 0;

    return;
}
