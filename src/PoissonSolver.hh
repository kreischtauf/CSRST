/***************************************************************************
                           PoissonSolver.hh
                         -------------------
    begin                : Tue Jun 28 2011
    copyright            : (C) 2011 by Christof Kraus
    email                : christof.kraus-csrst@my.mail.de
***************************************************************************/

#ifndef POISSONSOLVER_HH
#define POISSONSOLVER_HH

#define HAVE_MPI

#include <vector>

#ifdef HAVE_MPI
#    include "mpi.h"
#    include "Epetra_MpiComm.h"
#else
#    include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "ml_MultiLevelPreconditioner.h"
#include "ml_MultiLevelOperator.h"
#include "ml_epetra_utils.h"

#include "defs.hh"
#include "PartBunch.hh"
#include "Bend.hh"
#include "BoundaryCell.hh"
#include "Commands/BinaryVTKFile.hh"
#include "Commands/ASCIIVTKFile.hh"

/** \brief computes the electrostatic field due to a charge density
 *         inside the domain. Used at the beginning of a simulation
 *         where only the particles distribution in phase space is
 *         known. Implemented as described in
 *         <a href="http://dx.doi.org/10.1016/j.jcp.2010.02.022">A fast parallel Poisson solver on irregular domains applied to beam dynamics simulations</a>
 *         by A. Adelmann,P. Arbenz and Y. Ineichen.
 */
class PoissonSolver
{
public:
    PoissonSolver(PartBunch & bunch,
                  const Bend & bend,
                  const std::vector<BoundaryCell>& bc,
                  const Mesh_t & mesh,
                  const FieldLayout<DIM> & fl,
                  const BCType & periodicX = PEC,
                  const BCType & periodicY = PEC);
    ~PoissonSolver()
    { ;}

    /// actual computation
    void computeField(VField_Edge_t & EFD,
                      VField_Cell_t & HFD,
                      VField_Edge_t & JFD,
                      const double & dt,
                      const double & tol = 1e-8,
                      const int & maxIterations = 100);
private:
    void SetupMLList(Teuchos::ParameterList & MLList) const;

    Epetra_Map* getGlobalElements();

    /// building the stencil
    void StencilGeometry(Teuchos::RCP<Epetra_Vector> & RHS,
                         Teuchos::RCP<Epetra_CrsMatrix> & A);

    void linearStencil(std::vector<int> & Indices,
                       std::vector<double> & Values,
                       double & rhs,
                       size_t & NumEntries,
                       const Field<bool, DIM, Mesh_t, Vert> & isInside,
                       const size_t & gid);

    void cutoffStencil(std::vector<int> & Indices,
                       std::vector<double> & Values,
                       double & rhs,
                       size_t & NumEntries,
                       const Field<bool, DIM, Mesh_t, Cell> & isInside,
                       const size_t & gid);

    void initializeLHS(Teuchos::RCP<Epetra_Vector> & LHS) const;

    void shiftLHS(SField_t & rho,
                  Teuchos::RCP<Epetra_Vector> & LHS,
                  double tau = 1.0) const;

    void subtract(const SField_t & rho,
                  Teuchos::RCP<Epetra_Vector> & LHS) const;

    void add(const SField_t & rho,
             Teuchos::RCP<Epetra_Vector> & LHS) const;

    void initialGuessLHS(Teuchos::RCP<Epetra_Vector> & LHS) const;

    void fillRHS(const SField_t & rho,
                 Teuchos::RCP<Epetra_Vector> & RHS);

    void fillFields(const Teuchos::RCP<Epetra_Vector> & RHS,
                    VField_Edge_t & EFD,
                    VField_Cell_t & HFD,
                    const int & component = -1);

    void plotPotential(BinaryVtkFile & vtkFile,
                       const Teuchos::RCP<Epetra_Vector> & LHS,
                       const VField_Edge_t & EFD,
                       const int & component);

    bool isInside(const int & idx, const int & idy);

    void distance_bunch_to_boundary(std::vector<double> & dist,
                                    const double & x,
                                    const double & y);

    void print(const std::string & filename,
               Teuchos::RCP<Epetra_Vector> & RHS,
               Teuchos::RCP<Epetra_CrsMatrix> & A) const;

    void printRhoOverLine(const Teuchos::RCP<Epetra_Vector> & LHS,
                          const int & line_nr,
                          std::ostream & out);

    void printEyOverLine(const VField_Edge_t & EFD,
                         const int & line_nr,
                         std::ostream & out);

    void printHzOverLine(const VField_Cell_t & HFD,
                         const int & line_nr,
                         std::ostream & out);

    double _gamma;

    NDIndex<DIM> _lDom;
    NDIndex<DIM> _gDom;
    size_t _Nx;
    size_t _Ny;
    std::vector<double> _hr;
    BCType _bctX;
    BCType _bctY;

    PartBunch & _bunch;

    const Bend & _bend;
    const std::vector<BoundaryCell> & _boundaryCells;
    const Mesh_t & _mesh;
    const FieldLayout<DIM> & _FL;

    Timings::TimerRef _solveTimer;
    Timings::TimerRef _setupTimer;


};

inline
void PoissonSolver::SetupMLList(Teuchos::ParameterList & MLList) const {
    ML_Epetra::SetDefaults("SA", MLList);
    MLList.set("max levels", 8);
    MLList.set("increasing or decreasing", "increasing");

    // we use a V-cycle
    MLList.set("prec type", "MGV");

    // uncoupled aggregation is used (every processor aggregates
    // only local data)
    MLList.set("aggregation: type", "Uncoupled");

    // smoother related parameters
    MLList.set("smoother: type","Chebyshev");
    MLList.set("smoother: sweeps", 3);
    MLList.set("smoother: pre or post", "both");

    // on the coarsest level we solve with  Tim Davis' implementation of
    // Gilbert-Peierl's left-looking sparse partial pivoting algorithm,
    // with Eisenstat & Liu's symmetric pruning. Gilbert's version appears
    // as \c [L,U,P]=lu(A) in MATLAB. It doesn't exploit dense matrix
    // kernels, but it is the only sparse LU factorization algorithm known to be
    // asymptotically optimal, in the sense that it takes time proportional to the
    // number of floating-point operations.
    MLList.set("coarse: type", "Amesos-KLU");

    //FIXME: CHEBY COARSE LEVEL SOLVER
    // SEE PAPER FOR EVALUATION KLU vs. Chebyshev
    //MLList.set("coarse: sweeps", 10);
    //MLList.set("coarse: type", "Chebyshev");

    // turn on all output
//     if(verbose_m)
//         MLList.set("ML output", 101);
//     else
        MLList.set("ML output", 10);

    // try to optimize mem for xt3
    //MLList.set("low memory usage", true);

    // heuristic for max coarse size depending on number of processors
    int coarsest_size = std::max(Ippl::getNodes() * 10, 1024);
    MLList.set("coarse: max size", coarsest_size);
}

inline
Epetra_Map* PoissonSolver::getGlobalElements()
{
#ifdef HAVE_MPI
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm Comm;
#endif

    int NumMyElements = (_lDom[0].last() - _lDom[0].first() + 1) *
                        (_lDom[1].last() - _lDom[1].first() + 1);

    unsigned int lid = 0;
    std::vector<int> MyGlobalElements(NumMyElements);
    for (int j = _lDom[1].first(); j <= _lDom[1].last(); ++ j) {
        for (int i = _lDom[0].first(); i <= _lDom[0].last(); ++ i) {
            MyGlobalElements[lid ++] = i + j * _Nx;
        }
    }

    return new Epetra_Map(_Nx * _Ny, NumMyElements, &MyGlobalElements[0], 0, Comm);
}

#endif
