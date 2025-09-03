#ifndef NDTHETASOLVER_HPP
#define NDTHETASOLVER_HPP

#include <fstream>
#include <iostream>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/lac/solver_minres.h>

#include "NDProblem.hpp"
#include "WhiteGrayPartition.hpp"


using namespace dealii;

template<unsigned int DIM>
class NDThetaSolver
{
public:
  // Constructor. We provide the final time, time step Delta t and theta method
  // parameter as constructor arguments.
    NDThetaSolver(
      NDProblem<DIM> &problem_,
      double theta_,
      double deltat_,
      double T_,
      unsigned int &r_,
      const std::string &output_directory_ = "./",
      const std::string &output_filename_ = "output"
    )
    :
      problem(problem_)
    , diffusion_tensor(problem_.get_diffusion_tensor())
    , theta(theta_)
    , deltat(deltat_)
    , T(T_)
    , mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , pcout(std::cout, mpi_rank == 0)
    , r(r_)
    , output_directory(output_directory_)
    , output_filename(output_filename_)
    , mesh(MPI_COMM_WORLD)
  {
    if(theta < 0.0 || theta > 1.0)
      throw std::runtime_error("Theta parameter must be in the interval [0, 1].");
  }

  virtual ~NDThetaSolver() = default;

  // Initialization.
  virtual void setup();

  // Solve the problem.
  virtual void solve();

protected:
  // Assemble the tangent problem.
  virtual void assemble_system();

  // Solve the linear system associated to the tangent problem.
  virtual void solve_linear_system();

  // Solve the problem for one time step using Newton's method.
  virtual void solve_newton();

  // Output.
  void output(const unsigned int &time_step) const;

  // Problem references
  const NDProblem<DIM> &problem;
  const typename NDProblem<DIM>::DiffusionTensor &diffusion_tensor;

  // Theta method parameter.
  const double theta;

  // Current time and time step.
  double time;
  double deltat;

  // Final time.
  double T;

  // MPI parallel. /////////////////////////////////////////////////////////////

  // Number of MPI processes.
  const unsigned int mpi_size;

  // This MPI process.
  const unsigned int mpi_rank;

  // Parallel output stream.
  ConditionalOStream pcout;

  // Polynomial degree.
  const unsigned int r;

  // directory where the output files will be written
  std::string output_directory;

  // output filename
  std::string output_filename;

  // Jacobian matrix.
  TrilinosWrappers::SparseMatrix jacobian_matrix;

  // Mesh.
  parallel::fullydistributed::Triangulation<DIM> mesh;
  Triangulation<DIM> mesh_serial;

  // Finite element space.
  std::unique_ptr<FiniteElement<DIM>> fe;

  // Quadrature formula.
  std::unique_ptr<Quadrature<DIM>> quadrature;

  // DoF handler.
  DoFHandler<DIM> dof_handler;

  // DoFs owned by current process.
  IndexSet locally_owned_dofs;

  // DoFs relevant to the current process (including ghost DoFs).
  IndexSet locally_relevant_dofs;

  // Residual vector.
  TrilinosWrappers::MPI::Vector residual_vector;

  // Increment of the solution between Newton iterations.
  TrilinosWrappers::MPI::Vector delta_owned;

  // System solution (without ghost elements).
  TrilinosWrappers::MPI::Vector solution_owned;

  // System solution (including ghost elements).
  TrilinosWrappers::MPI::Vector solution;

  // System solution at previous time step.
  TrilinosWrappers::MPI::Vector solution_old;

  // Global concentration.
  double global_concentration;

  // Total domain volume.
  double total_domain_volume;

  private: 

  // Write the fiber field to the output file.
  void write_fiber_field_to_file() const;

  // Integrate the current solution on whole domain to compute the global concentration.
  void compute_current_global_concentration();

};

// Cranck-Nicolson solver
template<unsigned int DIM>
class NDCrankNicolsonSolver : public NDThetaSolver<DIM>
{
  public:
    NDCrankNicolsonSolver(NDProblem<DIM> &problem_,
                          double deltat_,
                          double T_,
                          unsigned int &r_,
                          const std::string &output_directory_ = "./",
                          const std::string &output_filename_ = "output")
      : NDThetaSolver<DIM>(problem_, 0.5, deltat_, T_, r_, output_directory_, output_filename_)
    {}

    virtual ~NDCrankNicolsonSolver() = default;
};

// Backward Euler solver
template<unsigned int DIM>
class NDBackwardEulerSolver : public NDThetaSolver<DIM>
{
  public:
    NDBackwardEulerSolver(
      NDProblem<DIM> &problem_,
      double deltat_,
      double T_,
      unsigned int &r_,
      const std::string &output_directory_ = "./",
      const std::string &output_filename_ = "output"
    )
      : NDThetaSolver<DIM>(problem_, 1.0, deltat_, T_, r_, output_directory_, output_filename_)
    {}

    virtual ~NDBackwardEulerSolver() = default;
};



#endif 