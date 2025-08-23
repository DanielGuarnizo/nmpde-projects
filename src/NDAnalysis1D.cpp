// #include "NDProblem.hpp"
// #include "InitialConditions.hpp"
// #include "FiberFields.hpp"
// #include "NDThetaSolver.hpp"
// #include "NDConfig.hpp"

// static NDConfig config_line = {
//     .dim = 1,
//     .T = 20.0,
//     .alpha = 1.8,
//     .deltat = 0.1,
//     .degree = 1,
//     .d_ext = 0.0,
//     .d_axn = 0.2,
//     .C_0 = 0.4,
//     .mesh = "../meshes/mesh-40.msh",
// };

// // static NDConfig config_line = {
// //     .dim = 1,
// //     .T = 7.0,
// //     .alpha = 1.0,
// //     .deltat = 0.01
// //     .degree = 1, // r
// //     .d_ext = 0.00,
// //     .d_axn = 0.2,
// //     .C_0 = 0.4,
// //     .mesh = "../meshes/mesh-40.msh",
// // };

// int main(int argc, char *argv[])
// {
//   Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

//   // NDConfig config = config_line;
//   NDConfig config;

//   config.parse(argc, argv);

//   ExponentialInitialCondition<1> initial_condition;
//   RadialFiberField<1> fiber_field;

//   NDProblem<1> problem(config.mesh, config.alpha, config.d_ext, config.d_axn, initial_condition, fiber_field);
//   problem.export_problem(std::string(config.output_dir) + config.output_filename + ".problem");

//   NDBackwardEulerSolver<1> solver(problem, config.deltat, config.T, config.degree, config.output_dir, config.output_filename);
//   solver.setup();
//   solver.solve();

//   return EXIT_SUCCESS;
// }


#include "NDProblem.hpp"
#include "InitialConditions.hpp"
#include "FiberFields.hpp"
#include "NDThetaSolver.hpp"
#include "NDConfig.hpp"

int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  // This is now a command-line driven program with no hard-coded defaults
  NDConfig config;
  config.parse(argc, argv);

  // --- THIS IS THE KEY CHANGE ---
  // Create our new "smart" initial condition.
  // It uses the mesh path provided via the command line to center itself.
  // CenteredExponentialInitialCondition<1> initial_condition(config.mesh);
  // in main()

  CenteredExponentialInitialCondition<1> initial_condition(config.mesh, config.C_0);
  // ExponentialInitialCondition<1> initial_condition;

  // The fiber field for a 1D problem is simple
  RadialFiberField<1> fiber_field;

  // Create the problem definition
  NDProblem<1> problem(config.mesh, config.alpha, config.d_ext, config.d_axn, initial_condition, fiber_field);
  problem.export_problem(std::string(config.output_dir) + config.output_filename + ".problem");

  // Create the numerical solver
  NDBackwardEulerSolver<1> solver(problem, config.deltat, config.T, config.degree, config.output_dir, config.output_filename);
  solver.setup();
  solver.solve();

  return EXIT_SUCCESS;
}