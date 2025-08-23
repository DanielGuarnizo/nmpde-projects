#ifndef INITIAL_CONDITIONS_HPP
#define INITIAL_CONDITIONS_HPP

#include "NDProblem.hpp"
#include <deal.II/grid/grid_in.h> // For reading meshes
#include <deal.II/grid/grid_tools.h> // For finding the bounding box

// small value to be set as initial condition in the region where the initial condition is 0 to avoid the solution to go below zero
//#define EPSILON 5e-3
#define EPSILON 0

using namespace dealii;

/**
 * @brief Constant initial condition.
 * 
 * @note If ray = 0 (default value), the initial condition is constant on the whole domain.
 * 
 * @tparam DIM Dimension of the problem.
 */

template<unsigned int DIM>
class ConstantInitialCondition: public NDProblem<DIM>::InitialConcentration
{
    public:
        virtual double value(const Point<DIM> &p, const unsigned int /*component*/ = 0) const override
        {
            if(ray == 0.0)
              return C_0;
            if(p.distance(origin) > ray)
                // return 0.0;
                return EPSILON;
            return C_0;
        }
      
      ConstantInitialCondition(double C_0_, Point<DIM> origin_, double ray_)
        : C_0(C_0_), origin(origin_), ray(ray_) {}
        
    private:
      double C_0;
      Point <DIM> origin;
      double ray;
};

template<unsigned int DIM>
class ExponentialInitialCondition: public NDProblem<DIM>::InitialConcentration
{
    public:
        virtual double value(const Point<DIM> &p, const unsigned int /*component*/ = 0) const override
        {
            double distance_from_origin = p.distance(origin);
            if(distance_from_origin > ray)
              // return 0.0;
              return EPSILON;
            return C_0*std::exp(-distance_from_origin*distance_from_origin/(2*sigma*sigma)) + EPSILON;
        }
      
      ExponentialInitialCondition(Point<DIM> origin_ = Point<DIM>(), double sigma_ = 0.1, double C_0_ = 0.4, double ray_ = 4)
        : C_0(C_0_), origin(origin_), ray(ray_), sigma(sigma_) {}
        
    private:
      double C_0;
      Point<DIM> origin;
      double ray;
      double sigma;
};

template <unsigned int DIM>
class SmoothBumpInitialCondition : public NDProblem<DIM>::InitialConcentration {
public:
    virtual double value(const Point<DIM>& p, const unsigned int /*component*/ = 0) const override {
        double distance_from_origin_sq = p.distance(origin) * p.distance(origin);
        double r_sq = distance_from_origin_sq / (ray * ray); // Normalized squared radius

        if (r_sq >= 1.0) {
            return EPSILON; // Or 0.0 if your solver handles it and you want strict zero
        } else {
            // Standard bump function form: C_0 * exp(1 - 1/(1 - r^2))
            // The exp(1) can be absorbed into C_0 if desired for a slightly cleaner look
            // or use exp(-1 / (1 - r_sq)) and adjust C_0 accordingly
            double bump_val = C_0 * std::exp(1.0 - 1.0 / (1.0 - r_sq));
            return bump_val + EPSILON; // Still good practice for numerical stability
            // Or: return std::max(bump_val, EPSILON);
        }
    }

    SmoothBumpInitialCondition(Point<DIM> origin_ = Point<DIM>(),
                               double C_0_ = 0.4,
                               double ray_ = 4.0) 
        : C_0(C_0_), origin(origin_), ray(ray_) {
        }

private:
    double C_0;      // Amplitude
    Point<DIM> origin;
    double ray;      // Radius of the bump's support
    // static const double EPSILON = 1e-9; // Or define globally
};

template<unsigned int DIM>
class QuadraticInitialCondition: public NDProblem<DIM>::InitialConcentration
{
    public:
        virtual double value(const Point<DIM> &p, const unsigned int /*component*/ = 0) const override
        {
            double distance_from_origin_squared = p.distance_square(origin);
            if(distance_from_origin_squared > ray_squared)
              // return 0.0;
              return EPSILON;
            return C_0*(1 - distance_from_origin_squared/(ray_squared));
        }

        QuadraticInitialCondition(double C_0_, Point<DIM> origin_, double ray_)
        : C_0(C_0_), origin(origin_), ray_squared(ray_ * ray_) {}

    private:
        double C_0;
        Point<DIM> origin;
        double ray_squared;
};

// NEW CLASS: An initial condition that automatically finds the center of the mesh.
// template <int DIM>
// class CenteredExponentialInitialCondition : public NDProblem<DIM>::InitialConcentration
// {
//     public:
//     CenteredExponentialInitialCondition(const std::string &mesh_file_path)
//     {
//         Triangulation<DIM> temp_tria;
//         GridIn<DIM> grid_in;
//         grid_in.attach_triangulation(temp_tria);
//         std::ifstream in_file(mesh_file_path);

//         if (!in_file)
//         {
//             throw std::runtime_error("Could not open mesh file: " + mesh_file_path +
//                                     " to determine domain center.");
//         }
//         grid_in.read_msh(in_file);

//         // --- FIX #1 ---
//         // Changed get_bounding_box to compute_bounding_box
//         const auto bounding_box = GridTools::compute_bounding_box(temp_tria);
//         center = bounding_box.center();

//         std::cout << "Initial condition will be centered at x = " << center[0] << std::endl;
//     }

//     virtual double value(const Point<DIM> &p,
//                         const unsigned int /*component*/) const override
//     {
//         const double C0 = 0.1;
//         const double sigma = 0.05;

//         // --- FIX #2 ---
//         // Changed distance_squared to distance_square
//         const double distance_squared = center.distance_square(p);

//         return C0 * std::exp(-distance_squared / (2.0 * sigma * sigma));
//     }

//     private:
//     Point<DIM> center;
// };
template <int DIM>
class CenteredExponentialInitialCondition : public NDProblem<DIM>::InitialConcentration
{
public:
  // --- CHANGE 1: The constructor now accepts the C_0 value ---
  CenteredExponentialInitialCondition(const std::string &mesh_file_path, double c0_from_config)
    : c0_peak_value(c0_from_config) // Store the value in our member variable
  {
    Triangulation<DIM> temp_tria;
    GridIn<DIM> grid_in;
    grid_in.attach_triangulation(temp_tria);
    std::ifstream in_file(mesh_file_path);

    if (!in_file)
      {
        throw std::runtime_error("Could not open mesh file: " + mesh_file_path);
      }
    grid_in.read_msh(in_file);

    const auto bounding_box = GridTools::compute_bounding_box(temp_tria);
    center = bounding_box.center();

    std::cout << "Initial condition will be centered at x = " << center[0] << std::endl;
    std::cout << "Peak initial concentration C_0 set to: " << c0_peak_value << std::endl;
  }

  virtual double value(const Point<DIM> &p,
                       const unsigned int /*component*/) const override
  {
    const double sigma = 0.05; // Width of the seed

    const double distance_squared = center.distance_square(p);

    // --- CHANGE 2: Use the stored C_0 value instead of a hard-coded one ---
    return c0_peak_value * std::exp(-distance_squared / (2.0 * sigma * sigma));
  }

private:
  Point<DIM> center;
  const double c0_peak_value; // A member variable to hold the C_0 value
};

#endif // INITIAL_CONDITIONS_HPP