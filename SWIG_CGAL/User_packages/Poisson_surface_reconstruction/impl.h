#ifndef SWIG_CGAL_POISSON_SURFACE_RECONSTRUCTION_IMPL_H
#define SWIG_CGAL_POISSON_SURFACE_RECONSTRUCTION_IMPL_H

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Poisson_implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/IO/read_points.h>

#include <deque>
#include <cstdlib>
#include <fstream>
#include <math.h>

// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------

// kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

// Simple geometric types
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef std::pair<Point, Vector> Point_with_normal;
typedef Kernel::Sphere_3 Sphere;
typedef std::deque<Point_with_normal> PointList;

// Poisson implicit function
typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;

// Surface mesher
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
typedef CGAL::Poisson_implicit_surface_3<Kernel, Poisson_reconstruction_function> Surface_3;

struct Counter {
  std::size_t i, N;
  Counter(std::size_t N)
    : i(0), N(N)
  {}

  void operator()()
  {
    i++;
    if(i == N){
      std::cerr << "Counter reached " << N << std::endl;
    }
  }

};

struct InsertVisitor {

  Counter& c;
  InsertVisitor(Counter& c)
    : c(c)
  {}

  void before_insertion()
  {
    c();
  }

};

template <class Polyhedron, class PointSet>
int reconstruction_poly(PointSet &point_set, Polyhedron& out)
{
    // Poisson options
    FT sm_angle = 20.0; // Min triangle angle (degrees).
    FT sm_radius = 100; // Max triangle size w.r.t. point set average spacing.
    FT sm_distance = 0.25; // Approximation error w.r.t. point set average spacing.
    std::string solver_name = "eigen"; // Sparse linear solver name.
    double approximation_ratio = 0.02;
    double average_spacing_ratio = 5;

    //for (int i=3; i+1<argc ; ++i)
    //{
    //  if (std::string(argv[i])=="-sm_radius")
    //    sm_radius = atof(argv[++i]);
    //  else if (std::string(argv[i])=="-sm_distance")
    //    sm_distance = atof(argv[++i]);
    //  else if (std::string(argv[i])=="-solver")
    //    solver_name = argv[++i];
    //  else if (std::string(argv[i])=="-approx")
    //    approximation_ratio = atof(argv[++i]);
    //  else if (std::string(argv[i])=="-ratio")
    //    average_spacing_ratio = atof(argv[++i]);
    //  else {
    //    std::cerr << "Error: invalid option " << argv[i] << "\n";
    //    return EXIT_FAILURE;
    //  }
    //}

    //***************************************
    // Loads mesh/point set
    //***************************************
    PointList points;

    //// Reads the point set file into points.
    //// Note: read_points() requires an iterator over points
    //// + property maps to access each point's position and normal.
    //std::string input_filename("C:\\Users\\pvasi\\Code\\6kurs\\diplom\\python\\image\\image.xyz");
    //if (!CGAL::IO::read_points(input_filename.c_str(), std::back_inserter(points),
    //                            CGAL::parameters::point_map(CGAL::make_first_of_pair_property_map(Point_with_normal()))
    //                                .normal_map(CGAL::make_second_of_pair_property_map(Point_with_normal()))))
    //{
    //    std::cerr << "Error: cannot read input file!" << input_filename << std::endl;
    //    return EXIT_FAILURE;
    //}

    const auto &pts = point_set.points();
    const auto &nrmls = point_set.normals();

    auto p = pts.begin();
    auto n = nrmls.begin();
    for (; p != pts.end() && n != nrmls.end(); p++, n++)
        points.push_back(std::make_pair(*p, *n));

    std::size_t nb_points = points.size();
    std::cerr << "Input: " << nb_points << " points" << std::endl;

    //***************************************
    // Checks requirements
    //***************************************
    if (nb_points == 0)
    {
      std::cerr << "Error: empty point set" << std::endl;
      return EXIT_FAILURE;
    }

    bool points_have_normals = (points.begin()->second != CGAL::NULL_VECTOR);
    //bool points_have_normals = points.has_normal_map();
    if (!points_have_normals)
    {
      std::cerr << "Input point set not supported: this reconstruction method requires oriented normals" << std::endl;
      return EXIT_FAILURE;
    }

    Counter counter(std::distance(points.begin(), points.end()));
    InsertVisitor visitor(counter);

    //***************************************
    // Computes implicit function
    //***************************************

    std::cerr << "Computes Poisson implicit function...\n";

    // Creates implicit function from the read points.
    // Note: this method requires an iterator over points
    // + property maps to access each point's position and normal.
    Poisson_reconstruction_function function(
                              points.begin(), points.end(),
                              CGAL::make_first_of_pair_property_map(Point_with_normal()),
                              CGAL::make_second_of_pair_property_map(Point_with_normal()),
                              visitor);
    //Poisson_reconstruction_function function(
    //                          points.begin(), points.end(),
    //                          points.point_map(),
    //                          points.normal_map(),
    //                          visitor);

    #ifdef CGAL_EIGEN3_ENABLED
    {
      if (solver_name == "eigen")
      {
        std::cerr << "Use Eigen 3\n";
        CGAL::Eigen_solver_traits<Eigen::ConjugateGradient<CGAL::Eigen_sparse_symmetric_matrix<double>::EigenType> > solver;
        std::cerr << "Compute implicit function\n";
        try
        {
            if (!function.compute_implicit_function(solver, visitor,
                                                    approximation_ratio,
                                                    average_spacing_ratio))
            {
              std::cerr << "Error: cannot compute implicit function" << std::endl;
              return EXIT_FAILURE;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what();
            return EXIT_FAILURE;
        }
      }
      else
      {
        std::cerr << "Error: invalid solver " << solver_name << "\n";
        return EXIT_FAILURE;
      }
    }
    #else
    {
      std::cerr << "Error: invalid solver " << solver_name << "\n";
      return EXIT_FAILURE;
    }
    #endif

    //***************************************
    // Surface mesh generation
    //***************************************

    std::cerr << "Surface meshing...\n";

    // Computes average spacing
    FT average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>
      (points, 6 /* knn = 1 ring */,
       CGAL::parameters::point_map(CGAL::make_first_of_pair_property_map(Point_with_normal())));

    // Gets one point inside the implicit surface
    Point inner_point = function.get_inner_point();
    FT inner_point_value = function(inner_point);
    if(inner_point_value >= 0.0)
    {
      std::cerr << "Error: unable to seed (" << inner_point_value << " at inner_point)" << std::endl;
      return EXIT_FAILURE;
    }

    // Gets implicit function's radius
    Sphere bsphere = function.bounding_sphere();
    FT radius = std::sqrt(bsphere.squared_radius());

    // Defines the implicit surface: requires defining a
    // conservative bounding sphere centered at inner point.
    FT sm_sphere_radius = 5.0 * radius;
    FT sm_dichotomy_error = sm_distance * average_spacing/ 1000.0; // Dichotomy error must be << sm_distance
    Surface_3 surface(function,
                      Sphere(inner_point, sm_sphere_radius * sm_sphere_radius),
                      sm_dichotomy_error / sm_sphere_radius);

    // Defines surface mesh generation criteria
    CGAL::Surface_mesh_default_criteria_3<STr> criteria(sm_angle,  // Min triangle angle (degrees)
                                                        sm_radius * average_spacing,  // Max triangle size
                                                        sm_distance * average_spacing); // Approximation error

    std::cerr         << "  make_surface_mesh(sphere center=("<<inner_point << "),\n"
                      << "                    sphere radius="<<sm_sphere_radius<<",\n"
                      << "                    angle="<<sm_angle << " degrees,\n"
                      << "                    triangle size="<<sm_radius<<" * average spacing="<<sm_radius*average_spacing<<",\n"
                      << "                    distance="<<sm_distance<<" * average spacing="<<sm_distance*average_spacing<<",\n"
                      << "                    dichotomy error=distance/"<<sm_distance*average_spacing/sm_dichotomy_error<<",\n"
                      << "                    Manifold_with_boundary_tag)\n";

    // Generates surface mesh with manifold option
    STr tr; // 3D Delaunay triangulation for surface mesh generation
    C2t3 c2t3(tr); // 2D complex in 3D Delaunay triangulation
    CGAL::make_surface_mesh(c2t3,                                 // reconstructed mesh
                            surface,                              // implicit surface
                            criteria,                             // meshing criteria
                            CGAL::Manifold_with_boundary_tag());  // require manifold mesh

    // Prints status
    std::cerr << "Surface meshing: " << tr.number_of_vertices() << " output vertices" << std::endl;

    if (tr.number_of_vertices() == 0)
      return EXIT_FAILURE;

    // Converts to polyhedron
    CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, out);

    return EXIT_SUCCESS;
}


#endif //SWIG_CGAL_POISSON_SURFACE_RECONSTRUCTION_IMPL_H
