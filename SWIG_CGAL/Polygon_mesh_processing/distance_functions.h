#ifndef SWIG_CGAL_PMP_DISTANCE_FUNCTIONS_H
#define SWIG_CGAL_PMP_DISTANCE_FUNCTIONS_H

typedef CGAL::Parallel_if_available_tag Concurrency_tag;

#include <SWIG_CGAL/Polyhedron_3/Polyhedron_3.h>
#include <SWIG_CGAL/Polyhedron_3/typedefs.h>
#include <CGAL/Polygon_mesh_processing/distance.h>

namespace CGAL_SWIG {
  typedef Polyhedron_3_wrapper< Polyhedron_3_,SWIG_Polyhedron_3::CGAL_Vertex_handle<Polyhedron_3_>,SWIG_Polyhedron_3::CGAL_Halfedge_handle<Polyhedron_3_>,SWIG_Polyhedron_3::CGAL_Facet_handle<Polyhedron_3_> > Polyhedron_3_SWIG_wrapper;


  double approximate_Hausdorff_distance(Polyhedron_3_SWIG_wrapper& tm1, Polyhedron_3_SWIG_wrapper& tm2)
  {
    return CGAL::Polygon_mesh_processing::approximate_Hausdorff_distance<Concurrency_tag>(tm1.get_data(), tm2.get_data());
  }

   double approximate_Hausdorff_distance(Polyhedron_3_SWIG_wrapper& tm1, Polyhedron_3_SWIG_wrapper& tm2, int number_of_points_per_area_unit, bool use_monte_carlo_sampling)
  {
    return CGAL::Polygon_mesh_processing::approximate_Hausdorff_distance<Concurrency_tag>(tm1.get_data(), tm2.get_data(), CGAL::Polygon_mesh_processing::parameters::number_of_points_per_area_unit(number_of_points_per_area_unit));//, CGAL::Polygon_mesh_processing::parameters::use_monte_carlo_sampling(use_monte_carlo_sampling));
  }

  double approximate_symmetric_Hausdorff_distance(Polyhedron_3_SWIG_wrapper& tm1, Polyhedron_3_SWIG_wrapper& tm2)
  {
    return CGAL::Polygon_mesh_processing::approximate_symmetric_Hausdorff_distance<Concurrency_tag>(tm1.get_data(), tm2.get_data());
  }

  double approximate_symmetric_Hausdorff_distance(Polyhedron_3_SWIG_wrapper& tm1, Polyhedron_3_SWIG_wrapper& tm2, int number_of_points_per_area_unit, bool use_monte_carlo_sampling)
  {
    return CGAL::Polygon_mesh_processing::approximate_symmetric_Hausdorff_distance<Concurrency_tag>(
    tm1.get_data(), tm2.get_data(),
    CGAL::Polygon_mesh_processing::parameters::number_of_points_per_area_unit(number_of_points_per_area_unit)
    );//, CGAL::Polygon_mesh_processing::parameters::use_monte_carlo_sampling(use_monte_carlo_sampling));
  }
} // end namespace CGAL_SWIG

#endif //SWIG_CGAL_PMP_DISTANCE_FUNCTIONS_H