// ------------------------------------------------------------------------------
// Copyright (c) 2013 GeometryFactory (FRANCE)
// Distributed under the Boost Software License, Version 1.0. (See accompany-
// ing file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
// ------------------------------------------------------------------------------

#ifndef SWIG_CGAL_POINT_SET_PROCESSING_3_H
#define SWIG_CGAL_POINT_SET_PROCESSING_3_H

#include <SWIG_CGAL/Kernel/Point_3.h>
#include <SWIG_CGAL/Kernel/Vector_3.h>
#include <SWIG_CGAL/Point_set_3/Point_set_3.h>

#include <CGAL/compute_average_spacing.h>
#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/jet_smooth_point_set.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/random_simplify_point_set.h>
#include <CGAL/remove_outliers.h>

#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

double compute_average_spacing(Point_set_3_wrapper<CGAL_PS3> point_set, int k)
{
  return CGAL::compute_average_spacing<Concurrency_tag> (point_set.get_data(), k);
}

void grid_simplify_point_set (Point_set_3_wrapper<CGAL_PS3> point_set, double epsilon)
{
  point_set.get_data().remove_from
    (CGAL::grid_simplify_point_set (point_set.get_data(), epsilon));
}

void jet_estimate_normals (Point_set_3_wrapper<CGAL_PS3> point_set, int k,
                           double neighbor_radius = 0.,
                           int degree_fitting = 2)
{
  point_set.get_data().add_normal_map();
  CGAL::jet_estimate_normals<Concurrency_tag>
    (point_set.get_data(), k,
     point_set.get_data().parameters().neighbor_radius (neighbor_radius).
     degree_fitting (degree_fitting));
}

void jet_smooth_point_set (Point_set_3_wrapper<CGAL_PS3> point_set, int k,
                           double neighbor_radius = 0.,
                           int degree_fitting = 2,
                           int degree_monge = 2)
{
  CGAL::jet_smooth_point_set<Concurrency_tag>
    (point_set.get_data(), k,
     point_set.get_data().parameters().neighbor_radius(neighbor_radius).
     degree_fitting(degree_fitting).
     degree_monge(degree_monge));
}

void mst_orient_normals (Point_set_3_wrapper<CGAL_PS3> point_set, int k,
                         double neighbor_radius = 0.,
                         typename Point_set_3_wrapper<CGAL_PS3>::Int_map
                         constrained_map = typename Point_set_3_wrapper<CGAL_PS3>::Int_map())
{
  if (constrained_map.is_valid())
    point_set.get_data().remove_from
      (CGAL::mst_orient_normals
       (point_set.get_data(), k,
        point_set.get_data().parameters().neighbor_radius (neighbor_radius).
        point_is_constrained_map(constrained_map.get_data())));
  else
    point_set.get_data().remove_from
      (CGAL::mst_orient_normals
       (point_set.get_data(), k,
        point_set.get_data().parameters().neighbor_radius (neighbor_radius)));
}

void pca_estimate_normals (Point_set_3_wrapper<CGAL_PS3> point_set, int k,
                           double neighbor_radius = 0.)
{
  point_set.get_data().add_normal_map();
  CGAL::pca_estimate_normals<Concurrency_tag>
    (point_set.get_data(), k, point_set.get_data().parameters().neighbor_radius(neighbor_radius));
}

void random_simplify_point_set (Point_set_3_wrapper<CGAL_PS3> point_set, double removed_percentage)
{
  point_set.get_data().remove_from
    (CGAL::random_simplify_point_set (point_set.get_data(), removed_percentage));
}

void remove_outliers (Point_set_3_wrapper<CGAL_PS3> point_set, int k,
                      double neighbor_radius = 0.,
                      double threshold_percent = 10.,
                      double threshold_distance = 0.)
{
  point_set.get_data().remove_from
    (CGAL::remove_outliers (point_set.get_data(), k,
                            point_set.get_data().parameters().neighbor_radius(neighbor_radius).
                            threshold_percent(threshold_percent).
                            threshold_distance(threshold_distance)));
}


#endif //SWIG_CGAL_POINT_SET_PROCESSING_3_H
