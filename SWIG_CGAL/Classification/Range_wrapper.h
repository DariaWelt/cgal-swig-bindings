// ------------------------------------------------------------------------------
// Copyright (c) 2020 GeometryFactory (FRANCE)
// Distributed under the Boost Software License, Version 1.0. (See accompany-
// ing file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
// ------------------------------------------------------------------------------

#ifndef SWIG_CGAL_CLASSIFICATION_RANGE_WRAPPER_H
#define SWIG_CGAL_CLASSIFICATION_RANGE_WRAPPER_H

#include <SWIG_CGAL/Point_set_3/Point_set_3.h>


#ifndef SWIG
struct Range_wrapper
{
  typename CGAL_PS3::template Property_range<int>::const_iterator begin;
  typename CGAL_PS3::template Property_range<int>::const_iterator end;

  Range_wrapper (typename Point_set_3_wrapper<CGAL_PS3>::Int_iterator ground_truth)
    : begin (ground_truth.get_cur()), end (ground_truth.get_end())
  { }

  int operator[] (const std::size_t& idx) const { return *(begin + idx); }
  std::size_t size() const { return (end - begin); }
};
#endif


#endif // SWIG_CGAL_CLASSIFICATION_RANGE_WRAPPER_H
