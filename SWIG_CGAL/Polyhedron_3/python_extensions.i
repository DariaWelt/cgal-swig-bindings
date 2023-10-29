// ------------------------------------------------------------------------------
// Copyright (c) 2011 GeometryFactory (FRANCE)
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// ------------------------------------------------------------------------------ 

%import "SWIG_CGAL/common.i"

//modifiers examples
%{
  #include "SWIG_CGAL/Polyhedron_3/modifier_example.h"
%}
%include "SWIG_CGAL/Polyhedron_3/Modifier_base.h"
%include "SWIG_CGAL/Common/triple.h"

%template(Polyhedron_3_Modifier_base)                  Modifier_base< Polyhedron_3_ >;

//simple modifiers
%template(Polyhedron_3_Modifier_1)                     Modifier_wrapper< Polyhedron_3_,Build_triangle<Polyhedron_3_::HalfedgeDS> >;
%template(Polyhedron_3_Modifier_2)                     Modifier_wrapper< Polyhedron_3_,Build_square<Polyhedron_3_::HalfedgeDS> >;

//Advanced modifier
%include "SWIG_CGAL//Common/triple.h"
%template(Integer_triple)  SWIG_CGAL::Triple<int,int,int>;
%define Triple_integer_range std::pair<Input_iterator_wrapper<Integer_triple,Integer_triple>,Input_iterator_wrapper<Integer_triple,Integer_triple> > %enddef

%{
  #include <SWIG_CGAL/Common/triple.h>
  typedef SWIG_CGAL::Triple<int,int,int> iInteger_triple;

  #include <SWIG_CGAL/Polyhedron_3/Modifier_base.h>
  #include <SWIG_CGAL/Kernel/Point_3.h>
  #include <CGAL/Polyhedron_3.h>
%}
SWIG_CGAL_declare_identifier_of_template_class(Integer_triple,SWIG_CGAL::Triple<int,int,int>)

//typemap for input iterators
SWIG_CGAL_input_iterator_typemap_in(Point_range,Point_3,Point_3,Point_3::cpp_base,SWIGTYPE_p_Point_3,"(LCGAL/Kernel/Point_3;)J",set_modifier_data)
SWIG_CGAL_input_iterator_typemap_in(Triple_integer_range,Integer_triple,Integer_triple,iInteger_triple,SWIGTYPE_p_SWIG_CGAL__TripleT_int_int_int_t,"(LCGAL/Polyhedron_3/Integer_triple;)J",set_modifier_data)

%template(Polyhedron_3_Modifier_3)                     Modifier_wrapper< Polyhedron_3_,Build_triangular_facets_from_point_range<Polyhedron_3_::HalfedgeDS> >;

%extend Modifier_wrapper< Polyhedron_3_,Build_triangular_facets_from_point_range<Polyhedron_3_::HalfedgeDS> >{
  void set_modifier_data(Point_range pt_range,Triple_integer_range int_range){
    std::copy (SWIG_CGAL::get_begin(pt_range),SWIG_CGAL::get_end(pt_range),$self->get_modifier_cpp_base().point_writer()); //copy points into the modifier
    std::copy (SWIG_CGAL::get_begin(int_range),SWIG_CGAL::get_end(int_range),$self->get_modifier_cpp_base().integer_triple_writer()); //copy triple of integer into the modifier
  }
}

