%module CGAL_Poisson_surface_reconstruction

%include "SWIG_CGAL/common.i"
Decl_void_type()

SWIG_CGAL_add_java_loadLibrary(CGAL_Poisson_surface_reconstruction)

%include "SWIG_CGAL/Common/Iterator.h"
%include "SWIG_CGAL/Kernel/typedefs.h"
%import  "SWIG_CGAL/Common/Macros.h"
%import  "SWIG_CGAL/Kernel/CGAL_Kernel.i"
%include  "SWIG_CGAL/Common/Output_iterator_wrapper.h"

//include files
%{
  #include <SWIG_CGAL/Kernel/typedefs.h>
  #include <SWIG_CGAL/Kernel/Point_3.h>
  #include <SWIG_CGAL/Polyhedron_3/all_includes.h>
  #include <SWIG_CGAL/Point_set_3/all_includes.h>
  #include <SWIG_CGAL/User_packages/Poisson_surface_reconstruction/impl.h>
  #include <SWIG_CGAL/Common/Output_iterator_wrapper.h>
%}

%types(Point_3*,Point_3); //needed so that the identifier SWIGTYPE_p_Point_3 is generated

%pragma(java) jniclassimports=%{
  import CGAL.Kernel.Point_3;
  import CGAL.Point_set_3.Point_set_3;
  import java.util.Iterator;
  import java.util.Collection;
  import CGAL.Polyhedron_3.Polyhedron_3;
%}

%pragma(java) moduleimports=%{
  import CGAL.Kernel.Point_3;
  import CGAL.Point_set_3.Point_set_3;
  import java.util.Iterator;
  import java.util.Collection;
  import CGAL.Polyhedron_3.Polyhedron_3;
%}

%define Integer_output_iterator boost::function_output_iterator<Container_writer<int,int> > %enddef
SWIG_CGAL_output_iterator_typemap_in(Integer_output_iterator,int,Integer,int,swig_types[0],"Ljava/lang/Integer;")

//import definitions of Polyhedron objects
%import "SWIG_CGAL/Polyhedron_3/CGAL_Polyhedron_3.i"
%import "SWIG_CGAL/Point_set_3/CGAL_Point_set_3.i"

//import Polyhedron_3 wrapper types
SWIG_CGAL_import_Polyhedron_3_SWIG_wrapper

void poisson_surface_reconstruction(const Point_set_3_wrapper<CGAL_PS3> &point_set, Polyhedron_3_SWIG_wrapper& output_polyhedron);

%{
  void poisson_surface_reconstruction(
    const Point_set_3_wrapper<CGAL_PS3> &point_set,
    Polyhedron_3_SWIG_wrapper& output_polyhedron)
  {
    reconstruction_poly(point_set.get_data(),
                        output_polyhedron.get_data());
  }
%}
