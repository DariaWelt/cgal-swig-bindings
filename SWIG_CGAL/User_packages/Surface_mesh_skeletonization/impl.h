#include <CGAL/Polyhedron_3.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <SWIG_CGAL/Kernel/Point_3.h>
#include <fstream>
#include <vector>

typedef Point_3 Point_3_wrapper;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron> Skeletonization;
typedef Skeletonization::Skeleton Skeleton;

typedef Skeleton::vertex_descriptor Skeleton_vertex;
typedef Skeleton::edge_descriptor Skeleton_edge;

typedef std::vector<Point> Polyline;
typedef std::vector<Polyline> Polylines;
typedef std::vector<Point_3_wrapper> Polyline_wrapper;
typedef std::vector<Polyline_wrapper> Polylines_wrapper;

//only needed for the display of the skeleton as maximal polylines
struct Display_polylines{
  const Skeleton &skeleton;
  Polylines &out;

  Display_polylines(const Skeleton &skeleton, Polylines &out)
    : skeleton(skeleton), out(out)
  {
  }

  void start_new_polyline()
  {
      out.push_back(Polyline());
  }

  void add_node(Skeleton_vertex v)
  {
    out.back().push_back(skeleton[v].point);
  }

  void end_polyline()
  {
  }
};

int skeletonization(const Polyhedron &tmesh, Polylines_wrapper &output_skeleton,
    Polylines_wrapper &output_correspondence)
{
  if (!CGAL::is_triangle_mesh(tmesh))
  {
    std::cout << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  Skeleton skeleton;

  CGAL::extract_mean_curvature_flow_skeleton(tmesh, skeleton);

  std::cout << "Number of vertices of the skeleton: " << boost::num_vertices(skeleton) << "\n";
  std::cout << "Number of edges of the skeleton: " << boost::num_edges(skeleton) << "\n";

  // Output all the edges of the skeleton.
  Polylines output_skeleton_unwrapped;
  Display_polylines display(skeleton, output_skeleton_unwrapped);
  CGAL::split_graph_into_polylines(skeleton, display);

  // Wrap points because linker error if do this in visitor, idk
  for (const auto &polyline : output_skeleton_unwrapped)
  {
    output_skeleton.push_back(Polyline_wrapper());
    for (const auto &point : polyline)
        output_skeleton.back().push_back(Point_3_wrapper(point));
  }

  // Output skeleton points and the corresponding surface points
  for (Skeleton_vertex v : CGAL::make_range(vertices(skeleton)))
    for (vertex_descriptor vd : skeleton[v].vertices)
    {
      output_correspondence.push_back(Polyline_wrapper());
      output_correspondence.back().push_back(Point_3_wrapper(skeleton[v].point));
      output_correspondence.back().push_back(Point_3_wrapper(get(CGAL::vertex_point, tmesh, vd)));
    }

  return EXIT_SUCCESS;
}

