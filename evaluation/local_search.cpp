#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>

#include "optimization.h"

using namespace std;

typedef CGAL::Simple_cartesian<double>  K;
typedef CGAL::Polygon_2<K>              Polygon_2;
typedef K::Segment_2                    Segment_2;
typedef K::Point_2                      Point_2;
typedef std::vector<Point_2>            Vertices;

typedef Polygon_2::Vertex_iterator VertexIterator;
typedef std::vector<VertexIterator> Vertex_Iterators;

// function that identifies if an edge is contained in (or lies on) a path of vertices
bool edge_in_path(Vertices path, Segment_2 edge){
    for (auto path_vertex = path.cbegin(); path_vertex != path.cend(); ++path_vertex){
        if (((Segment_2)(edge)).point(0) == *path_vertex || ((Segment_2)(edge)).point(1) == *path_vertex){
            return true;
        }
    }
    return false;
}

// function that takes a polygon, a vertex iterator and an int as arguments and returns a vector of points
// that represents a path in the polygon starting from vertex and containing i vertices
Vertices get_path(Polygon_2 &polygon, VertexIterator vertex, int i){
    int count = 0;
    Vertices path;

    for (auto new_vertex = vertex; new_vertex != polygon.vertices_end(); ++new_vertex){
        count++;
        path.push_back(*new_vertex);
        if (count == i){
            break;
        }
    }
    
    // in case we reach the end of the polygon before finding i vertices
    if (count < i){
        for (auto new_vertex = polygon.vertices_begin(); new_vertex != polygon.vertices_end(); ++new_vertex){
            count++;
            path.push_back(*new_vertex);
            if (count == i){
                break;
            }
        }
    }

    return path;
}

// function that takes a polygon, a path in this polygon and an edge of the polygon as arguments
// and erases this path from the current polygon by connecting the other vertices of the edges before and after
// the first and the last vertices of the path respectively with one another
// Then iterate the path backwards and add each vertex of it before the target point of the edge argument
Polygon_2 move(Polygon_2 &polygon, Vertices path, Segment_2 edge){
    Polygon_2 new_polygon = polygon;

    // erase path vertices
    if (path.size() >= 1){
        for (auto path_vertex = path.cbegin(); path_vertex != path.cend(); ++path_vertex){
            for (auto polygon_vertex = new_polygon.vertices_begin(); polygon_vertex != new_polygon.vertices_end() ; ++polygon_vertex){
                if (*path_vertex == *polygon_vertex){
                    new_polygon.erase(polygon_vertex);
                    break;
                }
            }
        }
        
    }

    // add them in their new position
    std::reverse(path.begin(), path.end());
    for (auto new_vertex = path.cbegin(); new_vertex != path.cend(); ++new_vertex){
        for (auto vertex = new_polygon.vertices_begin(); vertex != new_polygon.vertices_end(); ++vertex){
            if (*vertex == edge.target()){
                new_polygon.insert(vertex, *new_vertex);
                break;
            }
        }
    }

    // some times during the deletion of the path vertices a not simple polygon is produced,
    // which may turn to simple again while inserting the new vertices
    // in this case the orientation often changes, so we revert back to the original orientation
    if (polygon.orientation() != new_polygon.orientation()){
        new_polygon.reverse_orientation();
    }

    return new_polygon;
}

// function that takes a polygon, a path of this polygon and an edge of this polygon as arguments,
// calls move() as described above, decides if the movement made optimizes the result and returns
// the optimization result
// in order for a try to be successful, the new area must be greater than the old area (works for both
// maximization and minimazation problem due to signed area) and the new polygon must be simple
double try_move(Polygon_2 &polygon, Vertices path, Segment_2 edge){
    Polygon_2 new_polygon = move(polygon, path, edge);

    double new_area = new_polygon.area();
    double old_area = polygon.area();

    // if feasible solution, return the area improvement
    if (new_area > old_area && new_polygon.is_simple()){
        return new_area-old_area;
    }

    // else return 0
    return 0;
}

// local search core algorithm
//  for every edge to be replaced
//      for every path of max k length beginning with each one of the polygon vertices
//          check if the replacement produces a better result
//          keep track of the best result in a static variable along with the path and edge that result
//          to this
//  in the end apply the movement for the best result found
polygon_is_ok local_search(Polygon_2 init_polygon, int L, bool maximization_problem, double threshold, double max_time = 0){

    std::chrono::_V2::system_clock::time_point started, done;
    started = chrono::high_resolution_clock::now();

    // The polygons are oriented clockwise for area minimization and counterclockwise for area
    // maximization. This way, as signed areas are negative for clockwise orientation and positive for
    // counterclockwise orientation, solving the problem to maximize the area addresses both objective
    // functions. So, if we have a minimization problem we just reverse the orienation of the polygon.
    if (!maximization_problem){
        init_polygon.reverse_orientation();
    }

    Polygon_2 polygon = init_polygon;

    double max;
    do{
        
        done = chrono::high_resolution_clock::now();
        if ((chrono::duration_cast<chrono::milliseconds>(done-started).count() >= max_time )&& (max_time != 0)){
            //cout <<"local_search-timeout" << endl;
            if (!maximization_problem){
                polygon.reverse_orientation();
            }
            return polygon_is_ok{polygon, false};
        }

        // can be safely initialized to 0, as it will be compared to the return value of try_move(),
        // which will always be positive, for both maximization and minimization problems
        max = 0;
        Vertices best_path;
        Segment_2 best_edge;
        for (auto edge = polygon.edges_begin(); edge != polygon.edges_end(); ++edge){
            for (auto vertex = polygon.vertices_begin(); vertex != polygon.vertices_end(); ++vertex){
                done = chrono::high_resolution_clock::now();
                if ((chrono::duration_cast<chrono::milliseconds>(done-started).count() >= max_time )&& (max_time != 0)){
                    //cout <<"local_search-timeout" << endl;
                    if (!maximization_problem){
                        polygon.reverse_orientation();
                    }
                    return polygon_is_ok{polygon, false};
                }
                for (int i = 1; i <= L; i++){
                    Vertices path;
                    path = get_path(polygon, vertex, i);

                    // if current edge is contained in (or lies on) the path then continue for next k/L
                    // cause then the algorithm can not be applied
                    if (edge_in_path(path, *edge)){
                        continue;
                    }

                    // if current solution is the best one found so far
                    int current_solution;
                    if ((current_solution = try_move(polygon, path, *edge))){
                        if (current_solution > max){
                            best_path = path;
                            best_edge = *edge;
                            max = current_solution;
                        }
                    }
                }
            }
        }

        // if the best path has at least one vertex, then perform the replacement, else
        // there are no more better replacements to do, so break
        if (best_path.size()){
            polygon = move(polygon, best_path, best_edge);
        }

    }while(max >= threshold);

    // reverse back
    if (!maximization_problem){
        polygon.reverse_orientation();
    }

    return polygon_is_ok{polygon, true};
}