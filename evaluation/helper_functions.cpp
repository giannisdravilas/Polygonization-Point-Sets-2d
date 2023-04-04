#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <list>
#include <random>
#include <algorithm>
#include <cmath>
#include <stdio.h>

#include <CGAL/Polygon_with_holes_2.h>
#include <iostream>

#include <CGAL/Polygon_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/intersections.h>
#include <CGAL/random_selection.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Kd_tree_rectangle.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Line_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/Boolean_set_operations_2.h>

#include "polygonization.h"
#include "optimization.h"

typedef CGAL::Simple_cartesian<double>          K;
typedef K::Point_2                              Point_2;
typedef std::vector<Point_2>                    Points;
typedef CGAL::Polygon_2<K>                      Polygon_2;
typedef CGAL::Polygon_with_holes_2<K>           Polygon_with_holes_2;
typedef Polygon_2::Edge_const_iterator          EdgeIterator;
typedef Polygon_2::Vertex_iterator              VertexIterator;
typedef K::Segment_2                            Segment_2;
typedef std::vector<Segment_2>                  Segments;       
typedef CGAL::Search_traits_2<K>                Traits;
typedef CGAL::Kd_tree<Traits>                   Tree;
typedef CGAL::Dimension_tag<2>                  Tag_2;
typedef CGAL::Kd_tree_rectangle<double, Tag_2>  KdRectangle;
typedef CGAL::Fuzzy_iso_box<Traits>             Fuzzy_iso_box;
typedef CGAL::Line_2<K>                         Line_2;
typedef CGAL::Triangle_2<K>                     Triangle_2;

typedef CGAL::Boolean_tag<false>            Tag_false;

// function checking whether an edge of the polygon is visible from point p
bool isVisible(Segment_2 edge, Point_2 p, Polygon_2 polygon){

    // get start and end points of edge
    Point_2 start = (edge).point(0);
    Point_2 end = (edge).point(1);

    // create segments from the point to the start, end and middle of the polygon edge and check for visibility
    // without considering the rest of the polygon
    Segment_2 point_with_edge_start(p, start);

    bool result1 = do_intersect(edge, point_with_edge_start);

    Segment_2 point_with_edge_end(p, end);
    bool result2 = do_intersect(edge, point_with_edge_end);

    Point_2 middle((start.x()+end.x())/2, (start.y()+end.y())/2);
    Segment_2 point_with_edge_middle(p, middle);
    bool result3 = do_intersect(edge, point_with_edge_middle);

    Segments vector_point_to_start_end_middle;
    vector_point_to_start_end_middle.push_back(point_with_edge_start);
    vector_point_to_start_end_middle.push_back(point_with_edge_end);
    vector_point_to_start_end_middle.push_back(point_with_edge_middle);

    // check if point and edge overlap in segment, before checking anything else (e.g. same x coordinates)
    for (auto segment_to_check = vector_point_to_start_end_middle.cbegin(); segment_to_check != vector_point_to_start_end_middle.cend(); ++segment_to_check){
        
        const auto result = intersection(edge, *segment_to_check);

        if (result) {
            if (const Segment_2* intersection_segment = boost::get<Segment_2>(&*result)) {
                return false;
            }
        }
    }

    // if edge is visible from the point, now consider the rest of the polygon
    if (result1 && result2 && result3){
        bool red = true;

        // for every edge of the polygon other than the edge for which we are checking visibility
        for (auto other_edge = polygon.edges_begin(); other_edge != polygon.edges_end() ; ++other_edge){

            if (edge != *other_edge){

                // check if any segment from the point intersects with this edge 
                for (auto segment_to_check = vector_point_to_start_end_middle.cbegin(); segment_to_check != vector_point_to_start_end_middle.cend(); ++segment_to_check){

                    const auto result = intersection(*segment_to_check, *other_edge);

                    // if it intersects anywhere else than the start or the end of the edge, then the first mentioned
                    // edge is not visible
                    // we allow intersection in start and end of the edge, because this edge may be a neighbor of the
                    // first mentioned edge, i.e. having one common point
                    if (result) {
                        if (const Segment_2* intersection_segment = boost::get<Segment_2>(&*result)) {
                            red = false;
                            return false;
                        } else {
                            const Point_2* intersection_point = boost::get<Point_2 >(&*result);
                            if (*intersection_point != start && *intersection_point != end){
                                red = false;
                                return false;
                            }
                        }
                    }
                }
            }
        }
        if (red){
            return true;
        }
    }
    return false;
}

// Sort points by acd. x
bool sortPoints(Point_2 A, Point_2 B){
    if(A.x() == B.x()){
        return (A.y() < B.y());
    }
    return (A.x() < B.x());
}

// Compare functions for point sorting in incremental
// if 1st dimension is identical, sort according to the 2nd one, so that there are no cases
// with collinear points where the 3rd point is inbetween the 1st and the 2nd one
bool comparePointsAscendingX(Point_2 p1, Point_2 p2){
    if (p1.x() == p2.x()){
        return p1.y() < p2.y();
    }
    return (p1.x() < p2.x());
}

bool comparePointsDescendingX(Point_2 p1, Point_2 p2){
    if (p1.x() == p2.x()){
        return p1.y() > p2.y();
    }
    return (p1.x() > p2.x());
}

bool comparePointsAscendingY(Point_2 p1, Point_2 p2){
    if (p1.y() == p2.y()){
        return p1.x() < p2.x();
    }
    return (p1.y() < p2.y());
}

bool comparePointsDescendingY(Point_2 p1, Point_2 p2){
    if (p1.y() == p2.y()){
        return p1.x() > p2.x();
    }
    return (p1.y() > p2.y());
}