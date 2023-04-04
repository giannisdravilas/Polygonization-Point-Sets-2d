#include <CGAL/Polygon_2.h>
#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<double>  K;
typedef CGAL::Polygon_2<K>              Polygon_2;
typedef K::Point_2                      Point_2;
typedef K::Segment_2                    Segment_2;
typedef std::vector<Segment_2>          Segments;

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