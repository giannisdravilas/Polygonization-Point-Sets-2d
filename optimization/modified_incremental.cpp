#include <iostream>
#include <string>
#include <vector>
#include <random>

#include <CGAL/Polygon_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/convex_hull_2.h>

#include "polygonization.h"
#include "optimization.h"

typedef CGAL::Simple_cartesian<double>  K;
typedef CGAL::Polygon_2<K>              Polygon_2;
typedef K::Point_2                      Point_2;
typedef std::vector<Point_2>            Points;
typedef K::Segment_2                    Segment_2;
typedef std::vector<Segment_2>          Segments;
typedef Polygon_2::Vertex_iterator      VertexIterator;

using namespace std;

Polygon_2 modified_incremental_algorithm(string inputFile, string outputFile, string initializationType, int edgeSelection, Points pointsToIgnore){

    Points points; // Vector with all the points

    // Read From File
    ifstream in(inputFile.c_str());
    string line;
    int line_count = 0;
    double areaFromInput = 0;
    while (getline(in, line)){ // Add points to vector
        istringstream iss(line);
        line_count++;

        // get info from second line
        if(line_count == 2){
            char one;
            string two, three, four, areaFromInputString;
            if (!(iss >> one >> two >> three >> four >> areaFromInputString)){  // continue if line has wrong format
                continue;
            }

            // Find convex hull area from second line
            areaFromInputString.erase(remove(areaFromInputString.begin(), areaFromInputString.end(), '"'), areaFromInputString.end());
            areaFromInputString.erase(remove(areaFromInputString.begin(), areaFromInputString.end(), '}'), areaFromInputString.end());
            areaFromInput = stod(areaFromInputString);
            continue;
        }

        // read lines with points and add them to the vector
        int index, pointX, pointY;
        if (!(iss >> index >> pointX >> pointY)){ // Skip lines that do not include 3 numbers (index pointX pointY)
            continue;
        }
        points.push_back(Point_2(pointX,pointY));
    }

    // we initiliaze the clock after the reading of the points from the file,
    // so that only the polygon construction time is measured
    auto started = chrono::high_resolution_clock::now();

    // sort points according to initialization type
    if (initializationType == "1a"){
        sort(points.begin(), points.end(), comparePointsDescendingX);
    }else if (initializationType == "1b"){
        sort(points.begin(), points.end(), comparePointsAscendingX);
    } else if (initializationType == "2a"){
        sort(points.begin(), points.end(), comparePointsDescendingY);
    } else if (initializationType == "2b"){
        sort(points.begin(), points.end(), comparePointsAscendingY);
    }

    // check if the first 3 points are collinear
    // if so, swap the third one with the next one, until the first 3 points are no more collinear
    int swap_index = 3;
    while (collinear(points[0], points[1], points[2])){
        if (swap_index == points.size())
            exit(1);
        iter_swap(points.begin() + 2, points.begin() + swap_index);
        swap_index++;
    }

    // create the triangle from the first 3 points
    Polygon_2 polygon;
    polygon.push_back(points[0]);
    polygon.push_back(points[1]);
    polygon.push_back(points[2]);
    
    // if the polygon is clockwise oriented, then reverse it to counter-clockwise (so that its has the same
    // orientation as the convex hull that will be created later)
    if (polygon.is_clockwise_oriented()){
        polygon.reverse_orientation();
    }
    vector<Polygon_2> tempPolygonsForSubdivision;
    // for every point of the vector created from file data
    for (auto p = points.cbegin() + 3; p != points.cend(); ++p) {

        // create convex hull of the current polygonal chain (at the 1st iteration, it is a triangle)
        Points out;
        CGAL::convex_hull_2(polygon.begin(), polygon.end(), back_inserter(out));        
        Polygon_2 convex_hull(out.cbegin(), out.cend());
        
        // find red edges of the convex hull
        Segments red_edges;
        for (auto edge = convex_hull.edges_begin(); edge != convex_hull.edges_end() ; ++edge){
            if (isVisible(*edge, *p, convex_hull)){
                red_edges.push_back(*edge);
            }
        }
                
        // now find the visible polygonal chain edges, which are "behind" the red edges of the convex hull
        // we save these edges to a new vector, called visible_polygon_edges
        Segments visible_polygon_edges;

        // for every red edge of the convex hull get start and end points
        for (auto red_edge = red_edges.begin(); red_edge != red_edges.end() ; ++red_edge){

            Point_2 red_edge_start = (*red_edge).point(0);
            Point_2 red_edge_end = (*red_edge).point(1);

            // and iterate the polygonal chain to find its first edge that has the same start
            // as the red edge
            for (auto edge1 = polygon.edges_begin(); edge1 != polygon.edges_end() ; ++edge1){

                // when found
                if ((*edge1).point(0) == red_edge_start){
                    bool end_found = false;

                    // iterate until a polygonal chain edge is found with the same end as the red edge
                    // and check if each intermidiate edge of the polygonal chain is visible to the point
                    for (auto edge2 = edge1; edge2 != polygon.edges_end() ; ++edge2){
                        if (isVisible(*edge2, *p, polygon)){
                            visible_polygon_edges.push_back(*edge2);
                        }
                        if ((*edge2).point(1) == red_edge_end){
                            end_found = true;
                            break;
                        }
                    }

                    // if we have reached the end of the polygonal chain without finding an edge of it with
                    // the same end as the end of the convex hull red edge, we iterate the polygonal chain
                    // from the beginning
                    // this may happen if the polygonal chain starts/ends "behind" a red edge of the convex hull

                    if (end_found == false){
                        for (auto edge2 = polygon.edges_begin(); edge2 != edge1 ; ++edge2){
                            if (isVisible(*edge2, *p, polygon)){
                                visible_polygon_edges.push_back(*edge2);
                            }

                            // this time it is sure that we will find the end
                            if ((*edge2).point(1) == red_edge_end){
                                end_found = true;
                                break;
                            }
                        }
                    }
                    break;
                }
            }
        }
        
        // we will now determine which of the visible polygonal chain edges will be replaced
        Segment_2 edge_to_replace;
        
        if (edgeSelection == 1){

            int range = visible_polygon_edges.size()-1 - 0 + 1;
            int num = rand() % range + 0;
            edge_to_replace = visible_polygon_edges[num];
            while(((edge_to_replace.source() == points[0]) && (edge_to_replace.target() == points[1]))
                || ((edge_to_replace.target() == points[0]) && (edge_to_replace.source() == points[1]))){
                int num = rand() % range + 0;
                edge_to_replace = visible_polygon_edges[num];
            }

        // so that the area added between new point and edge to be replaced source and target is minimized
        }else if (edgeSelection == 2){

            // initialize min with the the first possible area and iterate to find min
            double min = area((visible_polygon_edges[0]).target(), (visible_polygon_edges[0]).source(), *p);
            edge_to_replace = visible_polygon_edges[0];
            for (auto visible_edge = visible_polygon_edges.begin(); visible_edge != visible_polygon_edges.end() ; ++visible_edge){
                double calculated_area = area((*visible_edge).target(), (*visible_edge).source(), *p);

                if((visible_edge->source() == points[0]) && (visible_edge->target() == points[1])){
                    continue;
                }
                if((visible_edge->target() == points[0]) && (visible_edge->source() == points[1])){
                    continue;
                }

                if (calculated_area < min){
                    min = calculated_area;
                    edge_to_replace = *visible_edge;
                }
            }

        // so that the area added between new point and edge to be replaced source and target is maximized
        } else if (edgeSelection == 3){

            // initialize max with -1 (no need to do it as in minimizing) and iterate to find max
            double max = -1;
            for (auto visible_edge = visible_polygon_edges.begin(); visible_edge != visible_polygon_edges.end() ; ++visible_edge){
                double calculated_area = area((*visible_edge).target(), (*visible_edge).source(), *p);

                if((visible_edge->source() == points[0]) && (visible_edge->target() == points[1])){
                    continue;
                }
                if((visible_edge->target() == points[0]) && (visible_edge->source() == points[1])){
                    continue;
                }

                if (calculated_area > max){
                    max = calculated_area;
                    edge_to_replace = *visible_edge;
                }
            }
        }
        
        // (Only in Subdivision Step), In the last step if we have selected an edge that cannot be change, get another one
        bool flag = false;
        if(p == (points.end()-1)){
            for (auto visible_edge = visible_polygon_edges.begin(); visible_edge != visible_polygon_edges.end() ; ++visible_edge){
                if(visible_edge->source() == *(p-1)){
                    edge_to_replace = *visible_edge;
                    flag = true;
                    break;
                }
            }
            //backtrack
            if(flag == false){
                cout << "ERROR, calling convex hull" << endl;
                return Polygon_2();
            }
        }
        
        // iterate polygonal chain and add the new point before the end of the edge that will be replaced
        for (VertexIterator vi = polygon.vertices_begin(); vi != polygon.vertices_end(); ++vi){
            if (*vi == edge_to_replace.target()){
                polygon.insert(vi, *p);
                break;
            }
        }
    }

    auto done = chrono::high_resolution_clock::now();

    // points[1] must be conencted with points[2]
    for (auto vertex = polygon.vertices_begin(); vertex != polygon.vertices_end() ; ++vertex){
        if((*vertex == points[0])&&(vertex != (polygon.vertices_end()-1))){
            if(*(vertex+1) != points[1]){
                cout << "ERROR, calling convex hull" << endl;
                return Polygon_2();
            }
        }else if((*vertex == points[0])&&(vertex == (polygon.vertices_end()-1))){
            if(*(polygon.vertices_begin()) != points[1]){
                cout << "ERROR, calling convex hull" << endl;
                return Polygon_2();
            }
        }
    }

    return polygon;
}