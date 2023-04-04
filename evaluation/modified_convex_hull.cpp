#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <random>

#include <CGAL/Polygon_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Triangle_2.h>

#include "polygonization.h"
#include "optimization.h"

typedef CGAL::Simple_cartesian<double>  K;
typedef CGAL::Polygon_2<K>              Polygon_2;
typedef K::Point_2                      Point_2;
typedef std::vector<Point_2>            Points;
typedef K::Segment_2                    Segment_2;
typedef Polygon_2::Edge_const_iterator  EdgeIterator;
typedef Polygon_2::Vertex_iterator      VertexIterator;
typedef std::vector<Segment_2>          Segments;       
typedef CGAL::Triangle_2<K>             Triangle_2; 

using namespace std;

Polygon_2 modified_convex_hull_algorithm(string inputFile, string outputFile, int edgeSelection, Points pointsToIgnore){
    
    srand((unsigned)time(NULL));    // Initialize time, for random selection of Edge (2nd algorithm)

    // -----------------------------INITIALIZE VECTORS-----------------------------
    Points pointVector; // Vector with all the points

    // Read From File
    ifstream in(inputFile.c_str());
    if(!in){
        cout << "Input-File: " << inputFile << " does not exist." <<endl;
        return Polygon_2();
    }
    string line;
    int line_count  = 0;
    double areaFromInput = 1;
    while(getline(in, line)){ // Add points to vector
        istringstream iss(line);

        // Second line. (Get area)
        line_count++;
        if(line_count == 2){
            char one;
            string two, three, four, areaFromInputString;
            if (!(iss >> one >> two >> three >> four >> areaFromInputString)){ // continue if line has wrong format
                continue;
            }

            // Find convex hull area from second line
            areaFromInputString.erase(remove(areaFromInputString.begin(), areaFromInputString.end(), '"'), areaFromInputString.end());
            areaFromInputString.erase(remove(areaFromInputString.begin(), areaFromInputString.end(), '}'), areaFromInputString.end());
            areaFromInput = stod(areaFromInputString);
            continue;
            
        }
        // Rest of the lines.
        int index, pointX, pointY;
        if (!(iss >> index >> pointX >> pointY)){ // Skip lines that do not include 3 numbers (index pointX pointY)
            continue;
        }
        pointVector.push_back(Point_2(pointX,pointY));
    }
 
    Points pointVectorConvex; // Vector with the points of the convex hull

    // Compute convex hull
    CGAL::convex_hull_2(pointVector.begin(), pointVector.end(), back_inserter(pointVectorConvex));

    // Remove the non-extreme Points(left and right) from the CH, if they are already in.(onyl in Subdivision Step)
    if((pointsToIgnore.size() == 2)){ // For first and last polygon
        Point_2 tempPoint1, tempPoint2;
        if(pointsToIgnore.at(0).y() < pointsToIgnore.at(1).y()){// First Polygon
            tempPoint1 = pointsToIgnore.at(0);
            tempPoint2 = pointsToIgnore.at(1);
            // Remove the non-extreme Point
            for (Points::iterator ptr = pointVectorConvex.begin(); ptr < pointVectorConvex.end(); ptr++){
                if(*ptr == tempPoint1){
                    pointVectorConvex.erase(ptr);
                    break;
                }
            }
            // Add the non-extreme Points after the extreme Point
            for (Points::iterator ptr = pointVectorConvex.begin(); ptr < pointVectorConvex.end(); ptr++){
                if(*ptr == tempPoint2){
                    pointVectorConvex.insert(ptr, tempPoint1);
                    break;
                }
            }
        }else{// Last Polygon
            tempPoint1 = pointsToIgnore.at(1);
            tempPoint2 = pointsToIgnore.at(0);
            // Remove the non-extreme Point
            for (Points::iterator ptr = pointVectorConvex.begin(); ptr < pointVectorConvex.end(); ptr++){
                if(*ptr == tempPoint1){
                    pointVectorConvex.erase(ptr);
                    break;
                }
            }
            // Add the non-extreme Points after the extreme Point
            for (Points::iterator ptr = pointVectorConvex.begin(); ptr < pointVectorConvex.end(); ptr++){
                if(ptr == pointVectorConvex.begin()){
                    if(*(pointVectorConvex.end()-1) == tempPoint2){
                        pointVectorConvex.insert(ptr, tempPoint1);
                        break;
                    }
                }else{
                    if(*(ptr-1) == tempPoint2){
                        pointVectorConvex.insert(ptr, tempPoint1);
                        break;
                    }
                } 
            }
        }
    }else if((pointsToIgnore.size() == 4)){ // For every other polygon
        // Remove the non-extreme Point
        for (Points::iterator ptr = pointVectorConvex.begin(); ptr < pointVectorConvex.end(); ptr++){
            if(*ptr == pointsToIgnore.at(1)){
                pointVectorConvex.erase(ptr);
                break;
            }
        }

        // Add the non-extreme Points after the extreme Point
        for (Points::iterator ptr = pointVectorConvex.begin(); ptr < pointVectorConvex.end(); ptr++){
            if(*ptr == pointsToIgnore.at(0)){
                if(ptr == pointVectorConvex.end()-1){
                    pointVectorConvex.insert(pointVectorConvex.begin(), pointsToIgnore.at(1));
                }else{
                    pointVectorConvex.insert(ptr+1, pointsToIgnore.at(1));
                }
                
                break;
            }
        }

        // Remove the non-extreme Point
        for (Points::iterator ptr = pointVectorConvex.begin(); ptr < pointVectorConvex.end(); ptr++){
            if(*ptr == pointsToIgnore.at(2)){
                //cout  << "(R)Delete: " << *ptr << endl;
                pointVectorConvex.erase(ptr);
                break;
            }
        }

        // Add the non-extreme Points after the extreme Point
        for (Points::iterator ptr = pointVectorConvex.begin(); ptr < pointVectorConvex.end(); ptr++){
            if(*ptr == pointsToIgnore.at(3)){
                pointVectorConvex.insert(ptr, pointsToIgnore.at(2));
                break;
            }
        }
    }

    // Make the Polygon (At the start, it has only the points of the Convex-Hull)
    Polygon_2 MyPolygon;
    for(Point_2 i: pointVectorConvex){
        MyPolygon.push_back(i);
    }
    
    //Sort all points(For set_difference to work)
    sort(pointVector.begin(), pointVector.end(), sortPoints); // Non-Convex
    sort(pointVectorConvex.begin(), pointVectorConvex.end(), sortPoints); // Convex

    // Remaining Points not yet in the Polygon
    Points diffVector; // Vector with the points of the convex hull
    set_difference(pointVector.begin(), pointVector.end(), pointVectorConvex.begin(), pointVectorConvex.end(), inserter(diffVector, diffVector.begin()));
    
    // -----------------------------ALGORITHM BASED ON CONVEX HULL-----------------------------
    // Until we have included all points in the Polygon:
    //      For every edge(E) in Polygon:
    //              Find nearest point(B)
    //      Select an edge (Select the edge for min or max area or random)
    //      Check if an edge is visible(from its nearest point B)
    //      Add Point(B) to the Polygon(remove edge(E) and insert two new ones)

    map<EdgeIterator, Point_2> mapOfNearestPointsLeft; // For edge selection for MIN and MAX erea
    bool isEdgeVisible = true; // Visibility of edge(Default: true)
    bool leftFlag = true;
    bool rightFlag = true;
    do{
        
        // Create a map of (EdgeIterator, NearestPoint)
        map<EdgeIterator, Point_2> mapOfNearestPoints;
        // For edge selection 2 and 3
        EdgeIterator selectedEdge; 

        // Edge And Point for selection
        Segment_2 selectedSegment;
        Point_2 selectedPoint;

        vector<EdgeIterator> tempEi;
        // Only if isVisible is true we will calculate the nearest points again
        if(isEdgeVisible == true){

            // Add to the map
            for (EdgeIterator ei = MyPolygon.edges_begin(); ei != MyPolygon.edges_end(); ++ei){
                
                // For every Edge not yet in the Polygon find its nearest Point
                int minDistance = CGAL::squared_distance(Segment_2(ei->source(), ei->target()), diffVector.at(0));
                Point_2 minPoint = diffVector.at(0); // Nearest Point
                for(Point_2 i: diffVector){
                    int tempMin = CGAL::squared_distance(Segment_2(ei->source(), ei->target()), i);
                    // If we have Points in the 'pointsToIgnore' vector
                    if(tempMin < minDistance){
                        minDistance = tempMin;
                        minPoint = i;
                    }
                }

                // Update map
                mapOfNearestPoints[ei] = minPoint;
            }


            // Update the map of nearest Points left
            mapOfNearestPointsLeft = mapOfNearestPoints;
        }

        // Pick Edge(random OR for max_area OR for min_area)
        if(edgeSelection == 1){ // Random
        
            auto item = mapOfNearestPointsLeft.begin();
            advance(item, rand() % mapOfNearestPointsLeft.size());
            selectedEdge = item->first;
            selectedSegment=Segment_2(selectedEdge->source(), selectedEdge->target());
            selectedPoint = item->second;
        }else if(edgeSelection == 2){ // MIN area
            double minArea = LLONG_MAX;

            // For every point and edge in "Nearest Points Map For Edge Selection" find the combination
            // that makes a Polygon with min area
            for (auto const& x : mapOfNearestPointsLeft){
                EdgeIterator selectedEdge2 = x.first;
                Segment_2 tempSegment = Segment_2(selectedEdge2->source(), selectedEdge2->target());
                Point_2 tempPoint = x.second;

                // Make triangle to compare the area
                Triangle_2 aTr(selectedEdge2->source(), selectedEdge2->target(), tempPoint);

                // Find min 
                double tempArea = MyPolygon.area() - abs(aTr.area());
                if(tempArea < minArea){
                    minArea = tempArea;
                    selectedSegment = tempSegment;
                    selectedPoint = tempPoint;
                    selectedEdge = selectedEdge2;
                }
            }
        }else if(edgeSelection == 3){ // MAX area
            double maxArea = -1;

            // For every point and edge in "Nearest Points Map For Edge Selection" find the combination
            // that makes a Polygon with max area
            for (auto const& x : mapOfNearestPointsLeft){
                EdgeIterator selectedEdge2 = x.first;
                Segment_2 tempSegment = Segment_2(selectedEdge2->source(), selectedEdge2->target());
                Point_2 tempPoint = x.second;

                // Make triangle to compare the area
                Triangle_2 aTr(selectedEdge2->source(), selectedEdge2->target(), tempPoint);

                // Find max 
                double tempArea = MyPolygon.area() - abs(aTr.area());
                if(tempArea > maxArea){
                    maxArea = tempArea;
                    selectedSegment = tempSegment;
                    selectedPoint = tempPoint;
                    selectedEdge = selectedEdge2;
                }
            }
        }
       
        // Check if Edge(selectedSegment) is visible through point(selectedPoint)
        isEdgeVisible = isVisible(selectedSegment, selectedPoint, MyPolygon);

        // Check if the slected segment can be changed(Only in Subdivision Step)
        if((count(pointsToIgnore.begin(), pointsToIgnore.end(), selectedSegment.target())>0) && 
            (count(pointsToIgnore.begin(), pointsToIgnore.end(), selectedSegment.source())>0)){
            isEdgeVisible = false;
        }
       
        // If not visible continue the loop.
        if(isEdgeVisible == false){
            // Erase the selected edge from the remaining pairs: <Edges, Point>
            mapOfNearestPointsLeft.erase(selectedEdge);
            // If there are no more visible edges, exit
            if(mapOfNearestPointsLeft.size() == 0){
                break;
            }
            continue;
        }
        
        // Add the Point to the Polygon
        for (VertexIterator vi = MyPolygon.vertices_begin(); vi != MyPolygon.vertices_end(); ++vi){
            if(*vi == selectedSegment.target()){
                MyPolygon.insert(vi, selectedPoint);
                break;
            }
        }
        
        // Remove Point from remaining points
        for(auto it = diffVector.begin(); it != diffVector.end(); ++it){
            if((it->x() == selectedPoint.x()) && (it->y() == selectedPoint.y())){
                diffVector.erase(it);
                break;
            }
        }
    }while(diffVector.size() > 0); // End loop when all points are in Polygon

    return MyPolygon;
}