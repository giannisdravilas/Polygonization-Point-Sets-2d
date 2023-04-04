#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <random>
#include <algorithm>
#include <cmath>
#include <stdio.h>
#include <iostream>

#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/intersections.h>
#include <CGAL/random_selection.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Fuzzy_iso_box.h>
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
typedef CGAL::Search_traits_2<K>                Traits;
typedef CGAL::Kd_tree<Traits>                   Tree;
typedef CGAL::Fuzzy_iso_box<Traits>             Fuzzy_iso_box;
typedef CGAL::Line_2<K>                         Line_2;
typedef CGAL::Triangle_2<K>                     Triangle_2;
typedef CGAL::Boolean_tag<false>                Tag_false;

using namespace std;

// Regular simulated annealing
Polygon_2 simulated_annealing_algorithm(string inputFile, string outputFile, bool maximization_problem, Polygon_2 startPolygon, int L, string annealing, Points pointsToIgnore){
    srand((unsigned)time(NULL)); // Initialize time, for random selection of R (For MetropolisCriterion)
    double R = (double) rand()/(double)RAND_MAX;
    // -----------------------------ALGORITHM BASED ON SIMULATED ANNEALING-----------------------------   
    // Error
    if(startPolygon.size()==0){
        Polygon_2();
    }
    // Our polygon needs to be CW
    if (startPolygon.is_counterclockwise_oriented()){
        startPolygon.reverse_orientation();
    }
    
    // Calculate starting energy
    double N = startPolygon.size();
    double areaP = abs(startPolygon.area());
    double areaCHp = getAreaOfConvexHull(startPolygon);
    double myEnergy1 = ((maximization_problem)? ((double)N*(1.0-(double)areaP/((double)areaCHp))): ((double)N*(double)areaP/((double)areaCHp)));   

    // Start Algorithm
    double T = 1.0;
    int noOfChanges = 0;
    int maxIterations = 10000;
    while((T >=0) || (maxIterations == 0)){
        Points alreadyCheckedpoints; // This vector is filled after localStep. It contains already checked Points   
        vector<pair<Point_2, Point_2>> alreadyCheckedpointsGlobal; // This vector is filled after globalStep. It contains already checked pairs of Points.
        Points setOfPointsReturned; // This vector is filled after a localStep(4 points) or globalStep(5 points).

        // Only if validity is true(Or error) continue
        int validity = 0;
        do{
            setOfPointsReturned.clear();
            if(annealing.compare("local") == 0){// Do local step
                validity = localStep(startPolygon, &alreadyCheckedpoints, &setOfPointsReturned);
            }else if(annealing.compare("global") == 0){// Do global step
                validity = globalStep(startPolygon, &alreadyCheckedpointsGlobal, &setOfPointsReturned, pointsToIgnore);

                // If we have subdivision Step, dont change the pq, qr edges
                if((count(pointsToIgnore.begin(), pointsToIgnore.end(), setOfPointsReturned[3])>0) && 
                    (count(pointsToIgnore.begin(), pointsToIgnore.end(), setOfPointsReturned[4])>0)){
                    validity = 0;
                }
                if((count(pointsToIgnore.begin(), pointsToIgnore.end(), setOfPointsReturned[0])>0) && 
                    (count(pointsToIgnore.begin(), pointsToIgnore.end(), setOfPointsReturned[1])>0)){
                    validity = 0;
                }
                if((count(pointsToIgnore.begin(), pointsToIgnore.end(), setOfPointsReturned[1])>0) && 
                    (count(pointsToIgnore.begin(), pointsToIgnore.end(), setOfPointsReturned[2])>0)){
                    validity = 0;
                }
            }
            // If we reached max iterations
            maxIterations--;
            if(maxIterations == 0){
                validity = -1;
            }
        }while(validity == 0);

        // If we have exited the loop because we didn't have any points left to check, exit the algorithm-loop.
        if(validity == -1){
            break;
        }

        // get new area
        double newArea = newAreaForTwoTriangles(annealing, setOfPointsReturned, startPolygon);

        // Calculate the new energy and the metropolis criterion.
        double myEnergy2 = ((maximization_problem)? ((double)N*(1.0-(double)newArea/((double)areaCHp))): ((double)N*(double)newArea/((double)areaCHp)));
        double DE = (double)myEnergy2 - (double)myEnergy1;
        double metropolisCriterion = (double)exp(-(double)DE/((double)T));

        // Only if DE<0 and if Metropolis criterion holds, apply the transition
        if((DE < 0)){
            if(metropolisCriterion > R){
                updatePolygon(&startPolygon, setOfPointsReturned,annealing);
                noOfChanges++;
                myEnergy1 = myEnergy2;
            } 
        }
        
        T = (double)T - (double)1.0/((double)L);// Update  "T" 
    }  

    // ------------------------------------------Write in Output-File------------------------------------------
    // reverse back
    if (startPolygon.is_clockwise_oriented()){
        startPolygon.reverse_orientation();
    }
    return startPolygon;
}

// Find rectangle for kDTree based search
Fuzzy_iso_box findRectagular(Polygon_2 startPolygon, Point_2 selectedPoint, Point_2* myA, Point_2* myB, Point_2* myC, Point_2* myD){

    // Find rectangle Vertexes(p,q,r,s)
    Point_2 A,B,C,D;
    for (VertexIterator currentPoint = startPolygon.vertices_begin(); currentPoint != startPolygon.vertices_end(); ++currentPoint){
        if(*currentPoint == selectedPoint){
            if(*currentPoint == *(startPolygon.vertices_end()-1)){
                A = *(currentPoint-1);
                B = *(currentPoint);
                C = *(startPolygon.vertices_begin());
                D = *(startPolygon.vertices_begin()+1);
            }else if(*(currentPoint+1) == *(startPolygon.vertices_end()-1)){
                A = *(currentPoint-1);
                B = *(currentPoint);
                C = *(currentPoint+1);
                D = *(startPolygon.vertices_begin());
            }else if(*currentPoint == *startPolygon.vertices_begin()){
                A = *(startPolygon.vertices_end()-1);
                B = *(currentPoint);
                C = *(currentPoint+1);
                D = *(currentPoint+2);  
            }else{
                A = *(currentPoint-1);
                B = *(currentPoint);
                C = *(currentPoint+1);
                D = *(currentPoint+2);
            }
            break;
        }
    }

    // Update returing Points(p,q,r,s)
    *myA = A;
    *myB = B;
    *myC = C;
    *myD = D;

    // Find max/min x and y of the four points. Used for calculating the rectangle
    Points temp = {A, B, C, D};
    Point_2 maxPointX(A), maxPointY(A), minPointX(A), minPointY(A);
    int maxX = A.x() , maxY = A.y(), minX = A.x(), minY = A.y();
    for(int i=0; i<4; i++){
        if(temp[i].x()>maxX){
            maxX=temp[i].x();
            maxPointX=temp[i];
        }
        if(temp[i].y()>maxY){
            maxY=temp[i].y();
            maxPointY=temp[i];
        }
        if(temp[i].x()<minX){
            minX=temp[i].x();
            minPointX=temp[i];
        }
        if(temp[i].y()<minY){
            minY=temp[i].y();
            minPointY=temp[i];
        }
        
    }   

    // Make the 4 lines that are parallel to the x and y axis 
    const Line_2 up(maxPointY, Point_2(maxPointY.x()+1, maxPointY.y())), right(maxPointX,Point_2(maxPointX.x(), maxPointX.y()+1));
    const Line_2 left(minPointX, Point_2(minPointX.x(), minPointX.y()+1)), down(minPointY, Point_2(minPointY.x()+1, minPointY.y()));

    // Find the two interestion Points of the previous 4 lines
    Point_2 lastUpRightPoint, lastDownLeftPoint;
    auto result1 = intersection(up, right);
    if (result1) {
        const Point_2* lUPoint = boost::get<Point_2>(&*result1);
        lastUpRightPoint = (*lUPoint);
    }
    const auto result2 = intersection(down, left);
    if (result2) {
        const Point_2* lDPoint = boost::get<Point_2>(&*result2);
        lastDownLeftPoint = (*lDPoint);
    }

    // Calculate and return the box
    return Fuzzy_iso_box(lastDownLeftPoint, lastUpRightPoint);
}

// Function for Local Step
int localStep(Polygon_2 startPolygon, Points* alreadyCheckedpoints, Points* setOf4Points){

    // Select Random Point
    Point_2 selectedPoint;
    random_selection(startPolygon.vertices_begin(), startPolygon.vertices_end(), 1, &selectedPoint);

    // If we have already checked this point, select another one.
    if(count(alreadyCheckedpoints->begin(), alreadyCheckedpoints->end(), selectedPoint)>0){
        return 0;
    }

    // Make kd-Tree
    Tree kdTree;
    kdTree.insert(startPolygon.vertices_begin(), startPolygon.vertices_end());
    
    // Find iso box for tree search and the for points(p,q,r,s)
    Point_2 A,B,C,D;
    Fuzzy_iso_box myS = findRectagular(startPolygon, selectedPoint, &A,&B,&C,&D);

    // Update returning Vector
    setOf4Points->push_back(Point_2(A));
    setOf4Points->push_back(Point_2(B));
    setOf4Points->push_back(Point_2(C));
    setOf4Points->push_back(Point_2(D));

    // Make kdTree search.
    list<Point_2>kdSearchRectPoints;
    kdTree.search(back_inserter(kdSearchRectPoints), myS);

    // Remove the points.
    kdSearchRectPoints.remove(A);
    kdSearchRectPoints.remove(B);
    kdSearchRectPoints.remove(C);
    kdSearchRectPoints.remove(D);

    // For every point located in the range, test whether one of the edges incident to the point intersects one of pr or qs
    bool validity = check_validity_for_local_step(Segment_2(A, C), Segment_2(B, D), A, D, startPolygon, kdSearchRectPoints);

    // Only if validity is false, update alreadyCheckedpoints nad check its size.
    if(validity == false){
        // Update returning vector
        alreadyCheckedpoints->push_back(Point_2(selectedPoint));
        if(alreadyCheckedpoints->size() == startPolygon.size()){
            return -1;
        }
    }

    return validity;
}

// Function for Global Step
int globalStep(Polygon_2 startPolygon, vector<pair<Point_2, Point_2>>* alreadyCheckedpointsGlobal, Points* setOf5Points, Points pointsToIgnore){

    // Select 2 Random Points
    Points temp;
    random_selection(startPolygon.vertices_begin(), startPolygon.vertices_end(), 2, back_inserter(temp));
    pair<Point_2, Point_2> points(temp[0], temp[1]);
    
    // If we have already checked one of the points, select 2 another ones.
    if(count(alreadyCheckedpointsGlobal->begin(), alreadyCheckedpointsGlobal->end(), points)>0){
        return false;
    }

    // Find q,r,p,t
    // q = point.first, r = q-1, p = q+1, s = point.second, t = s +1
    Point_2 R, P, Q = points.first;
    Point_2 T, S = points.second;
    for (VertexIterator currentPoint = startPolygon.vertices_begin(); currentPoint != startPolygon.vertices_end(); ++currentPoint){
        if(*currentPoint == Q){
            if(*currentPoint == *(startPolygon.vertices_end()-1)){
                P = *(currentPoint-1);
                R = *(startPolygon.vertices_begin());
            }else if(*currentPoint == *startPolygon.vertices_begin()){
                P = *(startPolygon.vertices_end()-1);
                R = *(currentPoint+1);    
            }else{
                P = *(currentPoint-1);
                R = *(currentPoint+1);
            }
        }
        if(*currentPoint == S){
            if(*currentPoint == *(startPolygon.vertices_end()-1)){
                T = *(startPolygon.vertices_begin());               
            }else{
                T = *(currentPoint+1);
            }
        }
    }
    
    // Update returning Vector
    setOf5Points->push_back(Point_2(R));
    setOf5Points->push_back(Point_2(Q));
    setOf5Points->push_back(Point_2(P));
    setOf5Points->push_back(Point_2(S));
    setOf5Points->push_back(Point_2(T));

    //Check vilidity
    Segment_2 pr(P, R), sq(S, Q), qt(Q, T);
    bool validity = check_validity_for_global_step(pr, sq, qt, startPolygon, P, R, S, Q, T);

    // Added edge cases if we checked all the points
    // Only if validity is false, update alreadyCheckedpointsGlobal nad check its size.
    if(validity == false){
        // Update returning vector
        alreadyCheckedpointsGlobal->push_back(pair<Point_2, Point_2>(points));
        if(alreadyCheckedpointsGlobal->size() == (startPolygon.size()-1)*startPolygon.size()){
            return -1;
        }
    }
    
    return validity;// Return validity
}

// Function for Subdivision Step
Polygon_2 simulated_annealing_algorithm_subdivision(string inputFile, string outputFile, bool maximization_problem, int L, string annealing, int edgeSelection, string startAlgorithm){

    vector<Polygon_2> vectorOfPolygons; 
    vector<pair<Segment_2, Segment_2>> vectorOftLowerHullEdge;
    // Make Small Polygons
    make_small_polygons_and_mark_segments(inputFile, &vectorOfPolygons, &vectorOftLowerHullEdge);
    std::chrono::_V2::system_clock::time_point started, done;

    // Solve the polygonization problem for each subset using global step.
    // And Force the polygon from convex_hull/incremental to have the initial right and left segments.
    list<Polygon_2> endPolygonsOfGLboalStep; // New list of Polgons, Required because we need new initial polygons.
    vector<Points> vectorpointsToIgnore;
    started = chrono::high_resolution_clock::now();// Start time

    double initial_polygon_area; // Sum area from all the polygons returned from convex_hull/incremental
    for(int i=0;i<vectorOfPolygons.size();i++){
        
        // Make Temp input file for convexHull/Incremental
        string myInputFile = make_input_file(i, vectorOfPolygons);

        cout  << "--------------------------------------  Start Polygon --------------------------------------" << i << endl;

        // Find the points from the marked segments
        Points pointsToIgnore = make_points_to_ignore(i, vectorOfPolygons.size(), vectorOftLowerHullEdge); 
        vectorpointsToIgnore.push_back(pointsToIgnore);

        // Call modified convex_hull_algorithm or incremental
        // The marked segments need to appear in the partial results, and they must not be modified
        Polygon_2 temptemppoly;
        int tempEdge = ((maximization_problem)? (3):(2));
        if(startAlgorithm == "convex_hull"){
            temptemppoly = modified_convex_hull_algorithm(myInputFile, outputFile, tempEdge, pointsToIgnore);
        }else if(startAlgorithm == "incremental"){
            temptemppoly = modified_incremental_algorithm(myInputFile, outputFile, "1b", tempEdge, pointsToIgnore);
            if(temptemppoly.size()==0){// If ERROR call convex-hull
                temptemppoly = modified_convex_hull_algorithm(myInputFile, outputFile, tempEdge, pointsToIgnore);
            }
        }
        initial_polygon_area = initial_polygon_area + abs(temptemppoly.area());// update area

        // Delete temp Input files for convexHull/Incremental
        if(remove(myInputFile.c_str()) != 0){
            perror( "Error deleting file" );
        }

        // Call global Step
        Polygon_2 tempPoly = simulated_annealing_algorithm(inputFile, to_string(i) +outputFile, maximization_problem, temptemppoly, L, "global", pointsToIgnore);
        if (tempPoly.is_clockwise_oriented()){
            tempPoly.reverse_orientation();
        }
        cout << fixed << "Difference in Area: " << abs(tempPoly.area()) - abs(temptemppoly.area())<< endl;
        endPolygonsOfGLboalStep.push_back(tempPoly);
        cout  << "--------------------------------------  End Polygon   --------------------------------------" << i << endl;
    }

    // Join Polygons
    list<Polygon_with_holes_2> symmR;
    CGAL::join(endPolygonsOfGLboalStep.begin(), endPolygonsOfGLboalStep.end(), back_inserter(symmR), Tag_false());
    Polygon_2 LastPolygon = symmR.front().outer_boundary();
    if (LastPolygon.is_clockwise_oriented()){
        LastPolygon.reverse_orientation();
    }
    
    // Remove pq and qr 
    for (VertexIterator iter = LastPolygon.begin(); iter!=LastPolygon.end(); ++iter){
        // For every vector of points to ignore
        for(int i=1; i < vectorpointsToIgnore.size();i++){
            // For every polygon erase one point
            Point_2 first = vectorpointsToIgnore.at(i).at(1);
            if((*iter == first)&&(*iter != *LastPolygon.begin())){
                LastPolygon.erase(iter-1);
            }else if((*iter == first)&&(*iter == *LastPolygon.begin())){
                LastPolygon.erase(LastPolygon.end()-1);
            }
        }
    }
    // Call Local Step
    Polygon_2 endPolygon = simulated_annealing_algorithm(inputFile, outputFile +".txt", maximization_problem, LastPolygon, L, "local", Points());
    // reverse back
    if (endPolygon.is_clockwise_oriented()){
        endPolygon.reverse_orientation();
    }
    done = chrono::high_resolution_clock::now(); // end timer
    // Write in output file
    outputForSA(inputFile, outputFile, endPolygon, maximization_problem, initial_polygon_area, started, done, annealing);
    return endPolygon;
}

// Check validity in startPolygon
bool check_validity_for_local_step(Segment_2 pr, Segment_2  qs, Point_2 A, Point_2 D, Polygon_2 startPolygon, list<Point_2> kdSearchRectPoints){

    // Check if pr and qs intersect each other.
    if (do_intersect(pr, qs)) { 
        return false; 
    }

    // Find the point before A and after D
    Point_2 pointBeforeA, pointAfterD;
    for (VertexIterator currentPoint = startPolygon.vertices_begin(); currentPoint != startPolygon.vertices_end(); ++currentPoint){
        if(*currentPoint == A){
            if(*currentPoint == *startPolygon.vertices_begin()){
                pointBeforeA = *(startPolygon.vertices_end()-1);
            }else{
                pointBeforeA = *(currentPoint-1);
            }
        }
        if(*currentPoint == D){
            if(*currentPoint == *(startPolygon.vertices_end()-1)){
                pointAfterD = *startPolygon.vertices_begin();
            }else{
                pointAfterD= *(currentPoint+1);
            }
        }
    } 

    // Check if qs intersect the segment before A and if pr intersect the segment after D.
    Segment_2 BeforeA(pointBeforeA, A), AfterD(D, pointAfterD);
    if (do_intersect(pr, AfterD) ||do_intersect(qs, BeforeA)) { // Check if they intersect each other.
        return false; 
    }

    // For every point x in rectangle, check if pr and qs intersect its segments
    for(auto it = kdSearchRectPoints.begin(); it != kdSearchRectPoints.end(); ++it){ // For each located Point
        Point_2 locatedPoint = *it;
        // Iterate the polygon until we find the point that we want.
        for (VertexIterator currentPoint = startPolygon.vertices_begin(); currentPoint != startPolygon.vertices_end(); ++currentPoint){
            if(*currentPoint == locatedPoint){
                if(*currentPoint == *(startPolygon.vertices_end()-1)){
                    Segment_2 xA(*currentPoint, *(currentPoint-1)), xY(*currentPoint, *startPolygon.vertices_begin());
                    if (intersection(pr, xA)||intersection(pr, xY)||intersection(qs, xA)||intersection(pr, xY)) {
                        return false;
                    }
                }else if(*currentPoint == *startPolygon.vertices_begin()){
                    Segment_2 xA(*currentPoint, *(startPolygon.vertices_end()-1)), xY(*currentPoint, *(currentPoint+1));
                    if (intersection(pr, xA)||intersection(pr, xY)||intersection(qs, xA)||intersection(pr, xY)) {
                        return false;
                    }
                }else{
                    Segment_2 xA(*currentPoint, *(currentPoint-1)), xY(*currentPoint, *(currentPoint+1));
                    if (intersection(pr, xA)||intersection(pr, xY)||intersection(qs, xA)||intersection(pr, xY)) {
                        return false;
                    }
                }
            }
        } 
    }
    return true;
}

// Write in output file for simulated annealing
void outputForSA(string inputFile, string outputFile, Polygon_2 startPolygon, bool maximization_problem, double initial_polygon_area, std::chrono::_V2::system_clock::time_point started, std::chrono::_V2::system_clock::time_point done, string annealing){
    
    cout << "output file is " << outputFile << endl;
    // Write in output file
    ofstream outputFileClass;
    outputFileClass.open(outputFile.c_str());
    outputFileClass << "Optimal Area Polygonization" << endl;

    // Write Verticles
    for (VertexIterator iter = startPolygon.begin(); iter!=startPolygon.end(); ++iter){
        outputFileClass << fixed << *iter << endl;
    }

    // Write Edges
    for (EdgeIterator ei = startPolygon.edges_begin(); ei != startPolygon.edges_end(); ++ei){
        outputFileClass << fixed << ei->source() <<" "<< ei->target() << endl;
    } 

    // Write rest information
    outputFileClass << "Algorithm: simulated_annealing_subdivision [" << ((maximization_problem)? "max]": "min]") << endl;
    outputFileClass << fixed << "Area_initial: " << initial_polygon_area << endl;
    double polygon_area = startPolygon.area();
    outputFileClass << fixed << "Area: " << polygon_area << endl;

    // Print for testing
    cout  << "Is Polygon Simple: " << startPolygon.is_simple() << endl;

    // Get starting polygon area(before convex_hull or incremtnal)
    double areaFromInput = 0;
    ifstream infile(inputFile);
    if (infile.good()){
        string sLine;
        getline(infile, sLine);
        getline(infile, sLine);
        areaFromInput = stod(std::string(&sLine[38], &sLine[sLine.size()-2]));
        cout << fixed << "Area From Input:  "<< setfill('0') << setw(17) << areaFromInput << endl;
    }
    cout << fixed << "Area Initial:     " << setfill('0') << setw(17) << initial_polygon_area << endl;
    cout << fixed << "Area End:         "<< setfill('0') << setw(17) << polygon_area << endl;
    cout << fixed << "Difference in Area: " << polygon_area - initial_polygon_area<< endl;
    outputFileClass << fixed << "Ratio_initial: " << initial_polygon_area/areaFromInput << endl;
    outputFileClass << fixed << "Ratio: " << polygon_area/initial_polygon_area << endl;
    outputFileClass << fixed << "Construction Time: " << chrono::duration_cast<chrono::milliseconds>(done-started).count() << " ms" << endl;
    outputFileClass.close();

    cout << fixed << "Construction Time: " << chrono::duration_cast<chrono::milliseconds>(done-started).count() << " ms" << endl;
    
    return;
}

// Find new are for simulated annealing
double newAreaForTwoTriangles(string annealing, Points setOfPointsReturned,Polygon_2 startPolygon){
    // Calculate the area of the new Polygon
    Triangle_2 aTr, bTr;
    if(annealing.compare("local") == 0){
        // Get the new triangles
        aTr = Triangle_2(setOfPointsReturned[0], setOfPointsReturned[2], setOfPointsReturned[1]);
        bTr = Triangle_2(setOfPointsReturned[2], setOfPointsReturned[3], setOfPointsReturned[1]);
    }else if(annealing.compare("global") == 0){
        // Get the new triangles
        aTr = Triangle_2(setOfPointsReturned[2], setOfPointsReturned[0], setOfPointsReturned[1]);
        bTr = Triangle_2(setOfPointsReturned[3], setOfPointsReturned[4], setOfPointsReturned[1]);
    }
    return (double)abs(startPolygon.area() + aTr.area() - bTr.area());
}

// Update polygon for simulated annealing
void updatePolygon(Polygon_2* startPolygon, Points setOfPointsReturned, string annealing){
    if(annealing.compare("local") == 0){
        // Remove q
        for (VertexIterator currentPoint = startPolygon->vertices_begin(); currentPoint != startPolygon->vertices_end(); ++currentPoint){
            if(*currentPoint == setOfPointsReturned[1]){
                startPolygon->erase(currentPoint);
                break;
            }
        }
        // Insert q before s.
        for (VertexIterator currentPoint = startPolygon->vertices_begin(); currentPoint != startPolygon->vertices_end(); ++currentPoint){
            if(*currentPoint == setOfPointsReturned[3]){
                startPolygon->insert(currentPoint, setOfPointsReturned[1]);
                break;
            }
        }
    }else if(annealing.compare("global") == 0){
        // Removing q and connect p and r
        for (VertexIterator currentPoint = startPolygon->vertices_begin(); currentPoint != startPolygon->vertices_end(); ++currentPoint){
            if(*currentPoint == setOfPointsReturned[1]){
                startPolygon->erase(currentPoint);
                break;
            }
        }
        // Insert q before t.
        for (VertexIterator currentPoint = startPolygon->vertices_begin(); currentPoint != startPolygon->vertices_end(); ++currentPoint){
            if(*currentPoint == setOfPointsReturned[4]){
                startPolygon->insert(currentPoint, setOfPointsReturned[1]);
                break;
            }
        }
    }
    return;
}

// Get are of convex hull from given Polygon
double getAreaOfConvexHull(Polygon_2 startPolygon){
    // Find Convex Hull of our input polygon.(Used for calculating energy)
    Points pointVectorConvex; // Vector with the points of the convex hull
    CGAL::convex_hull_2(startPolygon.vertices_begin(), startPolygon.vertices_end(),back_inserter(pointVectorConvex));
    Polygon_2 MyConvexHullPolygon;
    for(Point_2 i: pointVectorConvex){
        MyConvexHullPolygon.push_back(i);
    }
    return (double)abs(MyConvexHullPolygon.area());
}

// Visibility for global step
bool check_validity_for_global_step(Segment_2 pr, Segment_2 sq, Segment_2 qt, Polygon_2 startPolygon, Point_2 P, Point_2 R, Point_2 S, Point_2 Q, Point_2 T){

    /// If the new segemtns intersect eachother
    if (do_intersect(pr, sq) || do_intersect(pr, qt)) { 
        return false;
    }

    // Check if pr, sq or qt intersects any other segment
    for(const Segment_2& e  : startPolygon.edges()){
        // Pr 
        const auto result = intersection(pr, e);
        if (result) {
            if (const Point_2* intersection_point = boost::get<Point_2 >(&*result)) {
                if (*intersection_point != P && *intersection_point != R){
                    return false;
                }
            }else{
                return false;
            }
        }

        // Sq
        const auto result1 = intersection(sq, e);
        if (result1) {
            if (const Point_2* intersection_point = boost::get<Point_2 >(&*result1)) {
                if (*intersection_point != S && *intersection_point != Q){
                    return false; 
                }
            }else{
                return false;
            }
        }

        // Qt
        const auto result2 = intersection(qt, e);
        if (result2) {
            if (const Point_2* intersection_point = boost::get<Point_2 >(&*result2)) {
                if (*intersection_point != Q && *intersection_point != T){
                    return false; 
                }
            }else{
                return false;
            }
        }
    }

    return true;
}

// Make sub-Polygons for simulated annealing-subdivision
void make_small_polygons_and_mark_segments(string inputFile, vector<Polygon_2>* vectorOfPolygons, vector<pair<Segment_2, Segment_2>>* vectorOftLowerHullEdge){

    // Vector with all the points
    Points tempPoints; 
    ifstream in(inputFile.c_str());
    if(!in){
        cout << "Input-File: " << inputFile << " does not exist." <<endl;
        return;
    }

    // Read from inputFile
    string line;
    int line_count = 0;
    while(getline(in, line)){ // Add points to vector
        istringstream iss(line);

        // Second line. (Get area)
        line_count++;
        if(line_count == 2){
            continue;
        }

        // Rest of the lines.
        int index, pointX, pointY;
        if (!(iss >> index >> pointX >> pointY)){ // Skip lines that do not include 3 numbers (index pointX pointY)
            continue;
        }
        tempPoints.push_back(Point_2(pointX,pointY));
    }
    
    int m = tempPoints.size()/10; // Number of points in each subset
    double k = ceil((tempPoints.size()-1)/(m-1)); // Number of subsets

    // Sort points based on the lexicographic order (increasing x).
    sort(tempPoints.begin(), tempPoints.end(), sortPoints); 

    // Make 'k' small polygons, Every two subsets located in adjacent "slabs" share a point.
    // Also check if the 2 conditions on the LH edges are met.
    int countOfPoints = 0; // Number of points ready each time to be inserted in a polygon
    int maxPointsPerSet = m; // For polygon 0 its is = m, and for every other polygon its m-1

    // Iterators for making the polygons from startIter until iter
    VertexIterator startIter = tempPoints.begin(); 
    VertexIterator iter;
    for (iter = tempPoints.begin(); iter!=tempPoints.end(); ++iter){
        // If we have the min amount of points required
        if(countOfPoints >= maxPointsPerSet){
            Segment_2 segRight(Point_2(*(iter-2)), Point_2(*(iter-1)));
            Segment_2 segLeft(Point_2(*(iter-1)), Point_2(*(iter)));
            // If it is monotone, update the marked segements and make the polygon
            if((segRight.source().y() < segRight.target().y()) && (segLeft.source().y() > segLeft.target().y())){ 
                vectorOftLowerHullEdge->push_back(pair<Segment_2, Segment_2>(segRight, segLeft));
                vectorOfPolygons->push_back(Polygon_2(startIter, iter));
                startIter = iter-1; // Update starter iterator
                countOfPoints = 0;
                maxPointsPerSet = m-1; //After the first polyogn, update the max points number to m-1
            } 
        }
        countOfPoints++;
    }
    
    // Change the last Polygon
    Polygon_2 temp(startIter, iter);
    if(temp.size() < m/2){ // If the remaining points at the end are sufficiently few, they are appended to the last subset
        vectorOftLowerHullEdge->pop_back();
        int no = vectorOfPolygons->size()-1;
        vectorOfPolygons->at(no).insert(vectorOfPolygons->at(no).end(), startIter+1, iter);
    }else{ // Else, just
        Segment_2 segRight(Point_2(*(iter-2)), Point_2(*(iter-1)));
        Segment_2 segLeft(Point_2(*(iter-1)), Point_2(*(iter)));
        vectorOftLowerHullEdge->push_back(pair<Segment_2, Segment_2>(segLeft, segRight));
        vectorOfPolygons->push_back(Polygon_2(startIter, iter));
    }
    return;
}

// Find the points from the rightmost and leftmost segments for simulated annealing-subdivision
Points make_points_to_ignore(int i, int size, vector<pair<Segment_2, Segment_2>> vectorOftLowerHullEdge){
    Points pointsToIgnore; // Vectro with the points
    if(i==0){// If we are at the first polygon, check only right
        if(vectorOftLowerHullEdge.at(i).first.source().x() < vectorOftLowerHullEdge.at(i).first.target().x()){
            pointsToIgnore.push_back(vectorOftLowerHullEdge.at(i).first.source());
            pointsToIgnore.push_back(vectorOftLowerHullEdge.at(i).first.target());
        }else{
            pointsToIgnore.push_back(vectorOftLowerHullEdge.at(i).first.target());
            pointsToIgnore.push_back(vectorOftLowerHullEdge.at(i).first.source());
        }
    }else if(i == size-1){ // If we are at the last polygon, check only left
        if(vectorOftLowerHullEdge.at(i-1).second.source().x() < vectorOftLowerHullEdge.at(i-1).second.target().x()){
            pointsToIgnore.push_back(vectorOftLowerHullEdge.at(i-1).second.source());
            pointsToIgnore.push_back(vectorOftLowerHullEdge.at(i-1).second.target());
        }else{
            pointsToIgnore.push_back(vectorOftLowerHullEdge.at(i-1).second.target());
            pointsToIgnore.push_back(vectorOftLowerHullEdge.at(i-1).second.source());
        }
    }else{ // If we are at any other polygon, check only left and right
        if(vectorOftLowerHullEdge.at(i-1).second.source().x() < vectorOftLowerHullEdge.at(i-1).second.target().x()){
            pointsToIgnore.push_back(vectorOftLowerHullEdge.at(i-1).second.source());
            pointsToIgnore.push_back(vectorOftLowerHullEdge.at(i-1).second.target());
        }else{
            pointsToIgnore.push_back(vectorOftLowerHullEdge.at(i-1).second.target());
            pointsToIgnore.push_back(vectorOftLowerHullEdge.at(i-1).second.source());
        }
        if(vectorOftLowerHullEdge.at(i).second.source().x() < vectorOftLowerHullEdge.at(i).second.target().x()){
            pointsToIgnore.push_back(vectorOftLowerHullEdge.at(i).first.source());
            pointsToIgnore.push_back(vectorOftLowerHullEdge.at(i).first.target());
        }else{
            pointsToIgnore.push_back(vectorOftLowerHullEdge.at(i).first.target());
            pointsToIgnore.push_back(vectorOftLowerHullEdge.at(i).first.source());
        }
    }
    return pointsToIgnore;
}

// Make temp input file for convex Hull/incremental for simulated annealing-subdivision
string make_input_file(int i, vector<Polygon_2> vectorOfPolygons){
    // Make temp Inputfile for convex_hull_algorithm
    string myInputFile = "mytempfile" + to_string(i) + ".txt";
    ofstream outputFileClass;
    outputFileClass.open(myInputFile.c_str());
    outputFileClass  << "# "<< endl;
    outputFileClass  << "# parameters \"convex_hull\": {\"area\": \"" << 1 << "\"} "<< endl;
    int index = 0;
    for (VertexIterator iter = vectorOfPolygons.at(i).begin(); iter!=vectorOfPolygons.at(i).end(); ++iter){
        outputFileClass << index <<"    " <<(int)iter->x() << "    " << (int)iter->y()<< endl;
        index++;
    }
    outputFileClass.close();
    return myInputFile;
}