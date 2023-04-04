#include <string>

#include <CGAL/Polygon_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Kd_tree_rectangle.h>
#include <CGAL/Fuzzy_iso_box.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Point_2                              Point_2;
typedef std::vector<Point_2>                    Points;
typedef K::Segment_2                            Segment_2;
typedef CGAL::Search_traits_2<K>                Traits;
typedef CGAL::Fuzzy_iso_box<Traits>             Fuzzy_iso_box;

using namespace std;

Polygon_2 local_search(Polygon_2 init_polygon, int L, bool max, double threshold);
Polygon_2 modified_convex_hull_algorithm(string inputFile, string outputFile, int edgeSelection, Points pointsToIgnore);
Polygon_2 modified_incremental_algorithm(string inputFile, string outputFile, string initializationType, int edgeSelection, Points pointsToIgnore);
bool isVisible(Segment_2 edge, Point_2 p, Polygon_2 polygon);
Segment_2 findPolygonSegmentContainingPoint(Polygon_2 polygon, Point_2 p);

// Main functions for s.a.
Polygon_2 simulated_annealing_algorithm(string inputFile, string outputFile, bool maximization_problem, Polygon_2 startPolygon, int L, string annealing, Points pointsToIgnore);
int localStep(Polygon_2 startPolygon, Points* alreadyCheckedpoints, Points* setOf4Points);
int globalStep(Polygon_2 startPolygon,vector<pair<Point_2, Point_2>>* alreadyCheckedpointsGlobal, Points* setOf5Points, Points pointsToIgnore);
Polygon_2 simulated_annealing_algorithm_subdivision(string inputFile, string outputFile, bool maximization_problem, int L, string annealing, int edgeSelection, string startAlgorithm);

// Helping functions for s.a.
Fuzzy_iso_box findRectagular(Polygon_2 startPolygon, Point_2 selectedPoint, Point_2* myA, Point_2* myB, Point_2* myC, Point_2* myD);
bool check_validity_for_local_step(Segment_2 pr, Segment_2  qs, Point_2 AA, Point_2 DD, Polygon_2 startPolygon, list<Point_2> kdSearchRectPoints);
void outputForSA(string inputFile, string outputFile, Polygon_2 startPolygon, bool maximization_problem, double initial_polygon_area, std::chrono::_V2::system_clock::time_point started, std::chrono::_V2::system_clock::time_point done, string annealing);
double newAreaForTwoTriangles(string annealing, Points setOfPointsReturned,Polygon_2 startPolygon);
void updatePolygon(Polygon_2* startPolygon, Points setOfPointsReturned,string annealing);
double getAreaOfConvexHull(Polygon_2 startPolygon);
bool check_validity_for_global_step(Segment_2 pr, Segment_2 sq, Segment_2 qt, Polygon_2 startPolygon, Point_2 P, Point_2 R, Point_2 S, Point_2 Q, Point_2 T);
void make_small_polygons_and_mark_segments(string inputFile, vector<Polygon_2>* vectorOfPolygons, vector<pair<Segment_2, Segment_2>>* vectorOftLowerHullEdge);
Points make_points_to_ignore(int i, int size, vector<pair<Segment_2, Segment_2>> vectorOftLowerHullEdge);
string make_input_file(int i, vector<Polygon_2> vectorOfPolygons);