#include <string>

#include <CGAL/Polygon_2.h>
#include <CGAL/Simple_cartesian.h>

#include "optimization.h"

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;

using namespace std;

polygon_is_ok convex_hull_algorithm(string inputFile, string outputFile, int edgeSelection, double max_time);
polygon_is_ok incremental_algorithm(string inputFile, string outputFile, string initializationType, int edgeSelection, double max_time);
bool isVisible(Segment_2 edge, Point_2 p, Polygon_2 polygon);

// Sort points by acd. x
bool sortPoints(Point_2 A, Point_2 B);

// Compare functions for point sorting in incremental
// if 1st dimension is identical, sort according to the 2nd one, so that there are no cases
// with collinear points where the 3rd point is inbetween the 1st and the 2nd one
bool comparePointsAscendingX(Point_2 p1, Point_2 p2);

bool comparePointsDescendingX(Point_2 p1, Point_2 p2);

bool comparePointsAscendingY(Point_2 p1, Point_2 p2);

bool comparePointsDescendingY(Point_2 p1, Point_2 p2);