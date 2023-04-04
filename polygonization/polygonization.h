#include <string>

#include <CGAL/Polygon_2.h>
#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;

using namespace std;

int convex_hull_algorithm(string inputFile, string outputFile, int edgeSelection);
int incremental_algorithm(string inputFile, string outputFile, string initializationType, int edgeSelection);
bool isVisible(Segment_2 edge, Point_2 p, Polygon_2 polygon);
Segment_2 findPolygonSegmentContainingPoint(Polygon_2 polygon, Point_2 p);