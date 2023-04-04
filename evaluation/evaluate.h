#include <string>
#include <vector>

using namespace std;

double getCHAreaFromInputFile(string inputFile);
int getNumberOfPointsFromInputFile(string inputFile);
long long int calculateThreshold(double area);
void init_parameters_vector(vector<vector<string>>& parameters_vector);
void preprocessing(string inputFileInDir, double expInitCutOff, double expCutOff, long int threshold, vector<vector<string>>& parameters_vector);
void algorithms_comparison(vector<double>& scores, string inputFileInDir, vector<vector<string>> parameters_vector, long int cutOff, long int threshold, double convexHullArea);