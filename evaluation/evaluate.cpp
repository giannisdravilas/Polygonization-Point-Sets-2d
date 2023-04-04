#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <iomanip>
#include <dirent.h>
#include <cmath>
#include <limits>
#include <list>
#include <map>
#include <float.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>

#include "polygonization.h"
#include "optimization.h"
#include "evaluate.h"

typedef CGAL::Simple_cartesian<double>  K;
typedef CGAL::Polygon_2<K>              Polygon_2;

using namespace std;

int main(int argc, char *argv[]){

    // ------------------------------------------Check Parameters------------------------------------------
    string inputDir, outputFile;
    bool preprocess = false;
    for(int i = 1; i < argc; i++){
        if(!strcmp(argv[i],"-i")){
            inputDir = string(argv[i+1]);
        } else if(!strcmp(argv[i],"-o")){
            outputFile = string(argv[i+1]);
        }else if(!strcmp(argv[i],"-preprocess")){
            preprocess = true;
        }
    }
    cout <<"InputDir: " << inputDir << endl;
    cout <<"OutputFile: " << outputFile << endl;
    cout <<"Preprocess: " << preprocess << endl;

    // Open Dir
    struct dirent *entry;
    DIR *dp = opendir(inputDir.c_str());
    if (dp == NULL) {
        perror("opendir: Path does not exist or could not be read.");
        return -1;
    }
    
    // Map[Key: number of points, Value: vector of 20 scores (4 for each algorithm)]
    // scores layout is the following:
    //  min score | max score | min bound | max bound -> repeated 5 times, one for each algorithm
    map<int, vector<double>> scoresPerSize;
    int count_of_files = 0; // Counter

    // For every file in dir
    while ((entry = readdir(dp))){

        // Skip these file names
        if( (strcmp(entry->d_name, ".") == 0) || (strcmp(entry->d_name, "..") == 0) ){
            continue;
        }
        count_of_files++;

        // Get the name of the File
        string inputFileInDir = inputDir + "/" + entry->d_name;
        cout << endl << "=================================================================================================" << endl;
        cout << "(" << count_of_files << ") File-name: " << inputFileInDir << endl;
        cout << "=================================================================================================" << endl << endl;
        
        // Insert in map
        int noOfPoints = getNumberOfPointsFromInputFile(inputFileInDir);
        if(scoresPerSize.find(noOfPoints) == scoresPerSize.end()){ // If we can't find key in map, insert new element
            // Initialize vector with 20 places containing 0s (5 algorithms * 4 scores), except for the 4th place
            // of each algorithm, where max bound is initialized to double max.
            scoresPerSize[noOfPoints] = vector<double>(20, 0.0);
            for(int i=3; i < 20; i = i+4){
                scoresPerSize[noOfPoints][i] = DBL_MAX;
            }
        }
        vector<double> scores = scoresPerSize[noOfPoints];

        // Get Cut-Off
        long int cutOff = 500 * noOfPoints;

        // Cut-Off for Experiments in Preprocessing (total 10% of cutOff)
        double expInitCutOff = 0.05 * cutOff; expInitCutOff /= 4;
        double expCutOff = 0.05 * cutOff; expCutOff /= 92;

        // Get Area from Input File
        double convexHullArea = getCHAreaFromInputFile(inputFileInDir);

        // Calculate threshold
        long long int threshold = calculateThreshold(convexHullArea);
        
        // Make vector of parameters (2 dimensions, initalize 2nd dimension with 3 void positions)
        vector<vector<string>> parameters_vector;
        for(int i = 0; i < 10; i++){
            vector<string> temp; 
            temp.push_back(" ");temp.push_back(" ");temp.push_back(" ");
            parameters_vector.push_back(temp);
        }

        // Initialize parameters of algorithms' combinations
        init_parameters_vector(parameters_vector);

        // If we want preprocessing to be done for hyperparameters
        if(preprocess){
            preprocessing(inputFileInDir, expInitCutOff, expCutOff, threshold, parameters_vector);
        }

        // Call algorithms for final comparisons and update scores vector.
        algorithms_comparison(scores, inputFileInDir, parameters_vector, cutOff, threshold, convexHullArea);
        scoresPerSize[noOfPoints] = scores;
    }
    closedir(dp); // Close Dir 

    /////////////////////////////////////////////////////////////////////////
    ///////////////////////////   OUTPUT    /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    ofstream outfile;
    outfile.open(outputFile.c_str());
    int number_of_algs = 5;
    string algo_names[number_of_algs] = {"Local Search (Incremental)", "Simulated Annealing (Incremental)", "Local Search (Convex Hull)", "Simulated Annealing (Convex Hull)", "Simulated Annealing Subdivision"};
    
    // Headers
    outfile << "      ";
    for(int i = 0; i < number_of_algs; i++){
        outfile << "|| " << left << setw(46) << algo_names[i];
    }
    outfile << "||" << endl <<  "Size  ";
    for(int i = 0; i < number_of_algs; i++){
        outfile << setw(49) << "|| min score | max score | min bound | max bound";
    }
    outfile << "||" << endl;

    // Print Map with actual score values
    for (auto const& x : scoresPerSize){
        int noOfPoints = x.first;
        vector<double> scores = x.second;

        outfile << left << setw(6) << noOfPoints;

        for(int i=0;i<number_of_algs;i++){
            outfile <<  "|| " << left << setw(9) << scores[i*4] << " | " << left << setw(9) << scores[i*4 + 1] << " | " << left << setw(9) << scores[i*4 + 2] << " | " << left << setw(9) << scores[i*4 + 3] << " ";
        }
        outfile << "||" << endl;
    }
    return 0;
}

double getCHAreaFromInputFile(string inputFile){
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
            break;
        }
    }
    return areaFromInput;
}

int getNumberOfPointsFromInputFile(string inputFile){
    // Read From File
    ifstream in(inputFile.c_str());
    string line;
    int numberOfPoints = 0;
    while (getline(in, line)){ // Add points to vector
        istringstream iss(line);
        numberOfPoints++;
    }
    return numberOfPoints - 2;
}

// Threshold after which we assume further optimization would make no more difference
// Threshold is calculated as 2 orders of magnitude less than the area
long long int calculateThreshold(double area){
    int i = -2;
    while(area > 9){
        area = area / 10;
        i++;
    }
    return pow(10, i);
}

// initalize vector of parameters
void init_parameters_vector(vector<vector<string>>& parameters_vector){
    // Init parameters
    parameters_vector[0][0] = "1a";
    parameters_vector[0][1] = "5";
    parameters_vector[1][0] = "1a";
    parameters_vector[1][1] = "5000";
    parameters_vector[1][2] = "global";
    parameters_vector[2][0] = "5";
    parameters_vector[3][0] = "5000";
    parameters_vector[3][1] = "global";
    parameters_vector[4][0] = "1a";
    parameters_vector[4][1] = "5";
    parameters_vector[5][0] = "1a";
    parameters_vector[5][1] = "5000";
    parameters_vector[5][2] = "global";
    parameters_vector[6][0] = "5";
    parameters_vector[7][0] = "5000";
    parameters_vector[7][1] = "global";
    parameters_vector[8][0] = "5000";
    parameters_vector[8][1] = "global";
    parameters_vector[8][2] = "incremental";
    parameters_vector[9][0] = "5000";
    parameters_vector[9][1] = "global";
    parameters_vector[9][2] = "incremental";
}

// preprocessing for selection of best hyperparameters
// we run each combination of algorithms and hyperparameters only for a small amount of time and decide which is the best
// combination based on these preliminary results
void preprocessing(string inputFileInDir, double expInitCutOff, double expCutOff, long int threshold, vector<vector<string>>& parameters_vector){
    
    // Parameters
    string incrementalOptions[4] = {"1a", "1b", "2a", "2b"};
    string annealingT[2] = { "local", "global"};
    int annealingL[2] = {5000, 10000};
    int localSearchL[2] = {5, 10};
    double best_min_area_1 = DBL_MAX, best_min_area_2 = DBL_MAX;
    double best_max_area_1 = 0, best_max_area_2 = 0;

    Polygon_2 polygon;

    // MIN - incremental
    for(int i = 0; i < 4; i++){

        // Incremental
        polygon_is_ok init_polygon = incremental_algorithm(inputFileInDir, "temp.txt", incrementalOptions[i], 2, expInitCutOff);
        if(init_polygon.isOk){

            // Local_search
            for(int j = 0; j < 2; j++){
                polygon_is_ok returned_polygon = local_search(init_polygon.polygon, localSearchL[j], false, threshold, expCutOff);
                polygon = returned_polygon.polygon;

                // Update parameters
                if (polygon.area() < best_min_area_1){
                    best_min_area_1 = polygon.area();
                    parameters_vector[0][0] = incrementalOptions[i];
                    parameters_vector[0][1] = to_string(localSearchL[j]);
                }
            }

            
            // Simulated_annealing
            for(int j = 0; j < 2; j++){
                for(int k = 0; k < 2; k++){
                    polygon_is_ok returned_polygon = simulated_annealing_algorithm(inputFileInDir, "temp.txt", false, init_polygon.polygon, annealingL[j], annealingT[k], Points(), expCutOff);
                    polygon = returned_polygon.polygon;  

                    // Update parameters
                    if (polygon.area() < best_min_area_2){
                        best_min_area_2 = polygon.area();
                        parameters_vector[1][0] = incrementalOptions[i];
                        parameters_vector[1][1] = to_string(annealingL[j]);
                        parameters_vector[1][2] = annealingT[k];
                    }
                }
            }
        }
    }

    best_min_area_1 = DBL_MAX;
    best_min_area_2 = DBL_MAX;

    // MIN - convex_hull
    polygon_is_ok init_polygon = convex_hull_algorithm(inputFileInDir, "temp.txt", 2, expInitCutOff);
    if(init_polygon.isOk){

        // Local_search
        for(int j = 0; j < 2; j++){
            polygon_is_ok returned_polygon = local_search(init_polygon.polygon, localSearchL[j], false, threshold, expCutOff);
            polygon = returned_polygon.polygon; 

            // Update parameters
            if (polygon.area() < best_min_area_1){
                best_min_area_1 = polygon.area();
                parameters_vector[2][0] = to_string(localSearchL[j]);
            }
                
        }

        // Simulated_annealing
        for(int j = 0; j < 2; j++){
            for(int k = 0; k < 2; k++){
                polygon_is_ok returned_polygon = simulated_annealing_algorithm(inputFileInDir, "temp.txt", false, init_polygon.polygon, annealingL[j], annealingT[k], Points(), expCutOff);
                polygon = returned_polygon.polygon;  

                // Update parameters
                if (polygon.area() < best_min_area_2){
                    best_min_area_2 = polygon.area();
                    parameters_vector[3][0] = to_string(annealingL[j]);
                    parameters_vector[3][1] = annealingT[k];
                }   
            }
        }
    }
    
    // MAX - incremental
    for(int i = 0; i < 4; i++){

        // Incremental
        polygon_is_ok init_polygon = incremental_algorithm(inputFileInDir, "temp.txt", incrementalOptions[i], 3, expInitCutOff);
        if(init_polygon.isOk){

            // Local_search
            for(int j = 0; j < 2; j++){
                polygon_is_ok returned_polygon = local_search(init_polygon.polygon, localSearchL[j], true, threshold, expCutOff);
                polygon = returned_polygon.polygon;  

                // Update parameters
                if (polygon.area() > best_max_area_1){
                    best_max_area_1 = polygon.area();
                    parameters_vector[4][0] = incrementalOptions[i];
                    parameters_vector[4][1] = to_string(localSearchL[j]);
                }
            }

            // Simulated_annealing
            for(int j = 0; j < 2; j++){
                for(int k = 0; k < 2; k++){
                    polygon_is_ok returned_polygon = simulated_annealing_algorithm(inputFileInDir, "temp.txt", true, init_polygon.polygon, annealingL[j], annealingT[k], Points(), expCutOff);
                    polygon = returned_polygon.polygon;  

                    // Update parameters
                    if (polygon.area() > best_max_area_2){
                        best_max_area_2 = polygon.area();
                        parameters_vector[5][0] = incrementalOptions[i];
                        parameters_vector[5][1] = to_string(annealingL[j]);
                        parameters_vector[5][2] = annealingT[k];
                    }   
                }
            }
        } 
    }

    best_max_area_1 = 0;
    best_max_area_2 = 0;

    // MAX - convex_hull
    init_polygon = convex_hull_algorithm(inputFileInDir, "temp.txt", 3, expInitCutOff);
    if(init_polygon.isOk){

        // Local_search
        for(int j = 0; j < 2; j++){
            polygon_is_ok returned_polygon = local_search(init_polygon.polygon, localSearchL[j], true, threshold, expCutOff);
            polygon = returned_polygon.polygon; 

            // Update parameters
            if (polygon.area() > best_max_area_1){
                best_max_area_1 = polygon.area();
                parameters_vector[6][0] = to_string(localSearchL[j]);
            }   
        }
        
        // Simulated_annealing
        for(int j = 0; j < 2; j++){
            for(int k = 0; k < 2; k++){
                polygon_is_ok returned_polygon = simulated_annealing_algorithm(inputFileInDir, "temp.txt", true, init_polygon.polygon, annealingL[j], annealingT[k], Points(), expCutOff);
                polygon = returned_polygon.polygon; 

                // Update parameters
                if (polygon.area() > best_max_area_2){
                    best_max_area_2 = polygon.area();
                    parameters_vector[7][0] = to_string(annealingL[j]);
                    parameters_vector[7][1] = annealingT[k];
                }      
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////
    ////////////////// SIMULATED ANNEALING SUBDIVISION //////////////////////
    /////////////////////////////////////////////////////////////////////////

    string annealingTypeSubdivision[2] = {"incremental", "convex_hull"};

    best_min_area_1 = DBL_MAX;

    // MIN simulated_annealing_subdivision
    for(int j = 0; j < 2; j++){
        for(int k = 0; k < 2; k++){
            for(int l = 0; l < 2; l++){

                polygon_is_ok returned_polygon = simulated_annealing_algorithm_subdivision(inputFileInDir, "temp.txt", false, annealingL[j], annealingT[k], 2, annealingTypeSubdivision[l], expCutOff);
                polygon = returned_polygon.polygon;  

                // Update parameters
                if (polygon.area() < best_min_area_1){
                    best_min_area_1 = polygon.area();
                    parameters_vector[8][0] = to_string(annealingL[j]);
                    parameters_vector[8][1] = annealingT[k];
                    parameters_vector[8][2] = annealingTypeSubdivision[l];
                }
            }
        }
    }

    best_max_area_2 = 0;

    // MAX simulated_annealing_subdivision 
    for(int j = 0; j < 2; j++){
        for(int k = 0; k < 2; k++){
            for(int l = 0; l < 2; l++){
                polygon_is_ok returned_polygon = simulated_annealing_algorithm_subdivision(inputFileInDir, "temp.txt", true, annealingL[j], annealingT[k], 3, annealingTypeSubdivision[l], expCutOff);
                polygon = returned_polygon.polygon;  

                // Update parameters
                if (polygon.area() > best_max_area_2){
                    best_max_area_2 = polygon.area();
                    parameters_vector[9][0] = to_string(annealingL[j]);
                    parameters_vector[9][1] = annealingT[k];
                    parameters_vector[9][2] = annealingTypeSubdivision[l];
                }  
            }
        } 
    }
}

void algorithms_comparison(vector<double>& scores, string inputFileInDir, vector<vector<string>> parameters_vector, long int cutOff, long int threshold, double convexHullArea){
    
    /////////////////////////////////////////////////////////////////////////
    /////////////////////// COMPARE ALGORITHMS MIN //////////////////////////
    /////////////////////////////////////////////////////////////////////////

    polygon_is_ok init_polygon, polygon;
    double algorithm_score;
    
    // incremental_algorithm - local_search
    init_polygon = incremental_algorithm(inputFileInDir, "temp.txt", parameters_vector[0][0], 2, cutOff);
    if(init_polygon.isOk){
        polygon = local_search(init_polygon.polygon, stoi(parameters_vector[0][1]), false, threshold, cutOff); 
        if(polygon.isOk){
            algorithm_score = polygon.polygon.area() / convexHullArea;
        }else{
            algorithm_score = 1;
        }
    }else{
        algorithm_score = 1;
    }

    // Update scores vector
    scores[0] += algorithm_score;
    if(algorithm_score > scores[2]){
        scores[2] = algorithm_score;
    }

    // incremental_algorithm - simulated_annealing_algorithm
    init_polygon = incremental_algorithm(inputFileInDir, "temp.txt", parameters_vector[1][0], 2, cutOff);
    if(init_polygon.isOk){
        polygon = simulated_annealing_algorithm(inputFileInDir, "temp.txt", false, init_polygon.polygon, stoi(parameters_vector[1][1]), parameters_vector[1][2], Points(), cutOff);
        if(polygon.isOk){
            algorithm_score = polygon.polygon.area() / convexHullArea;
        }else{
            algorithm_score = 1;
        }
    }else{
        algorithm_score = 1;
    }

    // Update scores vector
    scores[4] += algorithm_score;
    if(algorithm_score > scores[6]){
        scores[6] = algorithm_score;
    }

    // convex_hull_algorithm - local_search
    init_polygon = convex_hull_algorithm(inputFileInDir, "temp.txt", 2, cutOff);
    if(init_polygon.isOk){
        polygon = local_search(init_polygon.polygon, stoi(parameters_vector[2][0]), false, threshold, cutOff); 
        if(polygon.isOk){
            algorithm_score = polygon.polygon.area() / convexHullArea;
        }else{
            algorithm_score = 1;
        }
    }else{
        algorithm_score = 1;
    }

    // Update scores vector
    scores[8] += algorithm_score;
    if(algorithm_score > scores[10]){
        scores[10] = algorithm_score;
    }

    // convex_hull_algorithm - simulated_annealing_algorithm
    init_polygon = convex_hull_algorithm(inputFileInDir, "temp.txt", 2, cutOff);
    if(init_polygon.isOk){        
        polygon = simulated_annealing_algorithm(inputFileInDir, "temp.txt", false, init_polygon.polygon, stoi(parameters_vector[3][0]), parameters_vector[3][1], Points(), cutOff);
        if(polygon.isOk){
            algorithm_score = polygon.polygon.area() / convexHullArea;
        }else{
            algorithm_score = 1;
        }
    }else{
        algorithm_score = 1;
    }

    // Update scores vector
    scores[12] += algorithm_score;
    if(algorithm_score > scores[14]){
        scores[14] = algorithm_score;
    }

    // simulated_annealing_algorithm_subdivision
    polygon = simulated_annealing_algorithm_subdivision(inputFileInDir, "temp.txt", false, stoi(parameters_vector[8][0]), parameters_vector[8][1], 2, parameters_vector[8][2], cutOff);
    if(polygon.isOk){
        algorithm_score = polygon.polygon.area() / convexHullArea;
    }else{
        algorithm_score = 1;
    }

    // Update scores vector
    scores[16] += algorithm_score;
    if(algorithm_score > scores[18]){
        scores[18] = algorithm_score;
    }
    

    /////////////////////////////////////////////////////////////////////////
    //////////////////////// COMPARE ALGORITHMS MAX /////////////////////////
    /////////////////////////////////////////////////////////////////////////

    // incremental_algorithm - local_search
    init_polygon = incremental_algorithm(inputFileInDir, "temp.txt", parameters_vector[4][0], 3, cutOff);
    if(init_polygon.isOk){
        polygon = local_search(init_polygon.polygon, stoi(parameters_vector[4][1]), true, threshold, cutOff); 
        if(polygon.isOk){
            algorithm_score = polygon.polygon.area() / convexHullArea;
        }else{
            algorithm_score = 0;
        }
    }else{
        algorithm_score = 0;
    }

    // Update scores vector
    scores[1] += algorithm_score;
    if(algorithm_score < scores[3]){
        scores[3] = algorithm_score;
    }

    // incremental_algorithm - simulated_annealing_algorithm
    init_polygon = incremental_algorithm(inputFileInDir, "temp.txt", parameters_vector[5][0], 3, cutOff);
    if(init_polygon.isOk){
        polygon = simulated_annealing_algorithm(inputFileInDir, "temp.txt", true, init_polygon.polygon, stoi(parameters_vector[5][1]), parameters_vector[5][2], Points(), cutOff);
        if(polygon.isOk){
            algorithm_score = polygon.polygon.area() / convexHullArea;
        }else{
            algorithm_score = 0;
        }
    }else{
        algorithm_score = 0;
    }

    // Update scores vector
    scores[5] += algorithm_score;
    if(algorithm_score < scores[7]){
        scores[7] = algorithm_score;
    }

    // convex_hull_algorithm - local_search
    init_polygon = convex_hull_algorithm(inputFileInDir, "temp.txt", 3, cutOff);
    if(init_polygon.isOk){
        polygon = local_search(init_polygon.polygon, stoi(parameters_vector[6][0]), true, threshold, cutOff); 
        if(polygon.isOk){
            algorithm_score = polygon.polygon.area() / convexHullArea;
        }else{
            algorithm_score = 0;
        }
    }else{
        algorithm_score = 0;
    }

    // Update scores vector
    scores[9] += algorithm_score;
    if(algorithm_score < scores[11]){
        scores[11] = algorithm_score;
    }

    // convex_hull_algorithm - simulated_annealing_algorithm
    init_polygon = convex_hull_algorithm(inputFileInDir, "temp.txt", 3, cutOff);
    if(init_polygon.isOk){
        polygon = simulated_annealing_algorithm(inputFileInDir, "temp.txt", true, init_polygon.polygon, stoi(parameters_vector[7][0]), parameters_vector[7][1], Points(), cutOff); 
        if(polygon.isOk){
            algorithm_score = polygon.polygon.area() / convexHullArea;
        }else{
            algorithm_score = 0;
        }
    }else{
        algorithm_score = 0;
    }

    // Update scores vector
    scores[13] += algorithm_score;
    if(algorithm_score < scores[15]){
        scores[15] = algorithm_score;
    }   

    // simulated_annealing_algorithm_subdivision
    polygon = simulated_annealing_algorithm_subdivision(inputFileInDir, "temp.txt", true, stoi(parameters_vector[9][0]), parameters_vector[9][1], 2, parameters_vector[9][2], cutOff);
    if(polygon.isOk){
        algorithm_score = polygon.polygon.area() / convexHullArea;
    }else{
        algorithm_score = 0;
    }

    // Update scores vector
    scores[17] += algorithm_score;
    if(algorithm_score < scores[19]){
        scores[19] = algorithm_score;
    }
}