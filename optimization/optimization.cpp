#include <iostream>
#include <string>
#include <cstring>

#include "polygonization.h"
#include "optimization.h"

#include <CGAL/Polygon_2.h>

using namespace std;

typedef CGAL::Simple_cartesian<double>  K;
typedef CGAL::Polygon_2<K>              Polygon_2;

int main(int argc, char *argv[]){
    // ------------------------------------------Check Parameters------------------------------------------

    string inputFile, outputFile, algorithmType, initializationAlgorithm, annealingType;
    int L;
    double threshold;
    bool maximization_problem = true;

    for(int i = 1; i < argc; i++){
        if(!strcmp(argv[i],"-i")){
            inputFile = string(argv[i+1]);
        } else if(!strcmp(argv[i],"-o")){
            outputFile = string(argv[i+1]);
        }else if(!strcmp(argv[i],"-algorithm")){
            algorithmType = string(argv[i+1]);
        }else if(!strcmp(argv[i],"-initialization_algorithm")){
            initializationAlgorithm = string(argv[i+1]);
        }else if(!strcmp(argv[i],"-L")){
            L = strtol(argv[i+1], NULL, 10);
        }else if(!strcmp(argv[i],"-min")){
            maximization_problem = false;
        }else if(!strcmp(argv[i],"-threshold")){
            threshold = stod(argv[i+1]);
        }else if(!strcmp(argv[i],"-annealing")){
            annealingType = string(argv[i+1]);
        }
    }

    // if input file does not exist, exit with error
    ifstream in(inputFile.c_str());
    if(!in){
        cout << "Input-File: " << inputFile << " does not exist." <<endl;
        exit(1);
    }

    Polygon_2 init_polygon;
    Polygon_2 polygon;

    std::chrono::_V2::system_clock::time_point started, done;

    // Call appropriate polygonization algorithm, providing the already received inputs
    if (algorithmType == "local_search"){

        if (initializationAlgorithm == "incremental"){
            if (maximization_problem){
                init_polygon = incremental_algorithm(inputFile, "temp.txt", "1a", 3);
            }else{
                init_polygon = incremental_algorithm(inputFile, "temp.txt", "1a", 2);
            }
        }else if (initializationAlgorithm == "convex_hull"){
            if (maximization_problem){
                init_polygon = convex_hull_algorithm(inputFile, "temp.txt", 3);
            }else{
                init_polygon = convex_hull_algorithm(inputFile, "temp.txt", 2);
            }
        }
        
        started = chrono::high_resolution_clock::now();
        polygon = local_search(init_polygon, L, maximization_problem, threshold);
        done = chrono::high_resolution_clock::now();

    }else if (algorithmType == "simulated_annealing"){

        if(annealingType == "subdivision"){
            polygon = simulated_annealing_algorithm_subdivision(inputFile, outputFile, maximization_problem, L, annealingType, 3, initializationAlgorithm);
            return 0;
        }else{
            if(initializationAlgorithm == "convex_hull"){
                if(maximization_problem){
                    init_polygon = convex_hull_algorithm(inputFile, "temp.txt", 3);
                }else{
                    init_polygon = convex_hull_algorithm(inputFile, "temp.txt", 2);
                }
            }else if(initializationAlgorithm == "incremental"){
                if(maximization_problem){
                    init_polygon = incremental_algorithm(inputFile, "temp.txt", "1a", 3);
                }else{
                    init_polygon = incremental_algorithm(inputFile, "temp.txt", "1a", 2);
                }
                
            }   

            started = chrono::high_resolution_clock::now();
            polygon = simulated_annealing_algorithm(inputFile, outputFile, maximization_problem, init_polygon, L, annealingType, Points());
            done = chrono::high_resolution_clock::now();
        } 
    }

    // output results to file
    ofstream outputFileClass(outputFile);
    outputFileClass << "Optimal Area Polygonization" <<endl;

    for (auto vertex = polygon.vertices_begin(); vertex != polygon.vertices_end() ; ++vertex){
        outputFileClass << fixed << *vertex << endl;
    }

    for (auto edge = polygon.edges_begin(); edge != polygon.edges_end() ; ++edge){
        outputFileClass << fixed << *edge << endl;
    }

    outputFileClass << "Algorithm: " << algorithmType << " [" << ((maximization_problem)? "max]": "min]") << endl;
    double initial_polygon_area = init_polygon.area();
    outputFileClass << "Area_initial: " << initial_polygon_area << endl;
    double polygon_area = polygon.area();
    outputFileClass << fixed << "Area: " << polygon_area << endl;
    
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
    outputFileClass << fixed << "Ratio: " << polygon_area/areaFromInput << endl;
    outputFileClass << fixed << "Construction Time: " << chrono::duration_cast<chrono::milliseconds>(done-started).count() << " ms" << endl;

    cout << "Construction Time: " << chrono::duration_cast<chrono::milliseconds>(done-started).count() << " ms" << endl;

    return 0;
}