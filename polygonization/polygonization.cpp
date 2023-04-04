#include <iostream>
#include <string>
#include <cstring>
#include "polygonization.h"

using namespace std;

int main(int argc, char *argv[]){
    // ------------------------------------------Check Parameters------------------------------------------

    string inputFile, outputFile, algorithmType, initializationType;
    int edgeSelection, onionInitialization;

    for(int i = 1; i < argc; i++){
        if(!strcmp(argv[i],"-i")){
            inputFile = string(argv[i+1]);
        } else if(!strcmp(argv[i],"-ο")){
            outputFile = string(argv[i+1]);
        }else if(!strcmp(argv[i],"-algorithm")){
            algorithmType = string(argv[i+1]);
        }else if(!strcmp(argv[i],"-edge_selection")){
            edgeSelection = strtol(argv[i+1], NULL,10);
        }else if(!strcmp(argv[i],"-initialization")){
            initializationType = string(argv[i+1]);
        }else if(!strcmp(argv[i],"-onion_initialization")){
            onionInitialization = strtol(argv[i+1], NULL,10);
        }
    }

    // Call appropriate polygonization algorithm, providing the already received inputs
    if (algorithmType == "incremental"){
        incremental_algorithm(inputFile, outputFile, initializationType, edgeSelection);
    }else if (algorithmType == "convex_hull"){
        convex_hull_algorithm(inputFile, outputFile, edgeSelection);
    }

    return 0;
}