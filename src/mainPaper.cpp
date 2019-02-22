#include "../include/NumericalCalculator.h"
#include <string>
#include <stdlib.h>

int main(int argc, char* argv[]){

    //files for input and output 
    string in;
    string out;


    /*
     * argv[1] = hardcoded smallest x
     * argv[2] = file id
     * argv[3] = single case
     * argv[4] = ndm
     * argv[n] = mxset[n-5]
     * */


    string id = argv[2];
    
    
    in= "C:/Users/AlexandrePoulin/Documents/relic/freezeout_v2.5/inputFiles/inputData"+id+".txt";
    out= "C:/Users/AlexandrePoulin/Documents/relic/freezeout_v2.5/outputFiles/outputOmega"+id +".txt";
   
    //initialize the numerical calculator
    NumericalCalculator numCalc(in,out);
    cout << "Done initializing." << endl;
    
    float maxMass = 0;

    for(int ind=5; ind<5+atoi(argv[4]); ind++){
        numCalc.inputReader->mxset[ind-5] = atof(argv[ind]);
        if(atof(argv[ind])>maxMass) maxMass = atof(argv[ind]);
    }

    numCalc.inputReader->changeMassScale(maxMass);
    
    numCalc.inputReader->changeSmallestX(atof(argv[1]));

    bool singleCase = atoi(argv[3]);
    if(singleCase){
        numCalc.solve(singleCase);
        for(int i=0; i<1; i++) cout << numCalc.omega[i] << endl;
    }
    else{
        numCalc.solve(false);
        numCalc.outputWriter->writeOmega(numCalc.omega);
    }
    return 1;

}
