#include "../include/NumericalCalculator.h"
#include <string>
#include <stdlib.h>


//argv= path, masses, g values, step Small, step Big, file number
int main(int argc, char* argv[]){

    //files for input and output 
    string in;
    string out;

    int ndm= int((argc-5)/2);
    
    string fileNumber(argv[2*ndm+3]);
    cout << fileNumber << endl;
    
    in= "C:/Users/AlexandrePoulin/Documents/relic/freezeout_v2.5/inputFiles/inputData_hidden" + fileNumber+ ".txt";
    out= "C:/Users/AlexandrePoulin/Documents/relic/freezeout_v2.5/outputFiles/outputOmega"  + fileNumber +".txt";
   
    //initialize the numerical calculator
    NumericalCalculator numCalc(in,out,ndm);
    cout << "Done initializing." << endl;

    double biggestMass = 0;
    
    for(int i = 0; i< ndm ;i++){ 
        if(atof(argv[i+1])>biggestMass) biggestMass=atof(argv[i+1]);

        cout << "i=" << i << "; mass: " << atof(argv[i+1])  << "; g:" << atof(argv[i+1+ndm]) << endl;
        numCalc.inputReader->mxset[i]=atof(argv[i+1]);
        numCalc.inputReader->gset[i]=atof(argv[i+1+ndm]);
    }

    numCalc.inputReader->step1=strtod(argv[2*ndm+1],NULL);
    numCalc.inputReader->step2=strtod(argv[2*ndm+2],NULL);
    numCalc.inputReader->changeMassScale(biggestMass);

    bool singleCase = bool(atoi(argv[2*ndm+4]));

    if(singleCase){
        //go through the solving process and print omega for both
        numCalc.solve(true);
        for(int i=0; i<ndm; i++) cout << numCalc.omega[i] << endl;
 
    }
    else{
        //go through the solving process and print omega for both
        numCalc.solve(false);
        numCalc.outputWriter->writeOmega(numCalc.omega);
    }
    return 1;

}
