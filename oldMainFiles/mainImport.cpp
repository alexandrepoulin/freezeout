#include "../include/NumericalCalculator.h"

int main(){

    //files for input and output 
    string in;
    string out;
    
    in= "C:/Users/AlexandrePoulin/Documents/relic/freezeout_v2.2/inputFiles/inputData.txt";
    out= "C:/Users/AlexandrePoulin/Documents/relic/freezeout_v2.2/outputFiles/outputData.txt";
    
    cout << "Starting" << endl;
    
    //initialize the numerical calculator
    NumericalCalculator numCalc(in,out);
    cout << "Done initializing." << endl;

    //go through the solving process and print omega for both
    numCalc.solve(true);
    cout << "Omega 0 = " << numCalc.omega[0] << endl;
    //cout << "Omega 1 = " << numCalc.omega[1] << endl;
    //cout << "Omega 2 = " << numCalc.omega[2] << endl;
    //cout << "Omega 3 = " << numCalc.omega[3] << endl;
    //cout << "Omega 4 = " << numCalc.omega[4] << endl;

    return 1;

}
