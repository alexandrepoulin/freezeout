#include "../include/NumericalCalculator.h"
#include <string>

int main(int argc, char* argv[]){

    //files for input and output 
    string in;
    string out;


    in= "C:/Users/AlexandrePoulin/Documents/relic/freezeout_v2.5/inputFiles/inputDataTest.txt";
    out= "C:/Users/AlexandrePoulin/Documents/relic/freezeout_v2.5/outputFiles/outputDataTest2.txt";
    
    //initialize the numerical calculator
    NumericalCalculator numCalc(in,out);
    cout << "Done initializing." << endl;
    

    //go through the solving process and print omega for both
    numCalc.solve(true);
    

    cout << numCalc.omega[0] << endl;
    cout << numCalc.omega[1] << endl;

    return 1;

}
