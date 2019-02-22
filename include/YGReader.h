#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <vector>
#include <iterator>

using namespace std;

//This class imports Yeq, gstar and gstarS
class YGReader
{

public:
    //constructor
    YGReader();
    //detructor
    ~YGReader(){};
   
    //These let you get the data
    double getYeq(double x);
    double getgstar(double T);
    double getgstarS(double T);
    double getgstarSDer(double T);

    //Reads through the files
    void read(bool impY);


protected:

    //index functions
    int YeqPreCalcInd(double x);
    int gindex(double T);
    
    //lenght of YeqPreCalc
    int lenPreCalc;
    
    //actual data
    vector<double> gstarData;
    vector<double> gstarSData;
    vector<double> gstarSDerData;
    vector<double> YeqPreCalc;
   
    //input streams
    ifstream inputg;
    ifstream inputgs;
    ifstream inputgsDer;
    ifstream inputY;

};
