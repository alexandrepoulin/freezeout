#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <set>
#include <math.h>
#include <algorithm>

#include "../include/ImportedSigma.h"

using namespace std;

//This object will be responsible for keeping track of the input from the input file
class InputReader
{

public:
    //Constructor
    InputReader(){};
    InputReader(string in, int nDM=-1);
    //Detructor
    ~InputReader(){};

    //lets you change the file you are working with
    void changeFile(string in, int nDM=-1);

    //some functions to help modify sigma v without having to import a new file
    //this is good if you want to loop alot
    void setZeroSigmav();
    void setSigvijkl(int i,int j, int k, int l, double val);
    void setSigvijklp(int i,int j, int k, int l, double val);

    //Get the cross section (regardless of s or p wave or neither)
    double getSigmav(double t, int ind);
    double getSigmav(double t, int i, int j, int k, int l);
   
    //Get the decay widths 
    double getGamma(int ind);
    double getGamma(int i, int k, int l);

    //Lets you change the massScale
    void changeMassScale(double m);

    //get the correct step size for the value of t
    double step(double t);

    // change the smallest possible starting X
    void changeSmallestX(double x);

    //find the total number of steps to take
    void findNSteps();

    //variables

    //options during the procedure
    bool usePreCalc; // use precalculated Yeq values
    bool waves; // use wave approximation (s and p wave)
    bool useW; // solve in terms of W=ln(Y)
    int method; //0=euler, 1=RK4, 2=LMM

    //some numerical constants which are helpful
    int Nx;
    int Nf;
    int Nf2;
    int Nf3;
    int Nf4;
    int Nsteps;
    const double PI = 3.14159265358979;
    const double PBTOGEVM2= 2.5680*pow(10,-9);
    const double MP = 1.2209 * pow(10,19);

    //some parameters that show up
    double massScale;
    double annPref;
    double step1; // used for x < stepChange
    double step2; // used for x > stepChange
    double smallestX; // the smallest we are willing to start

    //some lists
    vector<double> mxset; //mass list
    vector<double> nDiffSet; //chemical potential list
    vector<double> gset; //g list (degrees of freedom)
    vector<double> sigvijkl; //s wave cross sections
    vector<double> sigvijklp; //p wave cross sections
    vector<double> gammaijk; //decay widths
    set<int> nonZeroIndices; //indices which have a non-zero s or p wave cross section
    set<int> nonZeroIndicesGam; //indices which have a non-zero decay width
    string filename; //file name of inputData.txt
    string filenameInSig; //file name of imported sigma
    ImportedSigma inSig; //the imported sigma

    vector<vector<double>* > intervals; //Set of intervals in which to use the smaller step size
    vector<double> startingxVec;


protected:
 
    //setup the intervals for the step function. Pass in the starting x values
    void setupIntervals();

    //find the starting x values
    void startingX();   
    
    //Reads through the file
    void read(int nDM);
    
    //variables

    //this is the input stream
    ifstream input;

    //Used in startingX function
    double c= 0.001;
    

};
