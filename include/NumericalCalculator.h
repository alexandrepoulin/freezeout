#include <boost/math/special_functions/bessel.hpp>
#include <time.h>

#include "YGReader.h"
#include "InputReader.h"
#include "OutputWriter.h"


using namespace std;

//This object will be responsible to knowing how to do all of the numerical methods.
class NumericalCalculator
{
public:
    
    //Constructor
    NumericalCalculator();
    NumericalCalculator(string in, string out, int nDM=-1);
    //Detructor
    ~NumericalCalculator();

    //Let's you change the input and output files
    void newInput(string in);
    void newOutput(string out);

    //get the values of mu as a function of temperature
    double getMu(double t, int ind);
    //generates the Yeq if not using pre calculated
    double Yeq(double m, double mu, double T, double gstarS, double g);
    vector<double> genEqSet(double t);
    double genEqRatio(int i, int j, int k, int l, double t, double inI, double inJ);

    //This function gives you the Yi' from the equation Yi'=f(x,Y1,Y2,..Yn)
    //This is where ALL of the physics is located
    vector<double> f(vector<double> inset,vector<double> eqset,double t);

    //gives you the value for where you will start solving
    //vector<double> startingX();

    //The two methods for solving the differential equation
    void fnEuler(vector<double>& inset, vector<double> eqset,double t);
    void fnRK4(vector<double>& inset, vector<double> eqset,double t);
    void fnLMM(vector<double>& inset, vector<double> eqset, vector<vector<double>* >& prevF, int& prevFIndex, double t);
  
    //goes through and solves using the parameters in the inputReader
    //write determins whether to output text or write to a file (in case
    //you just want omega and will call it alot)
    void solve(bool write);

    //some parameters that popup 
    const double CMCBSMTOPB = 3.336* pow(10,25); //converts cm^3/s to pb
    const double YeqPref = 45.0/(4.0 * pow( 3.14159265358979,4));
    const double Yrad = 2.0 * 1.2020569 * YeqPref ;
    const double rhoC = 1.055 * pow(10,-5);
    const double s0 = 2970;

    //The list of resulting omega from solve()
    vector<double> omega;
    
    //io-streams
    InputReader* inputReader;
    OutputWriter* outputWriter;
    YGReader ygReader;
    

protected:

    //parameter that pops up in startingX()
    //const double c = 0.001;

};
