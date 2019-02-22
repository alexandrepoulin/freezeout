#include <iostream>
#include <string>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <algorithm>

using namespace std;

//This class is responsible for writing the output files.
class OutputWriter
{

public:
    //constructor
    OutputWriter(){};
    OutputWriter(string out,int nx);
    //detructor
    ~OutputWriter(){};

    //Closes the file if open
    void close();

    //writes to the file if open
    void write(double temp, vector<double> eqset, vector<double> inset, bool useW);
    
    //writes the omega values to the file if open
    void writeOmega(vector<double> omega);

    //Keeps track of the last valid Y received
    void keepY(bool useW, vector<double> inset);

    //useful to change output file if you run this many times
    void changeFile(string out);

    //some variables related to the number of DM and the last Y list
    int n;
    vector<double> lastY;

protected:

    //output stream
    ofstream writer;

};
