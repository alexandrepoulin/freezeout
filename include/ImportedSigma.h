#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <set>
#include <math.h>

using namespace std;

//This object will handle imported cross sections. The formatting for the imported file will 
//have all the cross sections values on one line along with the temperature as the first entry.
//They will be ordered by (ijkl) where (ijkl) denotes the process i+j->k+l. Not all values of
//i,j,k,l are allowed. Here are the conditions for a valid ijkl given Nx=n:
//  (i,j) and (k,l) are ordered pairs (i<j and k<l)
//  (i,j) != (k,l), or in other words i!=j, k!=l
//The ordering can be given by a few ways. Lets compare (ijkl) and (i',j',k',l') and the 
//condition where (ijkl) < (i'j'k'l'):
//  i<i' or if i==i', j<j' or if j==j', k<k' or if k==k', l<l'
//  i*n^3 + j*n^2 + k*n + l < i'*n^3 + j'*n^2 + k'*n + l'
//  Index(ijkl) < Index(i'j'k'l')
//The Index function will be explaind later once more things are defined. The number of valid
//(ijkl) is given by NTot=T_n*(T_(n+1)-1) where T_n is the nth triangle number given by n(n+1)/2.
//In terms of n, T_n*T_(n+2)/2 = n(n+1)(n+2)(n+3)/8.
//The Index function which takes (ijkl)->d for some 0<=d<Ntot, is found by counting how many
//valid indices are passed for each increase in the digit. It is therefore the sum of 4 terms:
//iTerm=1/8*(4+2n-i)*(i-1)*(6+2n*(4+n)-5i-2ni+i^2)
//jTerm=1/2(j-i)(n^2+i^2-2ni+5n-4i-j+5)
//kTerm= (k==i)? 0 : n-j+1+1/2*(k-1-i)*(2n+2+i-k)
//lTerm= (k==i)? l-j-1: l-k
//d=iTerm+jTerm+kTerm+lTerm
//
//To get ijkl from d, we first not that we always have Delta(iTerm)>Delta(jTerm)>Delta(kTerm)>Delta(lTerm)
//That is to say, changing i by 1 will generate a bigger difference in the iTerm, then changing
//j by 1 in jTerm, and so on. We thus solve iTerm==d for i, then floor. We then solve 
//jTerm==d-iTerm for j and floor, and so on.
//
//This function will give values from 0 to NTot-1. Now make sure the cross section file follows 
//these rules or nothing will make sense.
//
//The temperature will be assumed to have the range of x from 1 to 1000 starting from x=1, i.e.
//The nth entry in the temperature vector will be given by T[n]=T[0]*L/(L+n*999) where L is the
//length of the vector and T[0] is the massScale
//
//The VERY FIRST LINE of the file should be a number representing the number of DM species.

class ImportedSigma
{

public:
    //Constructor
    ImportedSigma(){};
    ImportedSigma(string file, int N, double m);
    //Destructor
    ~ImportedSigma(){};

    //opens the file and sets some variables
    void initialize(string file, int N, double m);

    //functions to get the thermally averaged cross section
    double getSigmav(int Tind, int ind,double T);
    double getSigmav(int Tind, int i, int j, int k, int l,double T);
    double getSigmav(double T, int ind);
    double getSigmav(double T, int i, int j, int k, int l);
    //functions to get the temperature from the index associated with that temperature
    double getTemp(int Tind);
    double getTemp(double T);

    //gives you a vector filled with i, j, k, l
    vector<int> ijklFromInd(int d);
    //gives you the index given i, j, k, l
    int indFromijkl(int i, int j, int k, int l);
    
    //the total number of entries in sigma
    int tot; 
    
    //mass scale
    double massScale;

    //This is for the starting x
    vector<double> smallestSigv;
    set<int> nonZeroIndicesNW; //indices which have a least one non-zero cross section


protected:

    //some index and helper functions
    int Tindex(double t);
    bool compare(int i, int j, int k, int l); //returns (ij)<(kl)
    void read(string file);
    
    //some variables for the combinatorics involved in the indexing
    int n;
    int len; // temp length
    double step; //step in x for the defined cross sections


    //input stream, cross sections and temperature
    ifstream input;
    vector< vector<double> > sigmav;
    vector<double> Temp;

};
