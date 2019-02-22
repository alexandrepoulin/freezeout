#include "../include/ImportedSigma.h"

//Constructor
ImportedSigma::ImportedSigma(string file, int N, double m){
    initialize(file,N,m);    
}

//Initializes the instance

void ImportedSigma::initialize(string file, int N, double m){

    cout << "Initializing ImportedSigma" << endl;
    massScale = m;
    n = N;
    tot=n*(n+1)*(n+2)*(n+3)/8;
    smallestSigv.resize(tot);
    //for(int i = 0; i<tot; i++) smallestSigv[i]=100000;
    //DEBUG
    for(int i = 0; i<tot; i++) smallestSigv[i]=0;
   
    read(file);
}

//functions to get the cross section
//This one does a linear extrapolation between the two points
double ImportedSigma::getSigmav(int Tind, int ind,double T){ 
    
    //This does a linear regression, and its not good enough
/*
    if(Tind==len-1) return sigmav[Tind][ind];
    double val1 = sigmav[Tind][ind];
    double val2 = sigmav[Tind+1][ind];
    double slope = (val1-val2)/ (Temp[Tind]-Temp[Tind+1]);
    double b = val1-slope*Temp[Tind];
    return (slope*T+b)>0?slope*T+b:0;
*/

    //we are going to fit a power law to two points.
    //First, let's make sure we are not at the end

    if(Tind>=len-1) return sigmav[Tind][ind];

    if(sigmav[Tind][ind]==0 ||  sigmav[Tind+1][ind]==0) return 0;

    //Next we shift the two rightmost point by the leftmost point
    //to put it at the origin

    double x1=log(massScale/Temp[Tind]);
    double x2=log(massScale/Temp[Tind+1]);
    double y1=log(sigmav[Tind][ind]);
    double y2=log(sigmav[Tind+1][ind]);

    //Do the regression
    double slope = (y1-y2)/(x1-x2);
    double b = y1-x1*slope;

    return exp(log(massScale/T)*slope+b);

}

//These all call the one that does the linear extrapolation
double ImportedSigma::getSigmav(int Tind, int i, int j, int k, int l, double T){ return getSigmav(Tind,indFromijkl(i,j,k,l),T);}
double ImportedSigma::getSigmav(double T, int ind){ return getSigmav(Tindex(T),ind,T);}
double ImportedSigma::getSigmav(double T, int i, int j, int k, int l){ return getSigmav(Tindex(T),indFromijkl(i,j,k,l),T);}

//Functions to access the temperature
double ImportedSigma::getTemp(int Tind){ return Temp[Tind];}
double ImportedSigma::getTemp(double T){ return getTemp(Tindex(T));}

//Some useful functions
vector<int> ImportedSigma::ijklFromInd(int d){
    
    //This is found from inverting the relationships from indFromijkl sequencially since they get consecutively smaller
    int i =floor(0.5*(3+2*n-sqrt(5+4*sqrt(-8*d+pow(1+3*n+n*n,2)))));
    int j =floor(0.5*(1-i+i*i+3*n-2*i*n+n*n-sqrt(1-8*d+6*n+11*n*n+6*pow(n,3)+pow(n,4))));
    int k =floor(0.5*(3+2*n-sqrt(1-8*d-pow(i,4)-4*j*j+4*n*(n+3)+pow(i,3)*(2+4*n)+4*j*(-1+n*(n+3))+2*i*(1+2*n)*(-1-2*j+n*(n+3))+i*i*(1+4*j-2*n*(5+3*n)))));

    int iTerm = (3+2*n-i)*i*(2+2*n*(3+n)-3*i-2*i*n+i*i)/8;
    int jTerm = (j-i)*(n*n+i*i-2*n*i+3*n-2*i-j+1)/2;
    int kTerm = k > i ? n-j+(k-1-i)*(2*n+2-i-k)/2 : 0;
    int tempD=d-iTerm-jTerm-kTerm; 
    int l = ( k == i ? 1+tempD+j :tempD+k);

    
    vector<int> ijkl ={i,j,k,l}; //c++98 warning, ignore
    return ijkl;
}


//get the index from the ijkl values starting from 0
//The i, j, k, l should start from 0, not 1
int ImportedSigma::indFromijkl(int i, int j, int k, int l){
    int iTerm = (3+2*n-i)*i*(2+2*n*(3+n)-3*i-2*i*n+i*i)/8;
    int jTerm = (j-i)*(n*n+i*i-2*n*i+3*n-2*i-j+1)/2;
    int kTerm = k > i ? n-j+(k-1-i)*(2*n+2-i-k)/2 : 0;
    int lTerm = k == i ? l-j-1:l-k;
    return iTerm+jTerm+kTerm+lTerm;
}

//index function for temperature
int ImportedSigma::Tindex(double t){
    if(t>Temp[0]) return 0; 
    int ind = floor((massScale/t-massScale/Temp[0])/step);
    if(ind>=len){
        cout << "Cross Section called with too small of temperature." <<  endl;
        return len-1;
    }
    return ind;

}

//returns (ij)<(kl)
bool ImportedSigma::compare(int i, int j, int k, int l){ return i<k || (i==k && j<l);}

//actually read in the cross section
void ImportedSigma::read(string file){
    
    input.open(file);

    if(!input){
        cout << "Error opening Cross Section file." << endl;
        return;
    }

    double temp;
    do{
        input >> temp;
        Temp.push_back(temp); // Record the temperature
        vector<double> tempVec; //list of cross sections
        for(int i=0; i<tot; i++){
            input >> temp; //grab a cross section
            if(smallestSigv[i]<temp) smallestSigv[i]=temp;

            if(temp!= 0. && nonZeroIndicesNW.find(i)==nonZeroIndicesNW.end()) nonZeroIndicesNW.insert(i);
            
            tempVec.push_back(temp);
        }
        sigmav.push_back(tempVec);
        input.ignore(256,'\n'); //go to next line
    }while(!input.eof());
   

    len=Temp.size();
    step = 999.0/(len-1);

}
