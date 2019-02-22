#include "../include/YGReader.h"

//Constructor
YGReader::YGReader(){
}

//functions to access the data
double YGReader::getYeq(double x){ return YeqPreCalc[YeqPreCalcInd(x)];}
double YGReader::getgstar(double T){ return gstarData[gindex(T)];}
double YGReader::getgstarS(double T){ return gstarSData[gindex(T)];}
//double YGReader::getgstarSDer(double T){ return gstarSDerData[gindex(T)];}
double YGReader::getgstarSDer(double T){ return 0.;}



//This function reads in the data
void YGReader::read(bool impY){
    
    cout << "Initializing YGReader" <<endl; 
    
    if(impY){
        //read the precalculated values of Y which take into account the changing gstar
        inputY.open("C:/Users/AlexandrePoulin/Documents/relic/freezeout_v2.5/YeqVals2.txt");
        if(!inputY) cout << "Error opening precalculated Yeq values. Defaulting to calculating them.";
        else{
            istream_iterator<double> start(inputY), end;
            YeqPreCalc = vector<double>(start,end);
        }
    }

    //len of Yeq vector
    lenPreCalc = YeqPreCalc.size(); 
 
    //Get the data for gstar and gstars and gstarsDer
    inputg.open("C:/Users/AlexandrePoulin/Documents/relic/freezeout_v2.5/gstarDataBig.txt");
    if(!inputg) cout << "Error opening gstar values.";
    else{
        istream_iterator<double> start(inputg), end;
        gstarData = vector<double>(start,end);
    }
    inputgs.open("C:/Users/AlexandrePoulin/Documents/relic/freezeout_v2.5/gstarSDataBig.txt");
    if(!inputgs) cout << "Error opening gstarS values.";
    else{
        istream_iterator<double> start(inputgs), end;
        gstarSData = vector<double>(start,end);
    }
    inputgsDer.open("C:/Users/AlexandrePoulin/Documents/relic/freezeout_v2.5/gstarSDerDataBig.txt");
    if(!inputgsDer) cout << "Error opening gstarS values.";
    else{
        istream_iterator<double> start(inputgsDer), end;
        gstarSDerData = vector<double>(start,end);
    }

    //close the file
    inputY.close();

}

//This function gives the index for the precalculated Y values
int YGReader::YeqPreCalcInd(double x){

    if(x<0.1){return 0;} //really high T, should not be in this case
    if(x>752.117) return lenPreCalc-1; //really low T, again, should not be in this case
    return ceil(25000000.0-2500000.0/x);

}

//This function gives the index for the precalculated gstar and gstarS values
int YGReader::gindex(double T)
{

    int val = round(10000.0/9*(log10(T)+5));
    if(val<0) return 0;
    if(val>=10001) return 10000;
    return val;

}


