#include "../include/NumericalCalculator.h"

int main(){
    
    string in = "C:/Users/AlexandrePoulin/Documents/relic/freezeout_v2.4/inputFiles/inputData.txt";
    string out = "C:/Users/AlexandrePoulin/Documents/relic/freezeout_v2.4/outputFiles/outputData.txt";
    string outOm = "C:/Users/AlexandrePoulin/Documents/relic/freezeout_v2.4/outputFiles/outputOmegaPDatag1.txt";
    
    
    cout <<"starting" << endl;
    NumericalCalculator numCalc(in,out);
    cout <<"numCalc allocated" << endl;

    bool test = true;
    if(test){
        cout << numCalc.inputReader->Nx << endl;
        numCalc.solve(true);
        cout << numCalc.omega[0]<<endl;
        cout << numCalc.omega[1]<<endl;
        cout << numCalc.omega[2]<<endl;
        cout << numCalc.omega[3]<<endl;
        cout << numCalc.omega[4]<<endl;
    }
    else{
        //double maxXsection = 5.5*pow(10,-26)*numCalc.CMCBSMTOPB;
        //double minXsection = 1.5*pow(10,-26)*numCalc.CMCBSMTOPB;
        double maxXsection = 100.0;
        double minXsection = 10.0;
        cout << "max = " << maxXsection << endl;
        cout << "min = " << minXsection << endl;

        double omega = 0.1120;//-0.0056;
        double omegaErr = 0.0001;

        ofstream writer;
        writer.open(outOm);
        if(!writer){
            cout << "Error opening omega writer" << endl;
            return -1;
        }

        int part = 50;

        double currentMass = 0;
        for(int i = 0; i < part+1; i++ ){
            cout <<"On i = " << i << " of " << part << endl;
            currentMass = pow(10,i/(double)part*5.0-1);
            //currentMass = 4.3+pow(10,i/(double)part-2.0);
            numCalc.inputReader->annPref/=numCalc.inputReader->massScale;
            numCalc.inputReader->mxset[0] = currentMass;
            numCalc.inputReader->massScale = currentMass;
            numCalc.inputReader->annPref*=numCalc.inputReader->massScale;
            
            double currentMax = maxXsection;
            double currentMin = minXsection;
            double currentAve = 0.0;
            double currentOmega = 0.0;
            while(abs(currentOmega-omega)>omegaErr){    
                currentAve = (currentMax+currentMin)*0.5;
                numCalc.inputReader->setSigvijklp(0,0,1,1,currentAve);
                numCalc.solve(false);
                currentOmega = numCalc.omega[0];

                //cout << currentMass << "    " << currentAve/numCalc.CMCBSMTOPB << "    " << currentOmega << endl;
                
                if(currentOmega < omega){
                    currentMax = currentAve;
                }
                else{
                    currentMin = currentAve;
                }
            }
            writer << currentMass << "    " << currentAve/numCalc.CMCBSMTOPB  << endl;

        }
        writer.close();
    }
    return 1;

}
