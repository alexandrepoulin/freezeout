#include "../include/OutputWriter.h"
#include <iomanip>

//Constructor
OutputWriter::OutputWriter(string out,int nx){

    n=nx;
    lastY.resize(n);
    writer.open(out);

    if(!writer) cout << "Error opening output file " << out << endl;
    
}

//lets you start writting to another file
void OutputWriter::changeFile(string out){

    close();
    writer.open(out);
    if(!writer) cout << "Error opening output file " << out << endl;


}

//checks to see if the stream is open, and if it is, close it
void OutputWriter::close(){ if(writer.is_open()) writer.close();}

//writes all the relevent information to the output file
void OutputWriter::write(double temp, vector<double> eqset, vector<double> inset, bool useW){
    
    if(!writer.is_open()){
        cout << "Trying to write with closed writer." << endl;
        return;
    }
    writer << temp;
    if(useW) for(int i=0; i < n; i++) writer << "    " << exp(eqset[i]) << "    " << exp(inset[i]);
    else for(int i=0; i < n; i++) writer << "    " << eqset[i] << "    " << inset[i];
    writer << " " << endl;

}

void OutputWriter::writeOmega(vector<double> omega){
    if(!writer.is_open()){
        cout << "Trying to write with closed writer." << endl;
        return;
    }
    for(int i=0; i < n-1; i++) writer <<setprecision(7) << omega[i] << "    ";
    writer << " " << endl;
}

//Keeps tracks of the last valid Y received
void OutputWriter::keepY(bool useW,vector<double> inset){

    transform(inset.begin(),inset.end(),lastY.begin(), [=](double x)->double{return useW? exp(x):x;});

}
