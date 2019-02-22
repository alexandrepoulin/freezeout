#include "../include/InputReader.h"

//constructor
InputReader::InputReader(string in, int nDM){

    filename = in;
    read(nDM);
}

//allows to change the file
void InputReader::changeFile(string in, int nDM){
    
    filename = in;
    read(nDM);

}

//Sets both sigvijkl and sigvijklp to have a bunch of zeros
void InputReader::setZeroSigmav(){
    sigvijkl.clear();
    sigvijklp.clear();
    nonZeroIndices.clear();
    sigvijkl.resize(Nf4);
    sigvijklp.resize(Nf4);
    for(int i = 0; i <Nf4; i++){
        sigvijkl[i]=0;
        sigvijklp[i]=0;
    }

}

//changes only the i j -> k l sigvijkl 
void InputReader::setSigvijkl(int i, int j, int k, int l, double val){
    
    sigvijkl[Nf3*i+Nf2*j+Nf*k+l] = val;
    if(val!=0) nonZeroIndices.insert(Nf3*i+Nf2*j+Nf*k+l);
    else nonZeroIndices.erase(Nf3*i+Nf2*j+Nf*k+l); 

}

//changes only the i j -> k l sigvijklp 
void InputReader::setSigvijklp(int i, int j, int k, int l, double val){
    
    sigvijklp[Nf3*i+Nf2*j+Nf*k+l] = val;
    if(val!=0) nonZeroIndices.insert(Nf3*i+Nf2*j+Nf*k+l);
    else nonZeroIndices.erase(Nf3*i+Nf2*j+Nf*k+l); 

}

//Get the cross section (regardless of method of cross section)
double InputReader::getSigmav(double t,int ind){
    
    if(waves){
        return sigvijkl[ind]+t/massScale*sigvijklp[ind];
    }
    return inSig.getSigmav(t,ind);

}
double InputReader::getSigmav(double t,int i,int j, int k, int l){
    
    return inSig.getSigmav(t,i,j,k,l);

}

//Get the decay width 
double InputReader::getGamma(int ind){
    
    return gammaijk[ind];

}
double InputReader::getGamma(int i,int k, int l){
    
    int ind = Nf2 * i + Nf * k + l; 
    return getGamma(ind);

}

void InputReader::changeMassScale(double m){
    massScale = m;
    annPref = PBTOGEVM2 * sqrt(PI/45) * MP * massScale;
    inSig.massScale = m;
    startingX();
    setupIntervals();
    findNSteps();
}

//give the value of the step size given a bvalue of t
double InputReader::step(double t){

    double x = massScale/t;
    for(int i =0; i< intervals.size(); i++) if(x> intervals[i]->at(0) && x< intervals[i]->at(1)) return step1;
    return step2;

}

//change the smallest possible starting X

void InputReader::changeSmallestX(double x){
    smallestX = x;
    startingX();
    setupIntervals();
    findNSteps();
}

//Find how many steps will occur duging the numerical integration
void InputReader::findNSteps(){
    
    double smallLen = 0;

    for(int i=0; i< intervals.size(); i++) smallLen += intervals[i]->at(1)-intervals[i]->at(0);

    Nsteps = smallLen/step1 + (1000-intervals[0]->at(0)-smallLen)/step2;

}

//Setup the intervals for the smaller step size for use in the step function
void InputReader::setupIntervals(){
    
    intervals.clear();
    vector<double> copy(startingxVec.begin(),startingxVec.end());
    double minSize = 25;

    if(copy.size()==1){
        intervals.push_back(new vector<double>(2,copy[0]));
        intervals[0]->at(1)+=minSize;
    }
    else{
        sort(copy.begin(),copy.end());
        
        vector<double>* tempVec = new vector<double>(2, copy[0]);
        tempVec->at(1)+=minSize;
        for(int i = 1; i < copy.size(); i++){
            if(copy[i]<= tempVec->at(1)) tempVec->at(1)=copy[i]+minSize;
            else{
                intervals.push_back(tempVec);
                tempVec = new vector<double>(2,copy[i]);
                tempVec->at(1)+=minSize;
            }
        }
    }

}

//Gives the starting value of x where the distinction is made between Yeq and Y
void InputReader::startingX() {

	//for the waves only
	set<int>::iterator it;

	//setup some arrays
	double sigs[Nf];
	double sigp[Nf];
	for (int i = 0; i < Nf; i++) {
		sigs[i] = 0;
		sigp[i] = 0;
	}
	//If we only have p wave, then the formula is slightly different
	bool nonZeroSwave = false;

	if(waves) {

		//The case for waves

		//looping over the non-zero cross section
		for (it = nonZeroIndices.begin(); it != nonZeroIndices.end(); ++it) {

			//getting the indices
			//using unsigned int because it makes division slightly faster;
			unsigned int l = (*it) % Nf;
			unsigned int k = (((*it) % Nf2) - l) / Nf;
			unsigned int j = (((*it) % Nf3) - ((*it) % Nf2)) / Nf2;
			unsigned int i = ((*it) - ((*it) % Nf3)) / Nf3;

			//calculating the values for s wave, converting it to GeV^-2
			double vals = sigvijkl[*it] * PBTOGEVM2;

			//adding the contribution to the arrays for the s-wave
			if (i != Nx) sigs[i] += vals;
			if (j != Nx && i != j) sigs[j] += vals;
			if (vals != 0) nonZeroSwave = true;

			//calculating the values for p wave, converting it to GeV^-2
			double valp = sigvijklp[*it] * PBTOGEVM2;

			//adding the contribution to the arrays for the p-wave
			if (i != Nx) sigp[i] += valp;
			if (j != Nx && i != j) sigp[j] += valp;

		}
	}
	else {

		//The case for no waves

		//looping over the indices
		for (int p = 0; p < inSig.tot; p++) {
			//getting the indices
			vector<int> ijkl = inSig.ijklFromInd(p);
			int l = ijkl[3];
			int k = ijkl[2];
			int j = ijkl[1];
			int i = ijkl[0];

			//using the smallest values from each channel and assuming s wave for simplicity
			double val = inSig.smallestSigv[p] * PBTOGEVM2;

			//DEBUG
			//cout << inSig.smallestSigv[p] << " " << i << " " << j << " " << k << " " <<l << endl;

			if (i != Nx) sigs[i] += val;
			if (j != Nx && i != j) sigs[j] += val;
			if (k != Nx) sigs[k] += val;
			if (l != Nx && l != k) sigs[l] += val;
			nonZeroSwave = true;
		}
	}
    
    //find the smallest x we need to for each particle except SM
	vector<double> answer;
	double val = 0;
	if (nonZeroSwave) {
		for (int i = 0; i < Nf - 1; i++) {
			//loop over the channels and calculate x (has s-wave)
			answer.push_back(
                    max(smallestX,
                    massScale/mxset[i]*0.4*(1+0.1*log(massScale/mxset[i]))*log((2 + c)*c*0.145*0.264 * 10 * MP * mxset[i]*sigs[i]) - 0.5*log(log((2 + c)*c*0.145*0.264 * 10 * MP* mxset[i]*sigs[i])) + log(1 + sigp[i] / sigs[i] / log(0.038 * 10 * MP* mxset[i] *sigs[i]))
                    )
                );
		}
	}
	else {
		//loop over the channels and calculate x (has no s-wave)
		for (int i = 0; i < Nf - 1; i++) {
			answer.push_back(
                    max(smallestX,massScale/mxset[i]*0.4*(1+0.1*log(massScale/mxset[i]))*log((2 + c)*c*0.145*0.264 * 10 * MP * mxset[i]*sigp[i]) - 0.5*log(log((2 + c)*c*0.145*0.264 * 10 * MP* mxset[i]*sigp[i]))
                        )
                    );
		}

	}
	//this will also include the SM case and since its bigger than 1000, it will always be the eqset
	answer.push_back(10000);
	startingxVec= vector<double>(answer.begin(),answer.end());

}

//sets all the variables from the input file
void InputReader::read(int nDM){

    //open the input stream
    input.open(filename);

    //Check to see if the file was successfully opened
    if(!input){
        cout << "Error Opening Input File" << endl;
        return;
    }
    
    //ignore first line in the file which is explanation
    input.ignore(256,'\n');
    input >> waves;

    input.ignore(256,'\n');
    input.ignore(256,'\n');
    input.ignore(256,'\n');
    
    //get the number of dark matter fields
    input >> Nx;
    if(nDM!=-1) Nx=nDM;
    Nf = Nx + 1; //total number of fields (including one SM)
    Nf2 = Nf * Nf;
    Nf3 = Nf2 * Nf;
    Nf4 = Nf2 * Nf2; // should be slightly faster than Nf3 * Nf

    
    //ignore the next few lines of explanation
    input.ignore(256,'\n');
    input.ignore(256,'\n');
    input.ignore(256,'\n');
    
    //get the masses of the dark matter fields
    mxset.resize(Nx+1);
    double maxMass = 0;
    for(int i = 0; i < Nx; i++){
        input >> mxset[i];
        if(mxset[i]>maxMass) maxMass = mxset[i];
    }
    mxset[Nx]=0;//for the standard model side
    massScale = maxMass;
   
    //ignore the next few lines of explanation
    input.ignore(256,'\n');
    input.ignore(256,'\n');
    input.ignore(256,'\n');
    
    //get the difference in number density of the dark matter fields
    nDiffSet.resize(Nx+1);
    for(int i = 0; i < Nx; i++) input >> nDiffSet[i];
    nDiffSet[Nx]=0;//for the standard model side

    //ignore the next few lines of explanation
    input.ignore(256,'\n');
    input.ignore(256,'\n');
    input.ignore(256,'\n');
    
    if(waves){
        input.ignore(256,'\n');
        input.ignore(256,'\n');
        input.ignore(256,'\n');
        
        
        //fill in the termally averaged cross sections 
        //termally averaged cross sections, the first two indices are incoming and last two are outgoing
        sigvijkl.resize( Nf4 ); //its faster to make one long array then to make a four dimensional array due to how fetching the memory works
        //keep track of which entries are non-zero
        for(int i = 0; i < Nf; i++){ //incoming particle
            for(int k = 0; k < Nf; k++){ // outgoing particle
                for(int j = 0; j < Nf; j++){ // incoming particle
                    for(int l = 0; l < Nf; l++){ //outgoing particle
                        input >> sigvijkl[(i * Nf3) + (j * Nf2) +  (k * Nf) + l];
                        if( sigvijkl[(i * Nf3) + (j * Nf2) +  (k * Nf) + l]!= 0){
                            nonZeroIndices.insert((i * Nf3) + (j * Nf2) +  (k * Nf) + l);
                        }
                    }
                }
                input.ignore(256,'\n'); //this is to go to the next line
            }    
        }    
        

        input.ignore(256,'\n');
        input.ignore(256,'\n');
        input.ignore(256,'\n');
        sigvijklp.resize( Nf4 ); //its faster to make one long array then to make a four dimensional array due to how fetching the memory works
        //keep track of which entries are non-zero
        for(int i = 0; i < Nf; i++){ //incoming particle
            for(int k = 0; k < Nf; k++){ // outgoing particle
                for(int j = 0; j < Nf; j++){ // incoming particle
                    for(int l = 0; l < Nf; l++){ //outgoing particle
                        input >> sigvijklp[(i * Nf3) + (j * Nf2) +  (k * Nf) + l];
                        if( sigvijklp[(i * Nf3) + (j * Nf2) +  (k * Nf) + l]!= 0){
                            nonZeroIndices.insert((i * Nf3) + (j * Nf2) +  (k * Nf) + l);
                        }
                    }
                }
                input.ignore(256,'\n'); //this is to go to the next line
            }    
        }    

        //skip some description
        input.ignore(256,'\n');
        input.ignore(256,'\n');
        input.ignore(256,'\n');
        
        gammaijk.resize( Nf3 ); //its faster to make one long array then to make a four dimensional array due to how fetching the memory works
        //keep track of which entries are non-zero
        for(int i = 0; i < Nf; i++){ // particle decaying
            for(int k = 0; k < Nf; k++){ // outgoing particle 1
                for(int j = 0; j < Nf; j++){ // outgoing particle 2
                    input >> gammaijk[(i * Nf2) + (j * Nf) +  k];
                    if( gammaijk[(i * Nf2) + (j * Nf) +  k]!= 0){
                        nonZeroIndicesGam.insert((i * Nf2) + (j * Nf) +  k);
                    }
                }
                input.ignore(256,'\n'); //this is to go to the next line
            }    
        }
        
        input.ignore(256,'\n');
        input.ignore(256,'\n');
    }
    else{
        //import the Cross section list
        getline(input,filenameInSig);
        filenameInSig.pop_back(); //string has a function called pop_back(), don't believe vims lies! This removes the '\n'
        inSig.initialize(filenameInSig,Nx, massScale);
        input.ignore(256,'\n');
        input.ignore(256,'\n');
    }

    //read the desired step size from the input file
    input >> step1;
    input >> step2;
 
    //skip some description
    input.ignore(256,'\n');
    input.ignore(256,'\n');
    input.ignore(256,'\n');
    
    input >> smallestX; //read the smallest possible starting value for x
    
    //skip some description
    input.ignore(256,'\n');
    input.ignore(256,'\n');
    input.ignore(256,'\n');
    
    input >> usePreCalc; //read in wether to use the precalculated Yeq
 
    //skip some description
    input.ignore(256,'\n');
    input.ignore(256,'\n');
    input.ignore(256,'\n');
    input >> method; //read in wether to use the precalculated

    //skip some description
    input.ignore(256,'\n');
    input.ignore(256,'\n');
    input.ignore(256,'\n');
    gset.resize(Nx+1);
    cout << "here" << endl; 
    for(int i = 0; i < Nx; i++) input >> gset[i];
    gset[Nx]=1; //For standard model values

    //skip some description
    input.ignore(256,'\n');
    input.ignore(256,'\n');
    input.ignore(256,'\n');
    input >> useW; //read in wether to use the precalculated


    //close the input File 
    input.close();
    
    annPref = PBTOGEVM2 * sqrt(PI/45) * MP * massScale;
    //if(!useRK4) annPref *= step;

    //Setup some of the things you need
    startingX();
    setupIntervals();
    findNSteps();

}
