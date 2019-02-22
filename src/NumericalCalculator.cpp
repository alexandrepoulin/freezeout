#include "../include/NumericalCalculator.h"

//Constructor
NumericalCalculator::NumericalCalculator(){
    cout << "Initializing InputReader" << endl;
    inputReader = new InputReader();
    cout << "Initializing outWriter" << endl;
    outputWriter = new OutputWriter();
    ygReader.read(inputReader->usePreCalc);
}
NumericalCalculator::NumericalCalculator(string in, string out, int nDM){
    cout << "Initializing InputReader" << endl;
    inputReader = new InputReader(in, nDM);
    cout << "Initializing outputWriter" << endl;
    outputWriter = new OutputWriter(out, inputReader->Nf);
    ygReader.read(inputReader->usePreCalc);
}

//Destructor
NumericalCalculator::~NumericalCalculator(){
    delete inputReader;
    outputWriter->close();
    delete outputWriter;
}

//Lets you change the input and output
void NumericalCalculator::newInput(string in){ inputReader->changeFile(in);}
void NumericalCalculator::newOutput(string out){ outputWriter->changeFile(out);

}

//get the value of mu as a function of temperature
double NumericalCalculator::getMu(double t, int ind){


    if(inputReader->nDiffSet[ind]==0.) return 0.;

    //for large x, the bessel function just gives 0 which is a problem.
    //we use the approximation that asinh(a*/K2(x)) approx asinh(b *exp(x)) = log(sqrt(b^2*exp(2x)+1)+b*exp(x)) approx +-(x+2log(abs(2b)))
    //the +- sign depends on the sign of b and b=sqrt(2x/pi)*a, and a is the rest of the stuff there
    if(t < inputReader->massScale / 500){
        return (inputReader->nDiffSet[ind]>0? 1:-1) * (inputReader->massScale + t * log(2*abs(inputReader->nDiffSet[ind]) * pow(inputReader->PI,2)/45 * ygReader.getgstarS(t)/inputReader->gset[ind]*pow(2*inputReader->PI*t/inputReader->mxset[ind],3./2)));
    }

    //exact value
      
    return t*asinh( inputReader->nDiffSet[ind] * 2 * pow(inputReader->PI,4) / 45 * ygReader.getgstarS(t) / inputReader->gset[ind] * pow(t/inputReader->massScale,2) / boost::math::cyl_bessel_k(2,inputReader->massScale/t) );

}

//Gives you the value of Yeq
double NumericalCalculator::Yeq(double m, double mu, double T, double gstarS, double g){


    if(T>inputReader->massScale/500){
        //Proper treatment
        return YeqPref * g/ gstarS * pow( m/T , 2) * boost::math::cyl_bessel_k(2, m/T) * exp(mu/T);
    }
    else{
        //use approximation
        return 45.0/(4*pow(3.14159265358979323,4)) * sqrt(3.14159265358979323/2) * g/ gstarS * pow( m/T , 3.0/2.0) * exp(-(m-mu)/T);
    }
}

//Generates the list of equilibrium values based on the temperature
vector<double> NumericalCalculator::genEqSet(double t){

    //temporary list
    vector<double> temp;
    temp.resize(inputReader->Nf);

    for(int i = 0; i < inputReader->Nx; i++){ 
        //make the new Yeq
        double val; 
        if(inputReader->usePreCalc){
            val = inputReader->gset[i]*ygReader.getYeq(inputReader->mxset[i]/t)* exp(getMu(t,i)/t);
        }
        else{
            val =Yeq(inputReader->mxset[i], getMu(t,i), t, ygReader.getgstarS(t), inputReader->gset[i]);
        }
        //If we use W, make sure that the large x behaviour works (i.e. no log(0) issues)
        if(inputReader->useW){
            temp[i] = isinf(log(val)) ? 3.0/2.0*log(inputReader->mxset[i]/t)-inputReader->mxset[i]/t+log(1.45 * ygReader.getgstar(t)/ ygReader.getgstarS(t)) 
            : log(val);
        }
        else{
            temp[i]=val;
        }
    }
    //set the radiation values
    if(inputReader->useW){
        temp[inputReader->Nx]=log(Yrad/ygReader.getgstar(t));
    }
    else{
        temp[inputReader->Nx]=Yrad/ygReader.getgstar(t);
    }

    return temp;
}

double NumericalCalculator::genEqRatio(int i, int j, int k, int l, double t, double inK, double inL){

    double prefactor = 1;
    double powerFactor= 1;
    double exponent = 0;
    if(i!= inputReader->Nx){
        prefactor *= YeqPref * inputReader->gset[i] / ygReader.getgstarS(t)*sqrt(3.14159265358979323/2);
        powerFactor *=inputReader->mxset[i];
        exponent -= inputReader->mxset[i]/t; 
    }
    else{
        prefactor *= Yrad/ygReader.getgstar(t);
    }
   
    if(j!= inputReader->Nx){
        prefactor *= YeqPref * inputReader->gset[j] / ygReader.getgstarS(t)*sqrt(3.14159265358979323/2);
        powerFactor *=inputReader->mxset[j];
        exponent -= inputReader->mxset[j]/t; 
    }
    else{
        prefactor *= Yrad/ygReader.getgstar(t);
    }
    
    if(k!= inputReader->Nx){
        prefactor /= YeqPref * inputReader->gset[k] / ygReader.getgstarS(t)*sqrt(3.14159265358979323/2);
        powerFactor /=inputReader->mxset[k];
        exponent += inputReader->mxset[k]/t; 
    }
    else{
        prefactor /= Yrad/ygReader.getgstar(t);
    }
   
    if(l!= inputReader->Nx){
        prefactor /= YeqPref * inputReader->gset[l] / ygReader.getgstarS(t)*sqrt(3.14159265358979323/2);
        powerFactor /=inputReader->mxset[l];
        exponent += inputReader->mxset[l]/t; 
    }
    else{
        prefactor /= Yrad/ygReader.getgstar(t);
    }

    double logRes = log(inK)+log(inL)+log(prefactor)+3./2.*log(powerFactor)+exponent;
    return exp(logRes);

}

//This is the f in Yi'=f(x,Y1,Y2,...,Yn). This is where the physics is
vector<double> NumericalCalculator::f(vector<double> inset, vector<double> eqset, double t){

    //DEBUG
    //cout <<"f called"<< endl;
    //temporary list
    vector<double> temp;
    temp.resize(inputReader->Nf); 
    
    //This is for the waves cases only
    set<int>::iterator it;

    //2 to 2 scattering
    if(inputReader->useW){
        if(inputReader->waves){

            //USES W AND WAVES
            
            //loop through the non-zero cross sections
            for(it=inputReader->nonZeroIndices.begin();it!=inputReader->nonZeroIndices.end();++it){
            
                //get the indices
                unsigned int l =(*it)%inputReader->Nf ; 
                unsigned int k =(((*it)%inputReader->Nf2) -  l)/inputReader->Nf; 
                unsigned int j =(((*it)%inputReader->Nf3) -  ((*it)%inputReader->Nf2))/inputReader->Nf2; 
                unsigned int i =((*it) -  ((*it)%inputReader->Nf3))/inputReader->Nf3; 
            
                //calculate the contribution
                double val = 0.5* ( 1.0 + ygReader.getgstarSDer(t)/3.0 ) *  inputReader->annPref * ygReader.getgstarS(t) / sqrt(ygReader.getgstar(t)) * pow(t,2) /pow(inputReader->massScale,2)* inputReader->getSigmav(t,*it)
                            * (exp(inset[i] + inset[j]) - exp(inset[k] + inset[l] + eqset[i] + eqset[j] - eqset[k] - eqset[l]));
                
                //if val < 0, then this means that the chemical potential of the final state is bigger than the initial state
                //so we would have the opposite process, in which case this should not contribute
                //if(val < 0) val =0;

                //putting in the values
                if(i!=inputReader->Nx) temp[i]-= exp(-inset[i])*val;
                if(j!=inputReader->Nx) temp[j]-= exp(-inset[j])*val;

                if(k!=inputReader->Nx) temp[k]+= exp(-inset[k])*val;
                if(l!=inputReader->Nx) temp[l]+= exp(-inset[l])*val;
            }
        }
        else{
            
            //USES W AND NO WAVES
            
            //loop through the non-zero cross sections
            for(it=inputReader->inSig.nonZeroIndicesNW.begin();it!=inputReader->inSig.nonZeroIndicesNW.end();++it){
                
                if(inputReader->getSigmav(t,*it)==0) continue;
                
                //get the indices
                vector<int> ijkl = inputReader->inSig.ijklFromInd((*it)); 
                int l = ijkl[3]; 
                int k = ijkl[2];
                int j = ijkl[1];
                int i = ijkl[0];
                
                //calculate the contribution
                double val =  ( 1.0 + ygReader.getgstarSDer(t)/3.0 ) *  inputReader->annPref * ygReader.getgstarS(t) / sqrt(ygReader.getgstar(t)) * pow(t,2) /pow(inputReader->massScale,2)* inputReader->getSigmav(t,(*it))
                            * (exp(inset[i] + inset[j]) - exp(inset[k] + inset[l] + eqset[i] + eqset[j] - eqset[k] - eqset[l]));

                //if val < 0, then this means that the chemical potential of the final state is bigger than the initial state
                //so we would have the opposite process, in which case this should not contribute
                //if(val < 0) val =0;
                
                //putting in the values
                if(i!=inputReader->Nx) temp[i]-= exp(-inset[i])*val;
                if(j!=inputReader->Nx) temp[j]-= exp(-inset[j])*val;

                if(k!=inputReader->Nx) temp[k]+= exp(-inset[k])*val;
                if(l!=inputReader->Nx) temp[l]+= exp(-inset[l])*val;
            }
        }
    }
    else{
        if(inputReader->waves){
            
            //USES Y AND WAVES
            
            //loop through the non-zero cross sections
            for(it=inputReader->nonZeroIndices.begin();it!=inputReader->nonZeroIndices.end();++it){
            
                //get the indices
                unsigned int l =(*it)%inputReader->Nf ; 
                unsigned int k =(((*it)%inputReader->Nf2) -  l)/inputReader->Nf; 
                unsigned int j =(((*it)%inputReader->Nf3) -  ((*it)%inputReader->Nf2))/inputReader->Nf2; 
                unsigned int i =((*it) -  ((*it)%inputReader->Nf3))/inputReader->Nf3; 
            
                double val = 0;
                
                
                /*
                //calculate the contribution (making sure that the large x behaviour does not give a 0/0)
                if(t> inputReader->massScale/600.0){
                    if(k!= inputReader->Nx && l !=inputReader->Nx){ 
                        val =  eqset[i] * (eqset[j] / eqset[k]) / eqset[l] == eqset[i] * (eqset[j] / eqset[k]) / eqset[l] ? 
                                 inset[i] * inset[j] - inset[k] * inset[l] * (eqset[i] * (eqset[j] / eqset[k]) / eqset[l]):
                                 inset[i] * inset[j] - inset[k] * inset[l] 
                                 * pow( inputReader->mxset[i] * inputReader->mxset[j] / (inputReader->mxset[k] * inputReader->mxset[l]), 3.0/2.0)  
                                 * exp(-((inputReader->mxset[i]+inputReader->mxset[j]-inputReader->mxset[k]-inputReader->mxset[l])/t));
                    }else if(k== inputReader->Nx && l!=inputReader->Nx){
                        val =  eqset[i] * (eqset[j] / eqset[k]) / eqset[l] == eqset[i] * (eqset[j] / eqset[k])/ eqset[l] ? 
                                 inset[i] * inset[j] - inset[k] * inset[l] * (eqset[i] * (eqset[j] / eqset[k]) / eqset[l]):
                                 inset[i] * inset[j] - inset[k] * inset[l] 
                                 * pow( inputReader->mxset[i] * inputReader->mxset[j] / (inputReader->mxset[l])/t, 3.0/2.0)  
                                 * exp(-((inputReader->mxset[i]+inputReader->mxset[j]-inputReader->mxset[l])/t))/eqset[k];
                    }else if(k!= inputReader->Nx && l==inputReader->Nx){
                        val =  eqset[i] * (eqset[j] / eqset[k]) / eqset[l] == eqset[i] * (eqset[j] / eqset[k]) / eqset[l] ? 
                                 inset[i] * inset[j] - inset[k] * inset[l] * (eqset[i] * (eqset[j] / eqset[k]) / eqset[l]):
                                 inset[i] * inset[j] - inset[k] * inset[l] 
                                 * pow( inputReader->mxset[i] * inputReader->mxset[j] / (inputReader->mxset[k])/t, 3.0/2.0)  
                                 * exp(-((inputReader->mxset[i]+inputReader->mxset[j]-inputReader->mxset[k])/t))/eqset[l];
                    }else{
                        val= inset[i] * inset[j] - inset[k] * inset[l] * (eqset[i] * (eqset[j] / eqset[k]) / eqset[l]);
                    
                    }
                }
                else val= inset[i] * inset[j] - inset[k] * inset[l] * genEqRatio(i,j,k,l,t);
                */
                
                if(t> inputReader->massScale/600){
                    double ratio = log(eqset[i])+log(eqset[j])-log(eqset[k])-log(eqset[l])+log(inset[k]) + log(inset[l]);
                    ratio=exp(ratio);
                    val= inset[i] * inset[j] -  ratio;
                }
                else{
                    val = inset[i] * inset[j] -genEqRatio(i,j,k,l,t,inset[k],inset[l]);
                }

                val *= ( 1.0 + ygReader.getgstarSDer(t)/3.0) *  inputReader->annPref * ygReader.getgstarS(t) / sqrt(ygReader.getgstar(t)) * pow(t,2) /pow(inputReader->massScale,2)* inputReader->getSigmav(t,*it);
                
                //if val < 0, then this means that the chemical potential of the final state is bigger than the initial state
                //so we would have the opposite process, in which case this should not contribute
                //if(val < 0) val =0;
               
         

                //putting in the values
                if(i!=inputReader->Nx) temp[i]-= val;
                if(j!=inputReader->Nx) temp[j]-= val;

                if(k!=inputReader->Nx) temp[k]+= val;
                if(l!=inputReader->Nx) temp[l]+= val;
            }
        }
        else{
            //USES Y AND NO WAVES
            
            //loop through the non-zero cross sections


            for(it=inputReader->inSig.nonZeroIndicesNW.begin();it!=inputReader->inSig.nonZeroIndicesNW.end();++it){
                
                if(inputReader->getSigmav(t,*it)==0) continue;
                //get the indices
                vector<int> ijkl = inputReader->inSig.ijklFromInd((*it)); 
                int l = ijkl[3]; 
                int k = ijkl[2];
                int j = ijkl[1];
                int i = ijkl[0];

                //DEBUG
                //cout<< p << " " << i << " " << j << " " << k << " " << l << " " << inputReader->getSigmav(t,p)<< endl;
                
                double val = 0;
                
                //calculate the contribution (making sure that the large x behaviour does not give a 0/0)
                
/*
                if(t>inputReader->massScale/600.0){
                    if(k!= inputReader->Nx && l !=inputReader->Nx){
                        val =  eqset[i] * (eqset[j] / eqset[k]) / eqset[l] == eqset[i] * (eqset[j] / eqset[k]) / eqset[l] ? 
                                 //inset[i] * inset[j] - inset[k] * inset[l] * (eqset[i] * (eqset[j] / eqset[k]) / eqset[l]):
                                 inset[i] * inset[j]*(1.0- (inset[k]/eqset[k]) * (inset[l]/eqset[l]) * (eqset[i]/inset[i]) * (eqset[j]/inset[j]) ):
                                 inset[i] * inset[j] - inset[k] * inset[l] 
                                 * pow( inputReader->mxset[i] * inputReader->mxset[j] / (inputReader->mxset[k] * inputReader->mxset[l]), 3.0/2.0)  
                                 * exp(-((inputReader->mxset[i]+inputReader->mxset[j]-inputReader->mxset[k]-inputReader->mxset[l])/t));
                    }else if(k== inputReader->Nx && l!=inputReader->Nx){
                        val =  eqset[i] * (eqset[j] / eqset[k]) / eqset[l] == eqset[i] * (eqset[j] / eqset[k]) / eqset[l] ? 
                                 inset[i] * inset[j] - inset[k] * inset[l] * (eqset[i] * (eqset[j] / eqset[k]) / eqset[l]):
                                 inset[i] * inset[j] - inset[k] * inset[l] 
                                 * pow( inputReader->mxset[i] * inputReader->mxset[j] / (inputReader->mxset[l])/t, 3.0/2.0)  
                                 * exp(-((inputReader->mxset[i]+inputReader->mxset[j]-inputReader->mxset[l])/t))/eqset[k];
                    }else if(k!= inputReader->Nx && l==inputReader->Nx){
                        val =  eqset[i] * (eqset[j] / eqset[k]) / eqset[l] == eqset[i] * (eqset[j] / eqset[k]) / eqset[l] ? 
                                 inset[i] * inset[j] - inset[k] * inset[l] * (eqset[i] * (eqset[j] / eqset[k]) / eqset[l]):
                                 inset[i] * inset[j] - inset[k] * inset[l] 
                                 * pow( inputReader->mxset[i] * inputReader->mxset[j] / (inputReader->mxset[k])/t, 3.0/2.0)  
                                 * exp(-((inputReader->mxset[i]+inputReader->mxset[j]-inputReader->mxset[k])/t))/eqset[l];
                    }else{
                        //val= inset[i] * inset[j] - inset[k] * inset[l] * (eqset[i] * (eqset[j] / eqset[k]) / eqset[l]);
                        val=inset[i] * inset[j]*(1.0- (inset[k]/eqset[k]) * (inset[l]/eqset[l]) * (eqset[i]/inset[i]) * (eqset[j]/inset[j]) );
                    
                    }
                }
                else val= inset[i] * inset[j] - inset[k] * inset[l] * genEqRatio(i,j,k,l,t);
                
*/
                if(t> inputReader->massScale/600){
                    double ratio = log(eqset[i])+log(eqset[j])-log(eqset[k])-log(eqset[l])+log(inset[k]) + log(inset[l]);
                    ratio=exp(ratio);
                    val= inset[i] * inset[j] -  ratio;
                }
                else{
                    val = inset[i] * inset[j] -genEqRatio(i,j,k,l,t,inset[k],inset[l]);
                }
                
                //DEBUG
                //cout << val<<" " << i << " " << j << " " <<k << " " << l << " " << inputReader->getSigmav(t, (*it))<<" " << 0.5 * ( 1.0 + ygReader.getgstarSDer(t)/3.0 ) *  inputReader->annPref * ygReader.getgstarS(t) / sqrt(ygReader.getgstar(t)) * pow(t,2) /pow(inputReader->massScale,2)<< endl;
                //val =  inset[i] * inset[j] - inset[k] * inset[l] * genEqRatio(i,j,k,l,t);

                //DEBUG
                //cout <<val <<endl;

                
                val *=  ( 1.0 + ygReader.getgstarSDer(t)/3.0 ) *  inputReader->annPref * ygReader.getgstarS(t) / sqrt(ygReader.getgstar(t)) * pow(t,2) /pow(inputReader->massScale,2)* inputReader->getSigmav(t,(*it));
                //DEBUG
                //cout <<val <<endl;

                //if val < 0, then this means that the chemical potential of the final state is bigger than the initial state
                //so we would have the opposite process, in which case this should not contribute
                //if(val < 0) val =0;

                //putting in the values
                if(i!=inputReader->Nx) temp[i]-= val;
                if(j!=inputReader->Nx) temp[j]-= val; 

                if(k!=inputReader->Nx) temp[k]+= val;
                if(l!=inputReader->Nx) temp[l]+= val;
            }
        }
    }

    //decays (only implemented for waves so far)
    if(inputReader->waves){
        if(inputReader->useW){

            //USES W AND WAVES
            
            //loop through the non-zero cross sections
            for(it=inputReader->nonZeroIndicesGam.begin();it!=inputReader->nonZeroIndicesGam.end();++it){
            
                //get the indices
                unsigned int k =(*it)%inputReader->Nf ; 
                unsigned int j =(((*it)%inputReader->Nf2) -  k)/inputReader->Nf; 
                unsigned int i =(((*it)) -  ((*it)%inputReader->Nf2))/inputReader->Nf2; 
            
                //calculate the contribution
                double val =  ( 1.0 + ygReader.getgstarSDer(t)/3.0 ) *  inputReader->annPref * ygReader.getgstarS(t) / sqrt(ygReader.getgstar(t)) * pow(t,2) /pow(inputReader->massScale,2)* inputReader->getGamma(*it)
                            * (exp(inset[i]) - exp(inset[j] + inset[k] + eqset[i] - eqset[j] - eqset[k]));
                
                //if val < 0, then this means that the chemical potential of the final state is bigger than the initial state
                //so we would have the opposite process, in which case this should not contribute
                //if(val < 0) val =0;

                //putting in the values
                if(i!=inputReader->Nx) temp[i]-= exp(-inset[i])*val;
                
                if(j!=inputReader->Nx) temp[j]+= exp(-inset[j])*val;
                if(k!=inputReader->Nx) temp[k]+= exp(-inset[k])*val;
            }
        }
        else{
            //USES Y AND WAVES
            
            //loop through the non-zero cross sections
            for(it=inputReader->nonZeroIndicesGam.begin();it!=inputReader->nonZeroIndicesGam.end();++it){
            
                //get the indices
                unsigned int k =(*it)%inputReader->Nf ; 
                unsigned int j =(((*it)%inputReader->Nf2) -  k)/inputReader->Nf; 
                unsigned int i =(((*it)) -  ((*it)%inputReader->Nf2))/inputReader->Nf2; 
            
                double val = 0;
                
                //calculate the contribution (making sure that the large x behaviour does not give a 0/0)
                if(k!= inputReader->Nx && j !=inputReader->Nx){
                    val =  eqset[i] / eqset[k] / eqset[j] == eqset[i] / eqset[k] / eqset[j] ? 
                             inset[i] - inset[k] * inset[j] * (eqset[i] / eqset[k]) / eqset[j]:
                             inset[i] - inset[k] * inset[j] 
                             * pow( t*inputReader->mxset[i] / (inputReader->mxset[k] * inputReader->mxset[j]), 3.0/2.0)  
                             * exp(-((inputReader->mxset[i]-inputReader->mxset[k]-inputReader->mxset[j])/t));
                }else if(k== inputReader->Nx && j!=inputReader->Nx){
                    val =  eqset[i] / eqset[k] / eqset[j] == eqset[i] / eqset[k] / eqset[j] ? 
                             inset[i] - inset[k] * inset[j] * (eqset[i] / eqset[k]) / eqset[j]:
                             inset[i] - inset[k] * inset[j] 
                             * pow( inputReader->mxset[i] / (inputReader->mxset[j]), 3.0/2.0)  
                             * exp(-((inputReader->mxset[i]-inputReader->mxset[j])/t))/eqset[k];
                }else if(k!= inputReader->Nx && j==inputReader->Nx){
                    val =  eqset[i] / eqset[k] / eqset[j] == eqset[i] / eqset[k] / eqset[j] ? 
                             inset[i] - inset[k] * inset[j] * (eqset[i] / eqset[k]) / eqset[j]:
                             inset[i] - inset[k] * inset[j] 
                             * pow( inputReader->mxset[i] / (inputReader->mxset[k]), 3.0/2.0)  
                             * exp(-((inputReader->mxset[i]-inputReader->mxset[k])/t))/eqset[j];
                }else{
                    val= inset[i] - inset[k] * inset[j] * (eqset[i] / eqset[k]) / eqset[j];
                
                }
                
                val *= ( 1.0 + ygReader.getgstarSDer(t)/3.0) *  inputReader->annPref * ygReader.getgstarS(t) / sqrt(ygReader.getgstar(t)) * pow(t,2) /pow(inputReader->massScale,2)* inputReader->getGamma(*it);
                
                //if val < 0, then this means that the chemical potential of the final state is bigger than the initial state
                //so we would have the opposite process, in which case this should not contribute
                //if(val < 0) val =0;
                

                //putting in the values
                if(i!=inputReader->Nx) temp[i]-= val;
                
                if(j!=inputReader->Nx) temp[j]+= val;
                if(k!=inputReader->Nx) temp[k]+= val;
            }

        }
    }
    
    
    //return the contributions
    return temp;
}


//Gives the starting value of x where the distinction is made between Yeq and Y
/*vector<double> NumericalCalculator::startingX(){
    
    //for the waves only
    set<int>::iterator it;
    
    //setup some arrays
    double sigs[inputReader->Nf]; 
    double sigp[inputReader->Nf];
    for(int i = 0; i < inputReader->Nf; i++){
        sigs[i]=0;
        sigp[i]=0;
    }
    //If we only have p wave, then the formula is slightly different
    bool nonZeroSwave= false;

    if(inputReader->waves){

        //The case for waves
        
        //looping over the non-zero cross section
        for(it=inputReader->nonZeroIndices.begin();it!=inputReader->nonZeroIndices.end();++it){
            
            //getting the indices
            //using unsigned int because it makes division slightly faster;
            unsigned int l =(*it)%inputReader->Nf ; 
            unsigned int k =(((*it)%inputReader->Nf2) -  l)/inputReader->Nf; 
            unsigned int j =(((*it)%inputReader->Nf3) -  ((*it)%inputReader->Nf2))/inputReader->Nf2; 
            unsigned int i =((*it) -  ((*it)%inputReader->Nf3))/inputReader->Nf3; 

            //calculating the values for s wave, converting it to GeV^-2
            double vals =  inputReader->sigvijkl[*it]* inputReader->PBTOGEVM2;
            
            //adding the contribution to the arrays for the s-wave
            if(i!=inputReader->Nx) sigs[i]+= vals;
            if(j!=inputReader->Nx && i!=j) sigs[j]+= vals;
            if(vals !=0) nonZeroSwave=true;
            
            //calculating the values for p wave, converting it to GeV^-2
            double valp = inputReader->sigvijklp[*it]*inputReader->PBTOGEVM2;
            
            //adding the contribution to the arrays for the p-wave
            if(i!=inputReader->Nx) sigp[i]+= valp;
            if(j!=inputReader->Nx && i!=j) sigp[j]+= valp;

        }
    }
    else{
        
        //The case for no waves
        
        //looping over the indices
        for(int p = 0; p < inputReader->inSig.tot; p++){
            //getting the indices
            vector<int> ijkl = inputReader->inSig.ijklFromInd(p); 
            int l = ijkl[3]; 
            int k = ijkl[2];
            int j = ijkl[1];
            int i = ijkl[0];

            //using the smallest values from each channel and assuming s wave for simplicity
            double val  = inputReader->inSig.smallestSigv[p]* inputReader->PBTOGEVM2;
           
            //DEBUG
            //cout << inputReader->inSig.smallestSigv[p] << " " << i << " " << j << " " << k << " " <<l << endl;

            if(i!=inputReader->Nx) sigs[i]+= val;
            if(j!=inputReader->Nx && i!=j) sigs[j]+= val;
            nonZeroSwave= true;
        }
    }
/ *
    //We want to start at the smallest x we need to, so we keep track of which channel will need the smallest x
    double minSoFar=pow(10,9);
    double val =0;
    if(nonZeroSwave){
        for(int i = 0; i < inputReader->Nf-1; i++){
            //loop over the channels and calculate x (has s-wave)
            val=log((2+c)*c*0.145*0.264*10* inputReader->MP * inputReader->massScale*sigs[i])-0.5*log(log((2+c)*c*0.145*0.264*10* inputReader->MP* inputReader->massScale*sigs[i]))+log(1+sigp[i]/sigs[i]/log(0.038*10* inputReader->MP* inputReader->massScale *sigs[i]));
            if(val<minSoFar) minSoFar = val;
        }
    }
    else{
        //loop over the channels and calculate x (has no s-wave)
        for(int i = 0; i < inputReader->Nf-1; i++){
            val=log((2+c)*c*0.145*0.264*10* inputReader->MP * inputReader->massScale*sigp[i])-0.5*log(log((2+c)*c*0.145*0.264*10* inputReader->MP* inputReader->massScale*sigp[i]));
            if(val<minSoFar) minSoFar = val;
        }
    
    }
    //Make sure we aren't starting before 1
    if(minSoFar < 1) minSoFar = 1/0.9;
    return minSoFar*0.9; // multiply by 0.9 to be super duper sure that we aren't missing stuff

* /

   //find the smallest x we need to for each particle except SM
    vector<double> answer;
    double val =0;
    if(nonZeroSwave){
        for(int i = 0; i < inputReader->Nf-1; i++){
            //loop over the channels and calculate x (has s-wave)
            answer.push_back(2.0*log((2+c)*c*0.145*0.264*10* inputReader->MP * inputReader->massScale*sigs[i])-0.5*log(log((2+c)*c*0.145*0.264*10* inputReader->MP* inputReader->massScale*sigs[i]))+log(1+sigp[i]/sigs[i]/log(0.038*10* inputReader->MP* inputReader->massScale *sigs[i])));
        }
    }
    else{
        //loop over the channels and calculate x (has no s-wave)
        for(int i = 0; i < inputReader->Nf-1; i++){
            answer.push_back(2.0*log((2+c)*c*0.145*0.264*10* inputReader->MP * inputReader->massScale*sigp[i])-0.5*log(log((2+c)*c*0.145*0.264*10* inputReader->MP* inputReader->massScale*sigp[i])));
        }
    
    }
    //this will also include the SM case and since its bigger than 1000, it will always be the eqset
    answer.push_back(10000); 
    return answer; 

}*/


//This will return the value of Y or W based on the current Y's or W's and current Yeq's or Weq's (uses Euler method) (using wave approximation)
void NumericalCalculator::fnEuler(vector<double>& inset, vector<double> eqset,double t){
    vector<double> func=f(inset,eqset,t);
    //Using Lambda Expression, ignore warning about c++11, we are assuming that.
    transform(inset.begin(),inset.end(),func.begin(),inset.begin(), [=](double y,double f){return y+inputReader->step(t)*f;});

}

//Uses RK4 (Runge Kutta 4) (loop up the equations on wikipedia or something)
void NumericalCalculator::fnRK4(vector<double>& inset, vector<double> eqset,double t){

    vector<double> tempin;
    vector<double> tempeq1 =genEqSet(t*inputReader->massScale/(inputReader->massScale+t*inputReader->step(t)/2.0));
    vector<double> tempeq2 =genEqSet(t*inputReader->massScale/(inputReader->massScale+t*inputReader->step(t)));
    tempin.resize(inputReader->Nf);

    vector<double> k1=f(inset,eqset,t);

    transform(inset.begin(),inset.end(),k1.begin(),tempin.begin(),[=](double y,double k)->double{return y+inputReader->step(t)/2.0*k;});
    vector<double> k2=f(tempin,tempeq1,t*inputReader->massScale/(inputReader->massScale+t*inputReader->step(t)/2.0));

    transform(inset.begin(),inset.end(),k2.begin(),tempin.begin(),[=](double y,double k)->double{return y+inputReader->step(t)/2.0*k;});
    vector<double> k3=f(tempin,tempeq1,t*inputReader->massScale/(inputReader->massScale+t*inputReader->step(t)/2.0));

    transform(inset.begin(),inset.end(),k3.begin(),tempin.begin(),[=](double y,double k)->double{return y+inputReader->step(t)*k;});
    vector<double> k4=f(tempin,tempeq2,t*inputReader->massScale/(inputReader->massScale+t*inputReader->step(t)));

    for(int i=0;i < inputReader->Nf;i++) inset[i]+=inputReader->step(t)/6.0*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i]);    

}   

//Uses Linear Multistep Method (Adams-bashforth order 5). The prevF variable will contain f(t_(n-1),y_(n-1),...) to f(t_(n-4),y_(n-4),...).
//The oldest will then be replaced with the newly calculated one. The variable prevFIndex will keep track of the oldest variable. 
//Taking this mod 4 will give (n-4) to (n-1) by incrementing the index.
void NumericalCalculator::fnLMM(vector<double>& inset, vector<double> eqset, vector<vector<double>* >& prevF, int& prevFIndex, double t){

    vector<double> newF = f(inset,eqset,t);

    for(int i=0;i < inputReader->Nf;i++) inset[i]+=inputReader->step(t)*(251./720. * prevF[prevFIndex]->at(i)-637./360. * prevF[(prevFIndex+1)%4]->at(i) +109./30. * prevF[(prevFIndex+2)%4]->at(i)-1387./360. * prevF[(prevFIndex+3)%4]->at(i)+1901./720. * newF[i]); 

    prevF[prevFIndex]=&newF;
    prevFIndex=(prevFIndex+1)%4;

}

//This function goes through and actually solves the equation. It also fills out the omega vector with the relic densities of each dark matter species.
void NumericalCalculator::solve(bool write){


    clock_t t1, t2; //This was there originally just for timing stuff
    t1=clock();

    vector<double> startingTVec;
    startingTVec.resize(inputReader->startingxVec.size());

    //DEBUG
    //inputReader->startingxVec[3]=10000;
    //inputReader->startingxVec[1]=10000;
    //inputReader->startingxVec[2]=10000;

    cout<<"Intervals for small step size"<< endl;
    for(int i=0; i<inputReader->intervals.size();i++){
        cout << inputReader->intervals[i]->at(0) << " to " << inputReader->intervals[i]->at(1)<< endl;
    }

    transform(inputReader->startingxVec.begin(),inputReader->startingxVec.end(),startingTVec.begin(), [=](double x){return inputReader->massScale/x;});

    double startingx = *min_element(inputReader->startingxVec.begin(),inputReader->startingxVec.end());
    double startingT = inputReader->massScale/startingx;
    if(write){
        cout << "Starting x val: " << startingx << endl;
        cout << "Starting T val: " << startingT << endl;
    }
   
    for(int i=0; i< inputReader->Nf;i++) cout <<i<< " starting x and t: "<< inputReader->startingxVec[i]<< " " << startingTVec[i] << endl;

    //temp and x will both be used
    double x0 = startingx; 
    double Temp = startingT;
    
    vector<double> eqset=genEqSet(Temp); //Equilibrium Y values at different temperatures
    vector<double> inset=eqset; //Actual Y values at different temperatures
    
    //When ever x is > nextSave, then we save.
    double nextSave =1.0;
    double saveInterval=1.0; //the interval in x between each save

    vector<vector<double>* > prevF;
    prevF.resize(4);
    for(int i=0;i<4;i++) prevF[i]=new vector<double>(inputReader->Nf,0);
    int prevFInd=0;

    cout << "method in use: "<< inputReader->method<< endl;

    //DEBUG
    cout<< inputReader->Nsteps<< endl;
    
    cout<<"starting main Loop"<<endl;
    //main loop for the numerical method
    for(int ti = 0; ti < inputReader->Nsteps; ti++){ //loop over temperature increments
         
        //DEBUG
        //cout<<ti<< " " << Temp << " " << x0 <<endl;

        //Create the eq set for the new x
        eqset=genEqSet(Temp);

    
        /*if(Temp<startingTVal){
            //if we are below the starting temperature, go through the calculation
            switch(inputReader->method){
                case 0: fnEuler(inset,eqset,Temp);
                        break;
                case 1: fnRK4(inset,eqset,Temp);
                        break;
                case 2: fnLMM(inset,eqset,prevF,prevFInd,Temp);
                        break;
            }
            //make sure to copy over the radiation values
            inset[inputReader->Nx] = eqset[inputReader->Nx];  
        }
        else{
            //if we are above the starting temperature, just set Y=Yeq or W=Weq
            inset=eqset;
        }*/

        switch(inputReader->method){
            case 0: fnEuler(inset,eqset,Temp);
                    break;
            case 1: fnRK4(inset,eqset,Temp);
                    break;
            case 2: fnLMM(inset,eqset,prevF,prevFInd,Temp);
                    break;
        }
        //for the fields that are still in equilibrium, we just make the inset = eqset (this will alwasy happen for SM)
        for(int i=0; i<inputReader->Nf; i++) if(Temp>startingTVec[i]) inset[i] = eqset[i];
        

        //check if NaN
        if( any_of(inset.begin(),inset.end(), [](double x){return x!=x;})){
            if(write) cout<< "NaN for a Y or W value on step ti =" << ti << "; temperature of "<< Temp<< endl;
            ti = inputReader->Nsteps;
            nextSave = 0;//Forces Save
        }

        //decide whether to write the point to the file or not
        if(write && x0>nextSave){
            outputWriter->write(Temp,eqset,inset,inputReader->useW);
            outputWriter->keepY(inputReader->useW,inset); // This is in case of a wild NaN near the end
            nextSave=x0+saveInterval;
            
        }
        
        //Increment stuff
        x0+=inputReader->step(Temp);
        Temp = inputReader->massScale/x0;
        if(x0 >1000) break;
    }
   
    //check if Y was NaN and if not, use that for omega
    if( all_of(inset.begin(),inset.end(), [](double x){return x==x;})){
        outputWriter->keepY(inputReader->useW,inset);
    }

    //calculate omega
    omega.resize(inputReader->Nx); 
    transform(outputWriter->lastY.begin(),outputWriter->lastY.end(),inputReader->mxset.begin(),omega.begin(), [=](double in1,double m)->double{ return m*s0*in1/rhoC;});

    
    //close the openned files
    if(write) outputWriter->close();

    t2=clock();
    //time taken
    if(write) cout << "Time Taken: " << ((float)t2-(float)(t1)) / CLOCKS_PER_SEC << endl; 

}
