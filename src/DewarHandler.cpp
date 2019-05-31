#include <iostream>
#include <math.h>
#include <cmath>
#include <algorithm>

#include "DewarHandler.h"

#include "pressureVals.h"
#include "etVals.h"
#include "stage2Vals.h"

#include <conio.h>

using namespace::std;

DewarHandler::DewarHandler(){
	initDewar();
}


void DewarHandler::initDewar(){

	initializeAsZero(dewarQuantEst.n, NOOFQUANTS);
	initializeAsZero(dewarQuantEst.m, NOOFQUANTS);
	initializeAsZero(dewarQuantEst.seqQuantEst, NOOFQUANTS);

	initializeAsZero(f1_series, WINDOWSIZE);
	initializeAsZero(f2_series, WINDOWSIZE);

	initializeAsZero(resS2AndModel_series, WINDOWSIZE);

	initializeAsZero(eqn_series, WINDOWSIZE);
	initializeAsZero(dist_series, WINDOWSIZE);

	initializeAsZero(alarms_series, WINDOWSIZE);

	initializeAsRef(P_series, pressureVals, WINDOWSIZE);
	initializeAsRef(ET_series, etVals, WINDOWSIZE);
	initializeAsRef(S2_series, stage2Vals, WINDOWSIZE);

	features.f1 = 0;
	features.f2 = 0;
	features.f4 = 0;
	features.f5 = 0;
	features.y0 = 0;

	stabilityLevel = 0;
	stabilityLabel = "Stable";
	sampleIndex = 0;
}



void DewarHandler::updateDewarObj(){

	updateData();
	updateStabilityLevel();
	updateDataPlot();
	updateResPlot();
}

void DewarHandler::updateData(){

		double centroid[5][11] = {	{2.30046776476642,2.04421545918589,1.90359863595407,1.76796015098872,1.61777484665221,1.46956740874188,1.30567087032978,1.15977351692190,0.953300304551104,0.598386737013665,0.0732726375841887},
									{-1.13787358751092,-1.66375010567253,-1.94719643813235,-2.14675064107972,-2.32067574041200,-2.47394482486541,-2.63286273846755,-2.82909459802170,-3.05419449375825,-3.26389945613897,-3.49732521768526},
									{-0.0828521335859643,-0.576107230721373,-0.834787206967767,-1.02999971721861,-1.19194016344773,-1.38266549783672,-1.55536747447899,-1.77055170658599,-2.01912846760733,-2.30097842377618,-2.75019427084803},
									{0.768636509980674,0.454651641983256,0.235561493882798,0.0465315518351459,-0.165677076625891,-0.391642627173204,-0.615682227945935,-0.825430457179592,-1.07125080489373,-1.31511461687061,-1.71209433354796},
									{1.90022766850592,1.56895899400586,1.29051197917614,0.990982347441981,0.745794365112853,0.537457518376578,0.285824061004898,0.0685003574505601,-0.202924862060755,-0.529224790570599,-1.09581654537788}	};

		normalizeFeatures();
//		cout<<features.f1<<endl<<features.f2<<endl<<features.f4<<endl<<features.f4<<endl<<features.f5;

		double theta_m_sol[5] = {5.15634525286613,0,16.4699636273839,5.52409305859874,7.36962642092408};
		double eqn = theta_m_sol[0]*features.f1*(1+features.f4) + theta_m_sol[4]*features.f5*(1-features.f4) + (theta_m_sol[3]*features.f4) + theta_m_sol[2];

		double  resS2AndModel = features.y0 - eqn;

		testQuantity = 0;

		if(resS2AndModel > 0){

			double tau[11] = {0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99};

			int N = 10;

            for(int k = 0; k < NOOFQUANTS; k++){
				seqQuantEstimator(resS2AndModel, dewarQuantEst.seqQuantEst[k], dewarQuantEst.n[k], dewarQuantEst.m[k], tau[k], N, k);
			}

			int K = 5;
            float minDistFromRef = 1000;

            for (int i = 0; i < K; i++){

				double normOfDiff = 0;

				for(int k = 2; k < NOOFQUANTS; k++)
					normOfDiff = normOfDiff + pow((dewarQuantEst.seqQuantEst[k] - centroid[i][k]), 2);

				normOfDiff = sqrt(normOfDiff);

				if(normOfDiff < minDistFromRef)
					minDistFromRef = normOfDiff;
			}

            testQuantity = sqrt(minDistFromRef);
		}

		else
            testQuantity = 0;

		float thresh = 1.5;

		bool alarm = 0;
        if(testQuantity > thresh)
          alarm = 1;
        else
          alarm = 0;

		shiftSeries(S2_series, WINDOWSIZE, S2);
		shiftSeries(P_series, WINDOWSIZE, P);
		shiftSeries(ET_series, WINDOWSIZE, ET);
		shiftSeries(f1_series, WINDOWSIZE, features.f1 + theta_m_sol[2]);
		shiftSeries(f2_series, WINDOWSIZE, features.f2 + theta_m_sol[2]);
		shiftSeries(resS2AndModel_series, WINDOWSIZE, resS2AndModel);
		shiftSeries(eqn_series, WINDOWSIZE, eqn);
		shiftSeries(dist_series, WINDOWSIZE, testQuantity);
		shiftSeries(alarms_series, WINDOWSIZE, alarm);

	}

void DewarHandler::updateStabilityLevel(){

		double estQuantStable[10] = {2.0791, 1.8741, 1.6982, 1.5225, 1.3223, 1.1245, 0.8948, 0.6539, 0.3670};
		double estQuantUnstable[10] = {5.7810, 4.4364, 3.6016, 3.0694, 2.7279, 2.4490, 2.1431, 1.8422, 1.5253, 0.9924};
		double estQuantSemiStable[10] = {4.1384, 3.2577, 2.7378, 2.3838, 2.1252, 1.8856, 1.6338, 1.3685, 1.0896, 0.6797};

		double normDistFromStable = 0;
		double normDistFromUnstable = 0;
		double normDistFromSemiStable = 0;

		double normOfDiff = 0;
		for(int k = 2; k < NOOFQUANTS; k++)
			normOfDiff = normOfDiff + pow((dewarQuantEst.seqQuantEst[k] -  estQuantStable[k-2]), 2);
		normDistFromStable  = sqrt(normOfDiff);

		normOfDiff = 0;
		for(int k = 2; k < NOOFQUANTS; k++)
			normOfDiff = normOfDiff + pow((dewarQuantEst.seqQuantEst[k] -  estQuantUnstable[k-2]), 2);
		normDistFromUnstable  = sqrt(normOfDiff);

		normOfDiff = 0;
		for(int k = 2; k < NOOFQUANTS; k++)
			normOfDiff = normOfDiff + pow((dewarQuantEst.seqQuantEst[k] -  estQuantSemiStable[k-2]), 2);
		normDistFromSemiStable  = sqrt(normOfDiff);


		float dist = 1000;

		if(normDistFromStable < dist){
		   stabilityLevel =  10;
		   stabilityLabel = "Stable";
           dist = normDistFromStable;
		}

        if(normDistFromUnstable < dist){
		   stabilityLevel =  1;
		   stabilityLabel = "Unstable";
           dist = normDistFromUnstable ;
		}

        if(normDistFromSemiStable < dist){
	   	   stabilityLevel =  5;
  		   stabilityLabel = "Semistable";
		   dist = normDistFromSemiStable;
		}

	}

void DewarHandler::updateDataPlot(){

	plotData(S2_series, "blue", eqn_series, "red", resS2AndModel_series, "green", dist_series, "orange");

}

void DewarHandler::updateResPlot(){
		
	double centroid1[11] =  {2.30046776476642,2.04421545918589,1.90359863595407,1.76796015098872,1.61777484665221,1.46956740874188,1.30567087032978,1.15977351692190,0.953300304551104,0.598386737013665,0.0732726375841887};
	double centroid2[11] =	{-1.13787358751092,-1.66375010567253,-1.94719643813235,-2.14675064107972,-2.32067574041200,-2.47394482486541,-2.63286273846755,-2.82909459802170,-3.05419449375825,-3.26389945613897,-3.49732521768526};
	double centroid3[11] =	{-0.0828521335859643,-0.576107230721373,-0.834787206967767,-1.02999971721861,-1.19194016344773,-1.38266549783672,-1.55536747447899,-1.77055170658599,-2.01912846760733,-2.30097842377618,-2.75019427084803};
	double centroid4[11] =	{0.768636509980674,0.454651641983256,0.235561493882798,0.0465315518351459,-0.165677076625891,-0.391642627173204,-0.615682227945935,-0.825430457179592,-1.07125080489373,-1.31511461687061,-1.71209433354796};
	double centroid5[11] =	{1.90022766850592,1.56895899400586,1.29051197917614,0.990982347441981,0.745794365112853,0.537457518376578,0.285824061004898,0.0685003574505601,-0.202924862060755,-0.529224790570599,-1.09581654537788};
	
	double tau[11] = {0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99};

	plotResData(tau, centroid1, centroid2, centroid3,centroid4,centroid5, dewarQuantEst.seqQuantEst);

}



void DewarHandler::plotData(double series1[], char* color1, float series2[], char* color2, double series3[], char* color3, float series4[], char* color4){

	  float x1[WINDOWSIZE], y1[WINDOWSIZE];
	  float x2[WINDOWSIZE], y2[WINDOWSIZE];
	  float x3[WINDOWSIZE], y3[WINDOWSIZE];
	  float x4[WINDOWSIZE], y4[WINDOWSIZE];
	  float x5[WINDOWSIZE], y5[WINDOWSIZE];
	  
	  for (int i = 0; i < WINDOWSIZE; i++){
		 
		  x1[i] = i;
		  y1[i] = (float)series1[i];

		  x2[i] = i;
		  y2[i] = series2[i];

		  x3[i] = i;
		  y3[i] = (float)series3[i];

		  x4[i] = i;
		  y4[i] = series4[i];
	  }
		 
	  g.metafl ("png");
	  g.setfil("DewarDataPlot.png");
	  g.filmod("delete");
	  g.scrmod ("revers");

	  g.disini ();

	  g.pagera ();
	  g.complx ();
	  g.axspos (450, 1800);
	  g.axslen (2200, 1200);

	  g.name   ("Sample Index","x");
	  g.name   ("Absolute Units","y");

	  g.labdig (-1, "x");
	  g.ticks  (9, "x");
	  g.ticks  (10, "y");

	  g.titlin ("Dewar Monitoring", 1);
	 
	  int ic;
	  ic = g.intrgb (0.95,0.95,0.95);
	  g.axsbgd (ic);
	  g.graf   (0.0, 1000, 0.0, 100.0, -10, 50, 0.0, 5.0);
	  g.setrgb (0.7, 0.7, 0.7);
	  g.grid   (1, 1);

	  g.color  ("fore");
	  g.height (50);
	  g.title  ();

	  g.color  (color1);
	  g.curve  (x1, y1, WINDOWSIZE);

	  g.color  (color2);
	  g.curve  (x2, y2, WINDOWSIZE);

	  g.color  (color3);
	  g.curve  (x3, y3, WINDOWSIZE);

	  g.color  (color4);
	  g.curve  (x4, y4, WINDOWSIZE);
  
	  g.errmod("all", "off");  
	  g.disfin ();
	  	 
	  
}

void DewarHandler::plotResData(double tau[], double centroid1[], double centroid2[], double centroid3[], double centroid4[], double centroid5[], float seqQuantEst[]){

	  float y[NOOFQUANTS];

	  float x1[NOOFQUANTS];
	  float x2[NOOFQUANTS];
	  float x3[NOOFQUANTS];
	  float x4[NOOFQUANTS];
	  float x5[NOOFQUANTS];
	  	
	  float quantEst[NOOFQUANTS];

	  for (int i = 0; i < NOOFQUANTS; i++){
		 
		  y[i] = (float)tau[i];

		  x1[i] = (float)centroid1[NOOFQUANTS - i - 1];
		  x2[i] = (float)centroid2[NOOFQUANTS - i - 1];
		  x3[i] = (float)centroid3[NOOFQUANTS - i - 1];
		  x4[i] = (float)centroid4[NOOFQUANTS - i - 1];
		  x5[i] = (float)centroid5[NOOFQUANTS - i - 1];

		  quantEst[i] = (float)seqQuantEst[NOOFQUANTS - i - 1];
	  }
		 
	  r.metafl ("png");
	  r.setfil("Distribution Plot.png");
	  r.filmod("delete");
	  r.scrmod ("revers");

	  r.disini ();

	  r.pagera ();
	  r.complx ();
	  r.axspos (450, 1800);
	  r.axslen (2200, 1200);

	  r.name   ("Residual","x");
	  r.name   ("Probability","y");

	  r.labdig (-1, "x");
	  r.ticks  (9, "x");
	  r.ticks  (10, "y");

	  r.titlin ("Cumulative Probability Distribution Plot", 1);
	 
	  int ic;
	  ic = g.intrgb (0.95,0.95,0.95);
	  r.axsbgd (ic);
	  r.graf   (-5.0, 10.0, -5.0, 1, 0, 1, 0.0, 0.1);
	  r.setrgb (0.7, 0.7, 0.7);
	  r.grid   (1, 1);

	  r.color  ("fore");
	  r.height (50);
	  r.title  ();

	  r.color  ("orange");
	  r.curve  (x1,y,NOOFQUANTS);

	  r.color  ("orange");
	  r.curve  (x2,y,NOOFQUANTS);
	
	  r.color  ("orange");
	  r.curve  (x3,y, NOOFQUANTS);

	  r.color  ("orange");
	  r.curve  (x4,y, NOOFQUANTS);

	  r.color  ("orange");
	  r.curve  (x5,y, NOOFQUANTS);

	  r.color  ("magenta");
	  r.curve  (quantEst,y,NOOFQUANTS);

	  r.errmod("all", "off");  
	  r.disfin ();


}

void DewarHandler::normalizeFeatures(){

		double f11[WINDOWSIZE];
		std::pair<double*, double*> minmax = std::minmax_element(std::begin(P_series), std::end(P_series));
		for(int i = 0; i < WINDOWSIZE; i++){
			f11[i] = (P_series[i] - *minmax.first) / (*minmax.second - *minmax.first);
			f11[i] = log(f11[i]) ;
		}
		for(int i = 0; i < WINDOWSIZE; i++){
			if(isInfinite(f11[i]))
				f11[i] = 2.2600e-04;
		}
		double f11_val = (P - min(*minmax.first, P)) / (max(*minmax.second, P) - min(*minmax.first, P));
		minmax = std::minmax_element(std::begin(f11), std::end(f11));
		features.f1 = (log(f11_val) - min(*minmax.first, f11_val)) / (max(*minmax.second, f11_val) - min(*minmax.first, f11_val));
		if(isInfinite(features.f1))
			features.f1 = 0;

		minmax = std::minmax_element(std::begin(ET_series), std::end(ET_series));
		features.f2 = (ET - *minmax.first) / (*minmax.second - *minmax.first);


		double f44[WINDOWSIZE];
		for(int i = 0; i < WINDOWSIZE; i++){
			f44[i] = 25 - ET_series[i];
		}
		minmax = std::minmax_element(std::begin(f44), std::end(f44));
		features.f4 = ((25 - ET) - min(*minmax.first, (25 - ET))) / (max(*minmax.second, (25 - ET)) - min(*minmax.first, (25 - ET)));


		double eqnpol_series[WINDOWSIZE];
		double coeff[4] = { 0.0000,   -0.0000,   -0.0000,   23.3071};
		double f55[WINDOWSIZE];
		for(int i = 0; i < WINDOWSIZE; i++){
			eqnpol_series[i] = coeff[0]*(i^3) + coeff[1]*(i^2) + coeff[2]*i + coeff[3];
			f55[i] = ET_series[i] - eqnpol_series[i];
		}
		double eqn_val = coeff[0]*(pow(ET,3)) + coeff[1]*(pow(ET,2)) + coeff[2]*ET + coeff[3];
		minmax = std::minmax_element(std::begin(f55), std::end(f55));
		features.f5 = ((ET - eqn_val) - min(*minmax.first, (ET - eqn_val))) / (max(*minmax.second, (ET - eqn_val)) - min(*minmax.first, (ET - eqn_val)));


		features.y0 = S2;
        if(features.y0 == 0)
           features.y0 = 10;

	}

void DewarHandler::seqQuantEstimator(double res, double quantEst, int n, int m, double tau, int N, int indexOfQuant){

	double s = 0.02;

	if(res >= quantEst)
		n = n + 1;
    else
		m = m + 1;


	if(n > N*tau){
		quantEst = quantEst + s;
		n = 0;
		m = 0;
	}

	else if(m > N*(1 - tau)){
		quantEst = quantEst - s;
		n = 0;
		m = 0;
	}

	else if(m + n == N){
		n = 0;
		m = 0;
	}


	dewarQuantEst.n[indexOfQuant] = n;
	dewarQuantEst.m[indexOfQuant] = m;
	dewarQuantEst.seqQuantEst[indexOfQuant] = quantEst;

}


void DewarHandler::initializeAsZero(int arrayOfVals[], int size){

	for(int i = 0; i < size; i++){
		arrayOfVals[i] = 0;
	}
}
void DewarHandler::initializeAsZero(float arrayOfVals[], int size){

	for(int i = 0; i < size; i++){
		arrayOfVals[i] = 0;
	}
}
void DewarHandler::initializeAsZero(bool arrayOfVals[], int size){

	for(int i = 0; i < size; i++){
		arrayOfVals[i] = 0;
	}
}
void DewarHandler::initializeAsZero(double arrayOfVals[], int size){

	for(int i = 0; i < size; i++){
		arrayOfVals[i] = 0;
	}
}
void DewarHandler::initializeAsRef(double arrayOfVals[], double arrayOfRefVals[], int size){
	for(int i = 0; i < size; i++){
		arrayOfVals[i] = arrayOfRefVals[i];
	}
}

void DewarHandler::shiftSeries(double arrayOfVals[], int size, double newVal){

		for(int i = 0; i < size-1; i++){
			arrayOfVals[i] = arrayOfVals[i+1];
		}

		arrayOfVals[size-1] = newVal;
	}
void DewarHandler::shiftSeries(float arrayOfVals[], int size, double newVal){

		for(int i = 0; i < size-1; i++){
			arrayOfVals[i] = arrayOfVals[i+1];
		}

		arrayOfVals[size-1] = newVal;
	}
void DewarHandler::shiftSeries(bool arrayOfVals[], int size, double newVal){

		for(int i = 0; i < size-1; i++){
			arrayOfVals[i] = arrayOfVals[i+1];
		}

		arrayOfVals[size-1] = newVal;
	}

int DewarHandler::isInfinite(double x)
{
    union { uint64 u; double f; } ieee754;
    ieee754.f = x;
    return ( (unsigned)(ieee754.u >> 32) & 0x7fffffff ) == 0x7ff00000 &&
           ( (unsigned)ieee754.u == 0 );
}
