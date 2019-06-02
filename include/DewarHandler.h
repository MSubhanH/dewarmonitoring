//============================================================================
// Name        : DewarHandler.h
// Author      : Muhammad Subhan Hameed
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <dislin.h>
#include "discpp.h"

#ifndef DEWARHANDLER_H
#define DEWARHANDLER_H

#define WINDOWSIZE 1000
#define NOOFQUANTS 11

#if defined _MSC_VER
  typedef unsigned __int64 uint64;
#else
  typedef uint64_t uint64;
#endif

class DewarHandler{

public:
	double S2;							// Stores current value of Stage 2 temperature
	double P;							// Stores current value of Pressure
	double ET;							// Stores current value of Environment temperature

	double S2_series[WINDOWSIZE];		// Stores Stage 2 values for complete WINDOWSIZE, initialized with predefined values
	double P_series[WINDOWSIZE];		// Stores Pressure values for complete WINDOWSIZE, initialized with predefined values
	double ET_series[WINDOWSIZE];		// Stores Environment values for complete WINDOWSIZE, initialized with predefined values

	float f1_series[WINDOWSIZE];		// Stores normalized values of Pressure for complete WINDOWSIZE, initialized with zero
	float f2_series[WINDOWSIZE];		// Stores normalized values of Environment temperature for complete WINDOWSIZE, initialized with zero

	double resS2AndModel_series[WINDOWSIZE];	// Stores values of residual between Stage 2 and model for complete WINDOWSIZE, initialized with zero

	float eqn_series[WINDOWSIZE];		// Stores values of model for complete WINDOWSIZE, initialized with zero
	float dist_series[WINDOWSIZE];		// Stores values of distance of current PD with closest reference PD for complete WINDOWSIZE, initialized with zero
	float testQuantity;					// Stores cuurent value of distance of current PD with closest reference PD

	bool alarms_series[WINDOWSIZE];		// Stores values of alarms for complete WINDOWSIZE, initialized with false (no alarms)

	int stabilityLevel;					// Stores current value of Stability Level (10 = Stable, 5 = Semistable or 1 = Unstable)
	char* stabilityLabel;				// Stores current value of Stability Label ('Stable', 'Semistable' or 'Unstable')

	int sampleIndex;					// Stores sample index

	struct SeqQuantEst{					
		int n[NOOFQUANTS];
		int m[NOOFQUANTS];
		float seqQuantEst[NOOFQUANTS];
	};
	SeqQuantEst dewarQuantEst;			// Stores current sequential qunatile estimate (11 quantiles in 0.01 to 0.99)

	Dislin g, r;						// Objects of DISLIN library to plot Dewar Data Monitoring Window and Dewar Probability Distributions Window

	struct Features{
		float f1;						// Normalized pressure
		float f2;						// Normalized environment temperature
		float f4;						// Normalized distance of environment temperature from threshold
		float f5;						// Normalized environment temperature without seasonal variation
		float y0;						// Stage 2 temperature
	};
	Features features;					// Stores the current values of normalized features for the model

	DewarHandler();						// Constructor -> calls initDewar()

	void initDewar();					// Initializes dewar properties 

	void updateDewarObj();				// Updates dewar object , called at every time sample in the main routine

	void updateData();					// Updates dewar data properties, called in updateDewarObj()
	void updateStabilityLevel();		// Updates dewar stability level and label, called in updateDewarObj()
	void updateDataPlot();				// Updates Dewar Data Monitoring Window, called in updateDewarObj()
	void updateResPlot();				// Updates Dewar Dewar Probability Distributions Window, called in updateDewarObj()

	void plotData(double series1[], char* color1, float series2[], char* color2, double series3[], char* color3, float series4[], char* color4, bool series5[], char* color5);	// Plots dewar data properties, called in updateDataPlot()
	void plotResData(double tau[], double centroid1[], double centroid2[], double centroid3[], double centroid4[], double centroid5[], float seqQuantEst[]);  // Plots dewar probaility distribution and reference distributions, called in updateResPlot()

	void normalizeFeatures();			// Normalizes features, called in updateData()
	void seqQuantEstimator(double res, double quantEst, int n, int m, double tau, int N, int indexOfQuant);  // Estimates dewar probability distribution (quantile estimate), called in updateData() 

	void initializeAsZero(int arrayOfVals[], int size);
	void initializeAsZero(float arrayOfVals[], int size);
	void initializeAsZero(bool arrayOfVals[], int size);
	void initializeAsZero(double arrayOfVals[], int size);
	void initializeAsRef(double arrayOfVals[], double arrayOfRefVals[], int size);

	void shiftSeries(double arrayOfVals[], int size, double newVal);
	void shiftSeries(float arrayOfVals[], int size, double newVal);
	void shiftSeries(bool arrayOfVals[], int size, double newVal);

	int isInfinite(double);

};

#endif
