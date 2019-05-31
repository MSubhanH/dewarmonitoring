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
	double S2;
	double P;
	double ET;

	double S2_series[WINDOWSIZE];
	double P_series[WINDOWSIZE];
	double ET_series[WINDOWSIZE];

	float f1_series[WINDOWSIZE];
	float f2_series[WINDOWSIZE];

	double resS2AndModel_series[WINDOWSIZE];

	float eqn_series[WINDOWSIZE];
	float dist_series[WINDOWSIZE];
	float testQuantity;

	bool alarms_series[WINDOWSIZE];

	int stabilityLevel;
	char* stabilityLabel;

	int sampleIndex;

	struct SeqQuantEst{
		int n[NOOFQUANTS];
		int m[NOOFQUANTS];
		float seqQuantEst[NOOFQUANTS];
	};
	SeqQuantEst dewarQuantEst;

	Dislin g, r;

	struct Features{
		float f1;
		float f2;
		float f4;
		float f5;
		float y0;
	};
	Features features;

	DewarHandler();

	void initDewar();

	void updateDewarObj();

	void updateData();
	void updateStabilityLevel();
	void updateDataPlot();
	void updateResPlot();

	void plotData(double series1[], char* color1, float series2[], char* color2, double series3[], char* color3, float series4[], char* color4);
	void plotResData(double tau[], double centroid1[], double centroid2[], double centroid3[], double centroid4[], double centroid5[], float seqQuantEst[]);

	void normalizeFeatures();
	void seqQuantEstimator(double res, double quantEst, int n, int m, double tau, int N, int indexOfQuant);

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
