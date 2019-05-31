//============================================================================
// Name        : DewarQuantileChangeDetection.cpp
// Author      : Muhammad Subhan Hameed
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <iomanip>
#include <conio.h>

#include <dislin.h>

#include "DewarHandler.h"

#include "pressureData10.h"
#include "etData10.h"
#include "stage2Data10.h"

using namespace std;

double getPressure(int);
double getStage2(int);
double getET(int);

int main(void){

	DewarHandler obj = DewarHandler();
	
	cout<<"Press ESC to run the dewar quantile based change detection algorithm"<<endl<<endl;

	cout<<std::left<<std::setfill(' ')<<std::setw(15)<<"Stage2";
	cout<<std::left<<std::setfill(' ')<<std::setw(15)<<"Model";
	cout<<std::left<<std::setfill(' ')<<std::setw(15)<<"Residual";
	cout<<std::left<<std::setfill(' ')<<std::setw(15)<<"TestQuantity";
	cout<<std::left<<std::setfill(' ')<<std::setw(15)<<"Alarm";
	cout<<std::left<<std::setfill(' ')<<std::setw(15)<<"StabilityLevel"<<endl;

	
	//cout<<std::left<<std::setfill(' ')<<std::setw(15)<<"F1";
	//cout<<std::left<<std::setfill(' ')<<std::setw(15)<<"F2";
	//cout<<std::left<<std::setfill(' ')<<std::setw(15)<<"F4";
	//cout<<std::left<<std::setfill(' ')<<std::setw(15)<<"F5";
	//cout<<std::left<<std::setfill(' ')<<std::setw(15)<<"Y0"<<endl;

	while(1){

		int keyBoardEntry = getch();

	
		if(keyBoardEntry == 27){

			// Replace the code for interface with sensor data here

			obj.S2 = getStage2(obj.sampleIndex);
			obj.P = getPressure(obj.sampleIndex);
			obj.ET = getET(obj.sampleIndex);

			////////////////////////////////////////////////////////

			obj.updateDewarObj();
			obj.sampleIndex = obj.sampleIndex + 1;


			cout<<std::left<<std::setfill(' ')<<std::setw(15)<<obj.S2;
			cout<<std::left<<std::setfill(' ')<<std::setw(15)<<obj.eqn_series[999];
			cout<<std::left<<std::setfill(' ')<<std::setw(15)<<obj.resS2AndModel_series[999];
			cout<<std::left<<std::setfill(' ')<<std::setw(15)<<obj.testQuantity;
			cout<<std::left<<std::setfill(' ')<<std::setw(15)<<obj.alarms_series[999];
			cout<<std::left<<std::setfill(' ')<<std::setw(15)<<obj.stabilityLabel<<endl;

			//cout<<std::left<<std::setfill(' ')<<std::setw(15)<<obj.features.f1;
			//cout<<std::left<<std::setfill(' ')<<std::setw(15)<<obj.features.f2;
			//cout<<std::left<<std::setfill(' ')<<std::setw(15)<<obj.features.f4;
			//cout<<std::left<<std::setfill(' ')<<std::setw(15)<<obj.features.f5;
			//cout<<std::left<<std::setfill(' ')<<std::setw(15)<<obj.features.y0<<endl;

		}

	}


	return 0;
}


double getPressure(int sampleIndex){

	double pressure;
	pressure = pressureData10[sampleIndex];

	return pressure;
}

double getStage2(int sampleIndex){

	double stage2;
	stage2 = stage2Data10[sampleIndex];

	return stage2;
}

double getET(int sampleIndex){

	double et;
	et = etData10[sampleIndex];

	return et;
}
