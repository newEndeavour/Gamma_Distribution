/*
  File:         Gamma_Distribution.cpp
  Version:      0.0.1
  Date:         30-Jan-2019
  Revision:     30-Jan-2019
  Author:       Jerome Drouin (jerome.p.drouin@gmail.com)

  Editions:	Please go to Gamma_Distribution.h for Edition Notes.

  https://en.wikipedia.org/wiki/Gamma_distribution
  A special thank you to:
  http://www.mymathlib.com/functions/probability/gamma_distribution.html


  Gamma_Distribution implements the Gamma Distribution with shape Alpha and 
  Rate Beta or Scale Theta. 

  Copyright (c) 2018-2019 Jerome Drouin  All rights reserved.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/


#include "Arduino.h"
#include "Gamma_Distribution.h"
#include <Gamma_Function.h>


// Constructor /////////////////////////////////////////////////////////////////
// Function that handles the creation and setup of instances
// Two Characterisations:
// Shape & Rate		: Alpha & Beta (Theta=0), Note: Beta = 1/Theta
// Shape & Scale	: Alpha & Theta (Beta=0), Note: Alpha = k
// The object enables two (mutually exclusive) characterisations:
//	- Shape (Alpha) & Rate (Beta) 	:  Gamma_Distribution(Shape>0,Rate>0,0)
//					 = Gamma_Distribution(Alpha>0,Beta>0,0) 
//	- Shape (Alpha) & Scale (Theta)	:  Gamma_Distribution(Shape>0,0,Scale)
//					 = Gamma_Distribution(Alpha>0,0,Theta>0) 
// Note: An attempt to construct an object with positive Alpha, Beta AND Theta 
// creates an indetermination and returns an error flag.		
//		
Gamma_Distribution::Gamma_Distribution(double _Alpha, double _Beta, double _Theta)
{

	// Object parameter's error handling
	error = 1;

	if (_Alpha<=0)	
		error = -2;				// Invalid definition domain
	else {
		if (_Beta>0 && _Theta>0) error = -3;	// Only 1 parameter must be >0, the other must be =0
		else {
			if (_Beta>0 && _Theta<=0) error = 1;
			if (_Beta<=0 && _Theta>0) error = 1;
			if (_Beta<=0 && _Theta<=0) error = -4;
		}
	}

	//Set initial values	
	Alpha			= _Alpha;		// 
	Beta			= _Beta;		// 
	Theta			= _Theta;		// 

}


// Public Methods //////////////////////////////////////////////////////////////
//Probability Density Function
double Gamma_Distribution::GetPDF(double x)
{
	if (error<0)
		return error;

	if (x <= 0.0) 
		return 0.0;
   	
	//Shape and Rate
	if (Beta>0) {
		if (Alpha<=Gamma_Function_Max_Arg())
			return pow(Beta,Alpha) * pow(x,Alpha-1) * exp(-Beta*x) / Gamma_Function(Alpha);
		else
			return exp( (Alpha - 1.0) * log(x) - x - Ln_Gamma_Function(Alpha) );
	}

	//Shape & Scale
	if (Theta>0) {		
		if (Alpha<=Gamma_Function_Max_Arg())
      			return pow(x,Alpha - 1.0) * exp(-x/Theta) / (pow(Theta,Alpha) * Gamma_Function(Alpha));
   		else 
			return exp( (Alpha - 1.0) * log(x) - x - Ln_Gamma_Function(Alpha) );
	}

}


//Cumulative Distribution Function
double Gamma_Distribution::GetCDF(double x)
{
	if (error<0)
		return error;
	
	if (x <= 0.0) 
		return 0.0;

	//Shape and Rate
	if (Beta>0) {
		if (Alpha<=Gamma_Function_Max_Arg())
			return Lower_Incomplete_Gamma_Function(x*Beta,Alpha) / Gamma_Function(Alpha);
		else
			return 1.0;
	}

	//Shape & Scale
	if (Theta>0) {		
		if (Alpha<=Gamma_Function_Max_Arg())
			return Lower_Incomplete_Gamma_Function(x/Theta,Alpha) / Gamma_Function(Alpha);
		else
			return 1.0;	
	}

}

//Mean
double 	Gamma_Distribution::GetMean(void)
{
	if (error<0)
		return error;

	//Shape and Rate
	if (Beta>0) {
		return Alpha/Beta;
	}

	//Shape & Scale
	if (Theta>0) {		
		return Alpha*Theta;
	}
}


//Variance
double 	Gamma_Distribution::GetVariance(void)
{
	if (error<0)
		return error;

	//Shape and Rate
	if (Beta>0) {
		return Alpha/pow(Beta,2);
	}

	//Shape & Scale
	if (Theta>0) {		
		return Alpha*pow(Theta,2);
	}
}


//Std Deviation
double 	Gamma_Distribution::GetStdDeviation(void)
{
	if (error<0)
		return error;

	//Shape and Rate
	if (Beta>0) {
		return sqrt(Alpha/pow(Beta,2));
	}

	//Shape & Scale
	if (Theta>0) {		
		return sqrt(Alpha*pow(Theta,2));
	}
}


//Skewness
double 	Gamma_Distribution::GetSkewness(void)
{
	if (error<0)
		return error;

	//Shape and Rate
	if (Beta>0) {
		return 2.0/sqrt(Alpha);
	}

	//Shape & Scale
	if (Theta>0) {		
		return 2.0/sqrt(Alpha);
	}
}


//Kurtosis
double 	Gamma_Distribution::GetKurtosis(void)
{
	if (error<0)
		return error;

	//Shape and Rate
	if (Beta>0) {
		return 6.0/Alpha;
	}

	//Shape & Scale
	if (Theta>0) {		
		return 6.0/Alpha;
	}
}



//Return Quantile z(P) from probability P
double Gamma_Distribution::GetQuantile(double p)
{
double Vm;
double Vh = 100;
double Vl = 0;
double Pr;
int i = 0;
double Eps;

	if (p <= 0.0) {
		return Vl;
	} else if (p >= 1.0) {
        	return Vh;
	} else {        
        	do 
		{
          		i++;
          		Vm = (Vh + Vl) / 2;
            
			Pr = GetCDF(Vm);
          		Eps = abs(Pr - p);
			        
          		//New Boundary selection
          		if (Pr > p) {
				Vh = Vm;
          		} else {
				Vl = Vm;
			}
            
        	} 
		while ((Eps > CONSTANT_EpsStop) && (i < 70));
	}
            
        if (i >= 70) {
            return -9999;
        } else {
            return Vm;
    	}

}


double Gamma_Distribution::GetAlpha(void)
{
	return Alpha;
}

double Gamma_Distribution::GetBeta(void)
{
	return Beta;
}

double Gamma_Distribution::GetTheta(void)
{
	return Theta;
}

double Gamma_Distribution::GetShape(void)
{
	return Alpha;
}

double Gamma_Distribution::GetRate(void)
{
	if (Beta>0)
		return Beta;
	
	return 1/Theta;
}

double Gamma_Distribution::GetScale(void)
{
	if (Theta>0)
		return Theta;
	
	return 1/Beta;
}


// Private Methods /////////////////////////////////////////////////////////////
// Functions only available to other functions in this library


// /////////////////////////////////////////////////////////////////////////////

