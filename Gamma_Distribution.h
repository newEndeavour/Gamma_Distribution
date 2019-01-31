/*
  File:         Gamma_Distribution.h
  Version:      0.0.1
  Date:         30-Jan-2019
  Revision:     30-Jan-2019
  Author:       Jerome Drouin (jerome.p.drouin@gmail.com)

  Editions:
  - 0.0.1	: First version
  - 0.0.2	: -

  References:
  https://github.com/newEndeavour/Gamma_Distribution
  https://en.wikipedia.org/wiki/Gamma_distribution
  A special thank you to:
  http://www.mymathlib.com/functions/probability/gamma_distribution.html
  
  Gamma_Distribution.h - Library for 'duino

  Gamma_Distribution implements the Gamma Distribution with shape Alpha and 
  Rate Beta or Scale Theta. 

  The object enables two (mutually exclusive) characterisations:
	- Shape (Alpha) & Rate (Beta) 	:  Gamma_Distribution(Shape>0,Rate>0,0)
					 = Gamma_Distribution(Alpha>0,Beta>0,0) 
	- Shape (Alpha) & Scale (Theta)	:  Gamma_Distribution(Shape>0,0,Scale)
					 = Gamma_Distribution(Alpha>0,0,Theta>0) 
  An attempt to construct an object with positive Alpha, Beta and Theta creates an 
  indetermination and returns an error flag.				
  See also: https://en.wikipedia.org/wiki/Gamma_distribution#Characterization_using_shape_α_and_rate_β


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


// ensure this library description is only included once
#ifndef Gamma_Distribution_h
#define Gamma_Distribution_h


#if ARDUINO >= 100
#include "Arduino.h"
#else
#include "WProgram.h"
#endif


//Gamma Constant Parameters
#define CONSTANT_EpsStop 0.0000001
    
//Number Pi
#define CONSTANT_Pi 3.14159265358979


// library interface description
class Gamma_Distribution
{
  // user-accessible "public" interface
  public:
  // methods
	Gamma_Distribution(double _Alpha, double _Beta, double _Theta);
	
	double 	GetPDF(double x);
	double 	GetCDF(double x);
	double 	GetQuantile(double p);

	double 	GetShape(void);
	double 	GetScale(void);
	double 	GetRate(void);

	double 	GetAlpha(void);
	double 	GetBeta(void);
	double 	GetTheta(void);

  // library-accessible "private" interface
  private:
  // variables
	int 	error;

	double	Alpha;			// 
	double	Beta;			// 
	double	Theta;			// 
	
  // methods
};

#endif
