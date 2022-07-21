//#This file is part of PyTransport.

//#PyTransport is free software: you can redistribute it and/or modify
//#it under the terms of the GNU General Public License as published by
//#the Free Software Foundation, either version 3 of the License, or
//#(at your option) any later version.

//#PyTransport is distributed in the hope that it will be useful,
//#but WITHOUT ANY WARRANTY; without even the implied warranty of
//#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//#GNU General Public License for more details.

//#You should have received a copy of the GNU General Public License
//#along with PyTransport.  If not, see <http://www.gnu.org/licenses/>.

// This file contains a prototype of the potential.h file of PyTransport -- it is edited by the PyTransScripts module

#ifndef POTENTIAL_H  // Prevents the class being re-defined
#define POTENTIAL_H


#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>

using namespace std;

// #Rewrite
// Potential file rewriten at Thu Jul 21 19:22:47 2022

class potential
{
private:
	int nF; // field number
	int nP; // params number which definFs potential


public:
	// flow constructor
	potential()
	{
// #FP
nF=2;
nP=4;

//        p.resize(nP);

// pdef

    }

    //void setP(vector<double> pin){
    //    p=pin;
    //}
	//calculates V()
	double V(vector<double> f, vector<double> p)
	{
		double sum ;

// Pot
  double x0 = std::pow(f[1], 2);
  double x1 = std::pow(p[0], 2);
  double x2 = f[0] - p[1];
  double x3 = (1.0/2.0)*std::pow(x2, 2);
  double x4 = std::pow(p[2]*x3 + (1.0/6.0)*p[3]*std::pow(x2, 3) + 1.0, 2);
  double x5 = 2*x1;
  sum=3*x0*x1*x4 - x0*x5*std::pow(p[2]*x2 + p[3]*x3, 2) - x4*x5/std::pow(f[0], 2);
         return sum;
	}

	//calculates V'()
	vector<double> dV(vector<double> f, vector<double> p)
	{
		vector<double> sum(nF,0.0);

// dPot
  double x0 = std::pow(p[0], 2);
  double x1 = f[0] - p[1];
  double x2 = std::pow(x1, 2);
  double x3 = (1.0/2.0)*p[2]*x2 + (1.0/6.0)*p[3]*std::pow(x1, 3) + 1.0;
  double x4 = std::pow(x3, 2);
  double x5 = std::pow(f[1], 2);
  double x6 = 2*f[0] - 2*p[1];
  double x7 = p[3]*x2;
  double x8 = p[2]*x6 + x7;
  double x9 = 2*x0;
  double x10 = p[2]*x1 + (1.0/2.0)*x7;

 sum[0]=3*x0*x3*x5*x8 - x10*x5*x9*(2*p[2] + p[3]*x6) - x3*x8*x9/std::pow(f[0], 2) + 4*x0*x4/std::pow(f[0], 3);

 sum[1]=-4*f[1]*x0*std::pow(x10, 2) + 6*f[1]*x0*x4;

		return sum;
	}

	// calculates V''
	vector<double> dVV(vector<double> f, vector<double> p)
	{
		vector<double> sum(nF*nF,0.0);

// ddPot
  double x0 = std::pow(p[0], 2);
  double x1 = f[0] - p[1];
  double x2 = std::pow(x1, 2);
  double x3 = (1.0/2.0)*p[2]*x2 + (1.0/6.0)*p[3]*std::pow(x1, 3) + 1.0;
  double x4 = std::pow(x3, 2);
  double x5 = x0*x4;
  double x6 = std::pow(f[0], -3);
  double x7 = 2*f[0];
  double x8 = -2*p[1] + x7;
  double x9 = p[2]*x8;
  double x10 = p[3]*x2;
  double x11 = x10 + x9;
  double x12 = std::pow(f[1], 2);
  double x13 = p[3]*x8;
  double x14 = 2*p[2] + x13;
  double x15 = (1.0/2.0)*x10;
  double x16 = x15 + (1.0/2.0)*x9;
  double x17 = 2/std::pow(f[0], 2);
  double x18 = x0*x3;
  double x19 = x11*x17;
  double x20 = 2*x0*x12*x14;
  double x21 = p[2]*x1 + x15;
  double x22 = 4*x0;
  double x23 = x21*x22;
  double x24 = 1.0/p[1];
  double x25 = std::pow(x21, 2)*x22;
  double x26 = 6*f[1]*x0*x11*x3 - f[1]*x14*x23 - x24*(6*f[1]*x0*x4 - f[1]*x25);

 sum[0]=-p[3]*x12*x23 + 3*x0*x11*x12*x16 + 8*x0*x11*x3*x6 + 3*x0*x12*x14*x3 - x0*x16*x19 - x14*x17*x18 - x20*(p[2] + (1.0/2.0)*x13) - 12*x5/std::pow(f[0], 4);

 sum[2]=x26;

 sum[1]=x26;

 sum[3]=x24*(3*x0*x11*x12*x3 + 4*x0*x4*x6 - x18*x19 - x20*x21)*std::exp(x24*x7) - x25 + 6*x5;

        return sum;
	}

	// calculates V'''
	vector<double> dVVV(vector<double> f, vector<double> p)
	{
        vector<double> sum(nF*nF*nF,0.0);
// dddPot
  double x0 = std::pow(p[0], 2);
  double x1 = f[0] - p[1];
  double x2 = std::pow(x1, 2);
  double x3 = (1.0/2.0)*p[2]*x2 + (1.0/6.0)*p[3]*std::pow(x1, 3) + 1.0;
  double x4 = std::pow(x3, 2);
  double x5 = std::pow(f[0], -3);
  double x6 = 2*f[0];
  double x7 = -2*p[1] + x6;
  double x8 = p[3]*x7;
  double x9 = 2*p[2] + x8;
  double x10 = p[2]*x7;
  double x11 = p[3]*x2;
  double x12 = x10 + x11;
  double x13 = (1.0/2.0)*x11;
  double x14 = (1.0/2.0)*x10 + x13;
  double x15 = std::pow(f[1], 2);
  double x16 = p[2] + (1.0/2.0)*x8;
  double x17 = p[3]*x0;
  double x18 = 2*x9;
  double x19 = x0*x12;
  double x20 = std::pow(f[0], -2);
  double x21 = 2*x20;
  double x22 = 4*x0;
  double x23 = p[3]*x22;
  double x24 = x20*x3;
  double x25 = x22*x9;
  double x26 = x14*x20;
  double x27 = x15*x16;
  double x28 = 8*x17;
  double x29 = std::pow(f[0], -4);
  double x30 = x19*x3;
  double x31 = 1.0/p[1];
  double x32 = p[2]*x1 + x13;
  double x33 = x25*x32;
  double x34 = -6*f[1]*x0*x12*x3 + f[1]*x33;
  double x35 = -x31*x34;
  double x36 = x22*std::pow(x32, 2);
  double x37 = x31*(-x31*(6*f[1]*x0*x4 - f[1]*x36) - x34);
  double x38 = -6*f[1]*x0*x12*x14 - 6*f[1]*x0*x3*x9 + f[1]*x16*x25 + f[1]*x28*x32;
  double x39 = -x35 - x37 - x38;
  double x40 = x0*x4;
  double x41 = x15*x32;
  double x42 = x0*x18;
  double x43 = 3*x0*x12*x15*x3 + 4*x0*x4*x5 - x21*x30 - x41*x42;
  double x44 = std::exp(x31*x6);
  double x45 = x31*x44;
  double x46 = x31*(-x36 + 6*x40 + x43*x45);
  double x47 = 3*x0*x12*x14*x15 + 8*x0*x12*x3*x5 + 3*x0*x15*x3*x9 - 2*x19*x26 - x23*x41 - x24*x42 - x27*x42 - 12*x29*x40;
  double x48 = 2*x37;
  double x49 = 6*x0*x12*x3 + x31*x44*x47 - x31*(6*x0*x4 - x36) - x33 - x46;

 sum[0]=6*p[3]*x0*x15*x3 + 12*x0*x12*x14*x5 + 3*x0*x12*x15*x16 + 6*x0*x14*x15*x9 + 12*x0*x3*x5*x9 - x15*x17*x18 - x16*x19*x21 - x23*x24 - x25*x26 - x27*x28 - 36*x29*x30 + 48*x0*x4/std::pow(f[0], 5);

 sum[4]=-x38 - x48;

 sum[2]=x39;

 sum[6]=x49;

 sum[1]=x39;

 sum[5]=x49;

 sum[3]=6*x30 - x33 + x45*x47 - 2*x46 + 2*x43*x44/std::pow(p[1], 2);

 sum[7]=x35*x44 + x44*x48;

        return sum;
	}

    int getnF()
    {
        return nF;
    }

    int getnP()
    {
        return nP;
    }

};
#endif
