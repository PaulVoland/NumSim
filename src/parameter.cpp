#include "parameter.hpp"
#include <iostream>
#include <cstdio>
#include <cstring>
#include <cstdlib>

using namespace std;
//------------------------------------------------------------------------------
/// Constructs a new parameter set with default lid driven cavity values
Parameter::Parameter() {
	_re      = 1000.0;
	_omega   = 1.7;
	_alpha   = 0.9;
	_dt      = 0.2;
	_tend    = 16.4;
	_eps     = 0.001;
	_tau     = 0.5;
	_itermax = 100;
	// Message of success
	cout << "Loaded default parameters." << endl;
}

/// Loads the parameter values from a specified file
// @param file filename including path information
void Parameter::Load(const char* file) {
	FILE* handle = fopen(file,"r");
	double inval;
	char name[20];
		while (!feof(handle)) {
		if (!fscanf(handle, "%s = %lf\n", name, &inval)) continue;
		if (strcmp(name,"re") == 0) _re = inval;
		else if (strcmp(name,"omg") == 0) _omega = inval;
		else if (strcmp(name,"alpha") == 0) _alpha = inval;
		else if (strcmp(name,"dt") == 0) _dt = inval;
		else if (strcmp(name,"tend") == 0) _tend = inval;
		else if (strcmp(name,"iter") == 0) _itermax = inval;
		else if (strcmp(name,"eps") == 0) _eps = inval;
		else if (strcmp(name,"tau") == 0) _tau = inval;
		else printf("Unknown parameter %s\n",name);
	}
	fclose(handle);
	cout << "Loaded parameters from file " << file << "." << endl;
}

/// Prints the parameter set
void Parameter::PrintVariables() {
	cout << "Showing used set of parameters..." 			  << endl;
	cout << "Re = " 							<< _re 	      << endl;
	cout << "omega = " 							<< _omega	  << endl;
	cout << "alpha = " 							<< _alpha	  << endl;
	cout << "dt = " 							<< _dt 		  << endl;
	cout << "tend = " 							<< _tend	  << endl;
	cout << "eps = "						    << _eps   	  << endl;
	cout << "tau = "						    << _tau 	  << endl;
	cout << "itermax = "						<< _itermax   << endl;
}

/// Getter functions for all parameters
const real_t&  Parameter::Re()      const {return _re;}
const real_t&  Parameter::Omega()   const {return _omega;}
const real_t&  Parameter::Alpha()   const {return _alpha;}
const real_t&  Parameter::Dt()      const {return _dt;}
const real_t&  Parameter::Tend()    const {return _tend;}
const real_t&  Parameter::Eps()     const {return _eps;}
const real_t&  Parameter::Tau()     const {return _tau;}
const index_t& Parameter::IterMax() const {return _itermax;}
//------------------------------------------------------------------------------