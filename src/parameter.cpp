/*
 * Copyright (C) 2015   Malte Brunn
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

//------------------------------------------------------------------------------
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
  _dt      = 0.0;
  _tend    = 50.0;
  _eps     = 0.001;
  _tau     = 0.5;
  _itermax = 100;
  _gx      = 0.0;
  _gy      = 0.0;
  _pr      = 0.0;
  _beta    = 0.0;
  _u_Dir   = 1.0;
  _v_Dir   = 0.0;
  _p_Dir   = 0.0;
  _T_h     = 0.0;
  _T_c     = 0.0;
  _k       = 1.0;

  // Message of success
  cout << "Loaded default parameters." << endl;
}
//------------------------------------------------------------------------------
/// Loads the parameter values from a specified file
// @param file filename including path information
void Parameter::Load(const char* file) {
  FILE *handle = fopen(file, "r");
  double inval;
  char name[200000];
  while (!feof(handle)) {
    if (!fscanf(handle, "%s = %lf\n", name, &inval))
      continue;
    if      (strcmp(name, "re")       == 0)
      _re      = inval;
    else if (strcmp(name, "omg")      == 0)
      _omega   = inval;
    else if (strcmp(name, "alpha")    == 0)
      _alpha   = inval;
    else if (strcmp(name, "dt")       == 0)
      _dt      = inval;
    else if (strcmp(name, "tend")     == 0)
      _tend    = inval;
    else if (strcmp(name, "eps")      == 0)
      _eps     = inval;
    else if (strcmp(name, "tau")      == 0)
      _tau     = inval;
    else if (strcmp(name, "iter")     == 0)
      _itermax = (uint32_t) inval;
    else if (strcmp(name, "gx")       == 0)
      _gx      = inval;
    else if (strcmp(name, "gy")       == 0)
      _gy      = inval;
    else if (strcmp(name, "pr")       == 0)
      _pr      = inval;
    else if (strcmp(name, "beta")     == 0)
      _beta    = inval;
    else if (strcmp(name, "u_D")      == 0)
      _u_Dir   = inval;
    else if (strcmp(name, "v_D")      == 0)
      _v_Dir   = inval;
    else if (strcmp(name, "p_D")      == 0)
      _p_Dir   = inval;
    else if (strcmp(name, "T_h")      == 0)
      _T_h     = inval;
    else if (strcmp(name, "T_c")      == 0)
      _T_c     = inval;
    else if (strcmp(name, "k")        == 0)
      _k       = inval;
    else
      printf("Unknown parameter %s\n", name);
  }
  fclose(handle);
  cout << "Loaded parameters from file " << file << "." << endl;
}
//------------------------------------------------------------------------------
/* Getter functions for all parameters
*/
const real_t&  Parameter::Re()      const {return _re;}
//------------------------------------------------------------------------------
const real_t&  Parameter::Omega()   const {return _omega;}
//------------------------------------------------------------------------------
const real_t&  Parameter::Alpha()   const {return _alpha;}
//------------------------------------------------------------------------------
const real_t&  Parameter::Dt()      const {return _dt;}
//------------------------------------------------------------------------------
const real_t&  Parameter::Tend()    const {return _tend;}
//------------------------------------------------------------------------------
const real_t&  Parameter::Eps()     const {return _eps;}
//------------------------------------------------------------------------------
const real_t&  Parameter::Tau()     const {return _tau;}
//------------------------------------------------------------------------------
const index_t& Parameter::IterMax() const {return _itermax;}
//------------------------------------------------------------------------------
const real_t&  Parameter::Gx()      const {return _gx;}
//------------------------------------------------------------------------------
const real_t&  Parameter::Gy()      const {return _gy;}
//------------------------------------------------------------------------------
const real_t&  Parameter::Pr()      const {return _pr;}
//------------------------------------------------------------------------------
const real_t&  Parameter::Beta()    const {return _beta;}
//------------------------------------------------------------------------------
const real_t&  Parameter::u_D()     const {return _u_Dir;}
//------------------------------------------------------------------------------
const real_t&  Parameter::v_D()     const {return _v_Dir;}
//------------------------------------------------------------------------------
const real_t&  Parameter::p_D()     const {return _p_Dir;}
//------------------------------------------------------------------------------
const real_t&  Parameter::T_H()     const {return _T_h;}
//------------------------------------------------------------------------------
const real_t&  Parameter::T_C()     const {return _T_c;}
//------------------------------------------------------------------------------
const real_t&  Parameter::k()     const {return _k;}
//------------------------------------------------------------------------------
