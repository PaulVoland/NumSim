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

#include "typedef.hpp"
//------------------------------------------------------------------------------
#ifndef __PARAMETER_HPP
#define __PARAMETER_HPP
//------------------------------------------------------------------------------
class Parameter {
public:
  /// Constructs a new Parameter set with default values
  // Driven Cavity parameters; see exercise sheet 1
  Parameter();

  /// Loads the parameter values from a file
  void Load(const char *file);

  /// Getter functions for all parameters
  const real_t  &Re()      const;
  const real_t  &Omega()   const;
  const real_t  &Alpha()   const;
  const real_t  &Dt()      const;
  const real_t  &Tend()    const;
  const real_t  &Eps()     const;
  const real_t  &Tau()     const;
  const index_t &IterMax() const;
  const real_t  &Gx()      const;
  const real_t  &Gy()      const;
  const real_t  &Pr()      const;
  const real_t  &Beta()    const;
  const real_t  &u_D()     const;
  const real_t  &v_D()     const;
  const real_t  &p_D()     const;
  const real_t  &T_H()     const;
  const real_t  &T_C()     const;

private:
  real_t  _re;
  real_t  _omega;
  real_t  _alpha;
  real_t  _dt;
  real_t  _tend;
  real_t  _eps;
  real_t  _tau;
  index_t _itermax;
  real_t  _gx;
  real_t  _gy;
  real_t  _pr;
  real_t  _beta;
  real_t  _u_Dir;
  real_t  _v_Dir;
  real_t  _p_Dir;
  real_t  _T_h;
  real_t  _T_c;
};
//------------------------------------------------------------------------------
#endif // __PARAMETER_HPP
