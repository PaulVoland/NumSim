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
#include "geometry.hpp"
#include "iterator.hpp"
//------------------------------------------------------------------------------
#ifndef __GRID_HPP
#define __GRID_HPP
//------------------------------------------------------------------------------
class Grid {
public:
  /// Constructs a grid based on a geometry
  Grid(const Geometry *geom);

  /// Constructs a grid based on a geometry with an offset
  // @param geom   Geometry information
  // @param offset distance of staggered grid point to cell's anchor point;
  //               (anchor point = lower left corner)
  Grid(const Geometry *geom, const multi_real_t &offset);

  /// Deletes the grid
  ~Grid();

  ///	Initializes the grid with a value
  void Initialize(const real_t &value);

  /// Write access to the grid cell at position [it]
  real_t &Cell(const Iterator &it);
  /// Read access to the grid cell at position [it]
  const real_t &Cell(const Iterator &it) const;

  /// Interpolate the value at an arbitrary position
  real_t Interpolate(const multi_real_t &pos) const;

  /// Computes the left-sided difference quotient in x-dim at [it]
  real_t dx_l(const Iterator &it) const;
  /// Computes the right-sided difference quotient in x-dim at [it]
  real_t dx_r(const Iterator &it) const;
  /// Computes the left-sided difference quotient in y-dim at [it]
  real_t dy_d(const Iterator &it) const;
  /// Computes the right-sided difference quotient in y-dim at [it]
  real_t dy_t(const Iterator &it) const;
  /// Computes the central difference quotient of 2nd order in x-dim at [it]
  real_t dxx(const Iterator &it) const;
  /// Computes the central difference quotient of 2nd order in y-dim at [it]
  real_t dyy(const Iterator &it) const;

  /// Computes u*du/dx with the donor cell method
  real_t DC_udu_x(const Iterator &it, const real_t &alpha) const;
  /// Computes v*du/dy with the donor cell method
  real_t DC_vdu_y(const Iterator &it, const real_t &alpha, const Grid *v) const;
  /// Computes u*dv/dx with the donor cell method
  real_t DC_udv_x(const Iterator &it, const real_t &alpha, const Grid *u) const;
  /// Computes v*dv/dy with the donor cell method
  real_t DC_vdv_y(const Iterator &it, const real_t &alpha) const;
  /// Computes u*dT/dx with the donor cell method
  real_t DC_udT_x(const Iterator &it, const real_t &alpha, const Grid *u) const;
  /// Computes v*dT/dy with the donor cell method
  real_t DC_vdT_y(const Iterator &it, const real_t &alpha, const Grid *v) const;


  /// Returns the maximal value of the grid
  real_t Max() const;
  /// Returns the maximal value of the interior points
  real_t InteriorMax() const;
  /// Returns the minimal value of the grid
  real_t Min() const;
  /// Returns the minimal value of the interior points
  real_t InteriorMin() const;
  /// Returns the absolute maximal value
  real_t AbsMax() const;

  /// Returns a pointer to the raw data
  real_t *Data();

  /// Getter for the grid offset
  const multi_real_t &getOffset() const;

  /// Returns a pointer to the Geometry
  const Geometry *getGeometry() const;

private:
  real_t *_data;
  multi_real_t _offset;
  const Geometry *_geom;
};
//------------------------------------------------------------------------------
#endif // __GRID_HPP
