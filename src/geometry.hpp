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
#include "iterator.hpp"
//------------------------------------------------------------------------------
#ifndef __GEOMETRY_HPP
#define __GEOMETRY_HPP
//------------------------------------------------------------------------------
/// Typedef for cell types
typedef enum {
  typeFluid,   // Standard fluid cell
  typeSolid,   // Simple wall, no slip, isolated against heat transfer
  typeIn,      // Simple inflow (forced velocity)
  typeInH,     // Horizontal inflow (parabolic)
  typeInV,     // Vertical inflow (parabolic)
  typeSlipH,   // Horizontal slip boundary
  typeSlipV,   // Vertical slip boundary
  typeOut,     // Outflow
  typeTDir_h,  // Dirichlet value for higher temperature (u,v,p treated as no slip)
  typeTDir_c   // Dirichlet value for lower temperature (u,v,p treated as no slip)

  typeCoupling // NEW  _interface to a coupled solid field


} CellType_t;
/// Typedef for cell boundary type (which boundary cells are fluid)
//       |   NO   |
//   ____|________|____
//       |        |
//   WE  |  CELL  |  EA
//   ____|________|____
//       |        |
//       |   SO   |
// Cells on the diagonals can be ignored
// Combinations not listed are invalid
typedef enum {
  cellNone = 0, // Cell is not surrounded by fluid, or is fluid cell
  cellN = 1,    // Cell N is fluid
  cellW = 2,    // Cell W is fluid
  cellNW = 3,   // Cells N & W are fluid
  cellS = 4,    // Cell S is fluid
  cellSW = 6,   // Cells S & W are fluid
  cellE = 8,    // Cell E is fluid
  cellNE = 9,   // Cells N & E are fluid
  cellSE = 12   // Cells S & E are fluid
} CellBoundary_t;
/// Typedef for cells
typedef struct {
  CellType_t type;      // Cell type
  CellBoundary_t fluid; // Surrounding fluid cells
  real_t factor;        // Scale factor for inital values
} Cell_t;
//------------------------------------------------------------------------------
class Geometry {
public:
  /// Constructs a default geometry:
  // driven cavity with 128 x 128 grid, no-slip boundary conditions
  // as shown below
  //
  //       u=1, v=0
  //    -------------
  //    |           |
  // u=0|           |u=0
  // v=0|           |v=0
  //    |           |
  //    |           |
  //    -------------
  //       u=0, v=0
  Geometry();
  /// Destructor
  ~Geometry();

  /// Loads a geometry from a file
  void Load(const char *file);

  /// Returns the total number of cells in each dimension
  const multi_index_t &TotalSize() const;
  /// Returns the total length of the domain
  const multi_real_t &TotalLength() const;
  /// Returns the overall meshwidth
  const multi_real_t &Mesh() const;
  /// Read access to the cell type field at position [it]
  const Cell_t &Cell(const Iterator &it) const;
  /// Returns the prescribed velocity values
  const multi_real_t &Velocity() const;
  /// Returns the prescribed pressure value
  const real_t &Pressure() const;
  /// Returns the prescribed temperature value
  const real_t &Temperature() const;
  ///
  const index_t &Num_Coupling() const;
  ///
  //##################################################
  const multi_real_t &getCouplingCellsNumbs() const;
  const multi_real_t &getHeatFluxs() const;
  void setHeatFluxs(const multi_real_t &HeatFlux) const;
  //##################################################

  /// Updates the velocity field u (parameter from .param used)
  void Update_U(Grid *u, const real_t &u_Dir) const;
  /// Updates the velocity field v (parameter from .param used)
  void Update_V(Grid *v, const real_t &v_Dir) const;
  /// Updates the pressure field p (parameter from .param used)
  void Update_P(Grid *p, const real_t &p_Dir) const;
  /// Updates the temperature field T (parameters from .param used)
  void Update_T(Grid *T, const real_t &T_h, const real_t &T_c) const;

private:
  Cell_t *_cell;
  index_t _num_coupling;

  multi_index_t _size;
  multi_real_t  _length;
  multi_real_t  _h;

  multi_real_t  _velocity;
  real_t        _pressure;
  real_t        _temperature;

//#####################################
  multi_real_t  _CouplingCellsNumbs // Dim N
  multi_real_t  _HeatFluxs // Dim N
//#####################################

  void UpdateCellDirichlet_U(Grid *u, const real_t &value, const Iterator &it) const;
  void UpdateCellDirichlet_V(Grid *v, const real_t &value, const Iterator &it) const;
  void UpdateCellNeumann(Grid *grid, const Iterator &it) const;
  void UpdateCellNeumann_P(Grid *grid, const Iterator &it) const;
  void UpdateCellDirichlet_T(Grid *T, const real_t &value, const Iterator &it) const;
};
//------------------------------------------------------------------------------
#endif // __GEOMETRY_HPP
