#include "grid.hpp"
#include <cmath>
//------------------------------------------------------------------------------
/// Constructs a grid based on a geometry
// @param geom information about the geometry
Grid::Grid(const Geometry* geom) {
  multi_real_t offset;
  offset[0] = 0;
  offset[1] = 0;
  // Allocate pointer for the raw data with fixed size (=length) given by geometry
  _data = new real_t[geom->Size()[0]*_geom->Size()[1]];
  _geom = geom;
  _offset = offset;
}

/// Constructs a grid based on a geometry with an offset
// @param geom   information about the geometry
// @param offset distance of staggered grid point to cell's anchor point;
//               (anchor point = lower left corner)
Grid::Grid(const Geometry* geom, const multi_real_t& offset) {
  // Allocate pointer for the raw data with fixed size (=length) given by geometry
  _data = new real_t[geom->Size()[0]*_geom->Size()[1]];
  _geom = geom;
  _offset = offset;
}

/// Destructor: deletes the grid
Grid::~Grid() {delete[] _data;}

/// Initializes the grid with a given value
// @param value initial value for the whole grid
void Grid::Initialize(const real_t& value) {
  for (index_t i = 0; i < _geom->Size()[0]*_geom->Size()[1]; i++) {
    _data[i] = value;
  }
}

/// Write access to the grid cell at position [it]
// @param it Iterator instance
real_t& Grid::Cell(const Iterator& it) {
  return _data[it]; // Uses cast command via Iterator::operator
}

/// Read access to the grid cell at position [it]
// @param it Iterator instance
const real_t& Grid::Cell(const Iterator& it) const {
  return _data[it]; // Uses cast command via Iterator::operator
}

/// Interpolates the value at an arbitrary position pos with bilinear interpolation
/// for the affiliated cells using iterators
// @param pos given position in [0,1]x[0,1]
real_t Grid::Interpolate(const multi_real_t& pos) const {
  // ToDo
}

/* Calculate differential operators
*/
/// Computes the left-sided difference quotient in x-dim at position [it]
// @param it Iterator instance
real_t Grid::dx_l(const Iterator& it) const {
  return (_data[it] - _data[it.Left()])/_geom->Mesh()[0];
}

/// Computes the right-sided difference quotient in x-dim at position [it]
// @param it Iterator instance
real_t Grid::dx_r(const Iterator& it) const {
  return (_data[it.Right()] - _data[it])/_geom->Mesh()[0];
}

/// Computes the left-sided difference quotient in y-dim at position [it]
// @param it Iterator instance
real_t Grid::dy_d(const Iterator& it) const {
  return (_data[it] - _data[it.Down()])/_geom->Mesh()[1];
}

/// Computes the right-sided difference quotient in y-dim at position [it]
// @param it Iterator instance
real_t Grid::dy_t(const Iterator& it) const {
  return (_data[it.Top()] - _data[it])/_geom->Mesh()[1];
}

/// Computes the central difference quotient of 2nd order in x-dim at [it]
// @param it Iterator instance
real_t Grid::dxx(const Iterator& it) const {
  return (_data[it.Right()] - 2*_data[it] + _data[it.Left()])/_geom->Mesh()[0]/_geom->Mesh()[0];
}

/// Computes the central difference quotient of 2nd order in y-dim at [it]
// @param it Iterator instance
real_t Grid::dyy(const Iterator& it) const {
  return (_data[it.Top()] - 2*_data[it] + _data[it.Down()])/_geom->Mesh()[1]/_geom->Mesh()[1];
}

/* Donor-Cell-method for the convective terms (discretized)
*/
/// Computes u*du/dx with the donor cell method
// @param it    Iterator instance
// @param alpha weighting factor between original 2nd order central difference scheme and Donor-Cell-method
real_t Grid::DC_udu_x(const Iterator& it, const real_t& alpha) const {
  real_t meanu_r = (_data[it] + _data[it.Right()])/2;
  real_t meanu_l = (_data[it.Left()] + _data[it])/2;
  real_t diffu_r = (_data[it] - _data[it.Right()])/2;
  real_t diffu_l = (_data[it.Left()] - _data[it])/2;
  return ((meanu_r*meanu_r) - (meanu_l*meanu_l) + alpha*(fabs(meanu_r)*diffu_r - fabs(meanu_l)*diffu_l))/_geom->Mesh()[0];
}

/// Computes v*du/dy with the donor cell method
// @param it    Iterator instance
// @param alpha weighting factor between original 2nd order central difference scheme and Donor-Cell-method
// @param v     grid instance for v
real_t Grid::DC_vdu_y(const Iterator& it, const real_t& alpha, const Grid* v) const {
  real_t meanv_r = (v->Cell(it) + v->Cell(it.Right()))/2;
  real_t meanv_rd = (v->Cell(it.Down()) + v->Cell(it.Down().Right()))/2;
  real_t meanu_t = (_data[it] + _data[it.Top()])/2;
  real_t meanu_d = (_data[it.Down()] + _data[it])/2;
  real_t diffu_t = (_data[it] - _data[it.Top()])/2;
  real_t diffu_d = (_data[it.Down()] - _data[it])/2;
  return ((meanv_r*meanu_t) - (meanv_rd*meanu_d) + alpha*(fabs(meanv_r)*diffu_t - fabs(meanv_rd)*diffu_d))/_geom->Mesh()[1];
}

/// Computes u*dv/dx with the donor cell method
// @param it    Iterator instance
// @param alpha weighting factor between original 2nd order central difference scheme and Donor-Cell-method
// @param u     grid instance for u
real_t Grid::DC_udv_x(const Iterator& it, const real_t& alpha, const Grid* u) const {
  real_t meanu_t = (u->Cell(it) + u->Cell(it.Top()))/2;
  real_t meanu_tl = (u->Cell(it.Left()) + u->Cell(it.Left().Top()))/2;
  real_t meanv_r = (_data[it] + _data[it.Right()])/2;
  real_t meanv_l = (_data[it.Left()] + _data[it])/2;
  real_t diffv_r = (_data[it] - _data[it.Right()])/2;
  real_t diffv_l = (_data[it.Left()] - _data[it])/2;
  return ((meanu_t*meanv_r) - (meanu_tl*meanv_l) + alpha*(fabs(meanu_t)*diffv_r - fabs(meanu_tl)*diffv_l))/_geom->Mesh()[0];
}

/// Computes v*dv/dy with the donor cell method
// @param it    Iterator instance
// @param alpha weighting factor between original 2nd order central difference scheme and Donor-Cell-method
real_t Grid::DC_vdv_y(const Iterator& it, const real_t& alpha) const {
  real_t meanv_t = (_data[it] + _data[it.Top()])/2;
  real_t meanv_d = (_data[it.Down()] + _data[it])/2;
  real_t diffv_t = (_data[it] - _data[it.Top()])/2;
  real_t diffv_d = (_data[it.Down()] - _data[it])/2;
  return ((meanv_t*meanv_t) - (meanv_d*meanv_d) + alpha*(fabs(meanv_t)*diffv_t - fabs(meanv_d)*diffv_d))/_geom->Mesh()[1];
}

 /// Returns the maximal value of the grid
real_t Grid::Max() const {
  real_t max = _data[0];
  for (index_t i = 1; i < _geom->Size()[0]*_geom->Size()[1]; i++) {
    if (_data[i] > max)
      max = _data[i];
  }
  return max;
}

/// Returns the maximal value of the interior points
real_t Grid::InteriorMax() const {
  InteriorIterator intit(_geom);
  real_t max = 0;
  while (intit.Valid()) {
    max = fmax(max, _data[it]);
    intit.Next();
  }
  return max;
}

/// Returns the minimal value of the grid
real_t Grid::Min() const {
  real_t min = _data[0];
  for (index_t i = 1; i < _geom->Size()[0]*_geom->Size()[1]; i++) {
    if (_data[i] <= min)
      min = _data[i];
  }
  return min;
}

/// Returns the absolute maximal value
real_t Grid::AbsMax() const {
  real_t max = this->InteriorMax();
  real_t min = this->Min();
  if ((max + min) > 0)
    return max;
  else
    return min;
}

/// Returns a pointer to the raw data
real_t* Grid::Data() {return _data;}
//------------------------------------------------------------------------------