#include "grid.hpp"
#include <cmath>
#include <algorithm>

using namespace std;
//------------------------------------------------------------------------------
/// Constructs a grid based on a geometry
// @param geom information about the geometry
Grid::Grid(const Geometry* geom) : _geom(geom) {
  _offset[0] = 0.0;
  _offset[1] = 0.0;
  // Allocate pointer for the raw data with fixed size (=length) given by geometry
  // (already included halo zone in geometry.cpp)
  _data = new real_t[geom->Size()[0]*geom->Size()[1]];
}
//------------------------------------------------------------------------------
/// Constructs a grid based on a geometry with an offset
// @param geom   information about the geometry
// @param offset distance of staggered grid point to cell's anchor point;
//               (anchor point = lower left corner)
Grid::Grid(const Geometry* geom, const multi_real_t& offset) : _geom(geom) {
    _offset = offset;
  // Allocate pointer for the raw data with fixed size (=length) given by geometry
  // (already included halo zone in geometry.cpp)
  _data = new real_t[geom->Size()[0]*geom->Size()[1]];
}
//------------------------------------------------------------------------------
/// Destructor: deletes the grid
Grid::~Grid() {delete[] _data;}
//------------------------------------------------------------------------------
/// Initializes the grid with a given value
// @param value initial value for the whole grid
void Grid::Initialize(const real_t& value) {
  fill_n(_data, geom->Size()[0]*geom->Size()[1], value);
}
//------------------------------------------------------------------------------
/// Write access to the grid cell at position [it]
// @param it Iterator instance
real_t& Grid::Cell(const Iterator& it) {
  return _data[it]; // Uses cast command via Iterator::operator
}
//------------------------------------------------------------------------------
/// Read access to the grid cell at position [it]
// @param it Iterator instance
const real_t& Grid::Cell(const Iterator& it) const {
  return _data[it]; // Uses cast command via Iterator::operator
}
//------------------------------------------------------------------------------
/// Interpolates the value at an arbitrary position pos with bilinear interpolation
/// for the affiliated cells using iterators
// @param pos given position in [0,1]x[0,1]
real_t Grid::Interpolate(const multi_real_t& pos) const {
  real_t x          = pos[0] - _offset[0];
  real_t y          = pos[1] - _offset[1];
  multi_real_t h    = _geom->Mesh();
  index_t _increm_y = _geom->Size()[0];
  // Instantiate indices and distances
  index_t i, j;
  real_t dx1, dx2, dy1, dy2;

  // find inner cell index for anchor cell in format
  // h(0)*[i,i+1) x h(1)*[j,j+1) (i,j = 0,...,_geom->Size()[0,1])-1) (inner numbering)
  // if x,y >= 0
  if (x < 0) {
    i = 0; // is in outer index format
    dx1 = x + h[0];
    dx2 = -x;
  } else {
    i = (index_t) (x/h[0]); // is in inner index format
    // calculate distances to anchor point
    dx1 = x - i*h[0];
    dx2 = (i + 1)*h[0] - x;
    i++; // convert to outer index format
  }
  if (y < 0) {
    j = 0; // is in outer index format
    dy1 = y + h[1];
    dy2 = -y;
  } else {
    j = (index_t) (y/h[1]); // is in inner index format
    // calculate distances to anchor point
    dy1 = y - j*h[1];
    dy2 = (j + 1)*h[1] - y;
    j++; // convert to outer index format
  }

  Iterator it = Iterator(_geom, i + j*_increm_y);
  // bilinear interpolation
  return 1.0/(h[0]*h[1])*(
    _data[it]               *dx2  *dy2 + 
    _data[it.Right()]       *dx1  *dy2 + 
    _data[it.Top()]         *dx2  *dy1 + 
    _data[it.Right().Top()] *dx1  *dy1 );
}
//------------------------------------------------------------------------------
/* Calculate differential operators
*/
/// Computes the left-sided difference quotient in x-dim at position [it]
// @param it Iterator instance
real_t Grid::dx_l(const Iterator& it) const {
  return (_data[it] - _data[it.Left()])/_geom->Mesh()[0];
}
//------------------------------------------------------------------------------
/// Computes the right-sided difference quotient in x-dim at position [it]
// @param it Iterator instance
real_t Grid::dx_r(const Iterator& it) const {
  return (_data[it.Right()] - _data[it])/_geom->Mesh()[0];
}
//------------------------------------------------------------------------------
/// Computes the left-sided difference quotient in y-dim at position [it]
// @param it Iterator instance
real_t Grid::dy_d(const Iterator& it) const {
  return (_data[it] - _data[it.Down()])/_geom->Mesh()[1];
}
//------------------------------------------------------------------------------
/// Computes the right-sided difference quotient in y-dim at position [it]
// @param it Iterator instance
real_t Grid::dy_t(const Iterator& it) const {
  return (_data[it.Top()] - _data[it])/_geom->Mesh()[1];
}
//------------------------------------------------------------------------------
/// Computes the central difference quotient of 2nd order in x-dim at [it]
// @param it Iterator instance
real_t Grid::dxx(const Iterator& it) const {
  return (_data[it.Right()] - 2*_data[it] + _data[it.Left()])/_geom->Mesh()[0]/_geom->Mesh()[0];
}
//------------------------------------------------------------------------------
/// Computes the central difference quotient of 2nd order in y-dim at [it]
// @param it Iterator instance
real_t Grid::dyy(const Iterator& it) const {
  return (_data[it.Top()] - 2*_data[it] + _data[it.Down()])/_geom->Mesh()[1]/_geom->Mesh()[1];
}
//------------------------------------------------------------------------------
/* Donor-Cell-method for the convective terms (discretized)
*/
/// Computes u*du/dx with the donor cell method
// @param it    Iterator instance
// @param alpha weighting factor between original 2nd order central difference scheme and Donor-Cell-method
real_t Grid::DC_udu_x(const Iterator& it, const real_t& alpha) const {
  real_t meanu_r = _data[it] + _data[it.Right()];
  real_t meanu_l = _data[it.Left()] + _data[it];
  real_t diffu_r = _data[it] - _data[it.Right()];
  real_t diffu_l = _data[it.Left()] - _data[it];
  return ((meanu_r*meanu_r) - (meanu_l*meanu_l) + alpha*(fabs(meanu_r)*diffu_r - fabs(meanu_l)*diffu_l))/(4.0*_geom->Mesh()[0]);
}
//------------------------------------------------------------------------------
/// Computes v*du/dy with the donor cell method
// @param it    Iterator instance
// @param alpha weighting factor between original 2nd order central difference scheme and Donor-Cell-method
// @param v     grid instance for v
real_t Grid::DC_vdu_y(const Iterator& it, const real_t& alpha, const Grid* v) const {
  real_t meanv_r  = v->Cell(it) + v->Cell(it.Right());
  real_t meanv_rd = v->Cell(it.Down()) + v->Cell(it.Down().Right());
  real_t meanu_t  = _data[it] + _data[it.Top()];
  real_t meanu_d  = _data[it.Down()] + _data[it];
  real_t diffu_t  = _data[it] - _data[it.Top()];
  real_t diffu_d  = _data[it.Down()] - _data[it];
  return ((meanv_r*meanu_t) - (meanv_rd*meanu_d) + alpha*(fabs(meanv_r)*diffu_t - fabs(meanv_rd)*diffu_d))/(4.0*_geom->Mesh()[1]);
}
//------------------------------------------------------------------------------
/// Computes u*dv/dx with the donor cell method
// @param it    Iterator instance
// @param alpha weighting factor between original 2nd order central difference scheme and Donor-Cell-method
// @param u     grid instance for u
real_t Grid::DC_udv_x(const Iterator& it, const real_t& alpha, const Grid* u) const {
  real_t meanu_t  = u->Cell(it) + u->Cell(it.Top());
  real_t meanu_tl = u->Cell(it.Left()) + u->Cell(it.Left().Top());
  real_t meanv_r  = _data[it] + _data[it.Right()];
  real_t meanv_l  = _data[it.Left()] + _data[it];
  real_t diffv_r  = _data[it] - _data[it.Right()];
  real_t diffv_l  = _data[it.Left()] - _data[it];
  return ((meanu_t*meanv_r) - (meanu_tl*meanv_l) + alpha*(fabs(meanu_t)*diffv_r - fabs(meanu_tl)*diffv_l))/(4.0*_geom->Mesh()[0]);
}
//------------------------------------------------------------------------------
/// Computes v*dv/dy with the donor cell method
// @param it    Iterator instance
// @param alpha weighting factor between original 2nd order central difference scheme and Donor-Cell-method
real_t Grid::DC_vdv_y(const Iterator& it, const real_t& alpha) const {
  real_t meanv_t = _data[it] + _data[it.Top()];
  real_t meanv_d = _data[it.Down()] + _data[it];
  real_t diffv_t = _data[it] - _data[it.Top()];
  real_t diffv_d = _data[it.Down()] - _data[it];
  return ((meanv_t*meanv_t) - (meanv_d*meanv_d) + alpha*(fabs(meanv_t)*diffv_t - fabs(meanv_d)*diffv_d))/(4.0*_geom->Mesh()[1]);
}
//------------------------------------------------------------------------------
/// Returns the maximal value of the grid
real_t Grid::Max() const {
  real_t max = _data[0];
  for (index_t i = 1; i < _geom->Size()[0]*_geom->Size()[1]; i++) {
    if (_data[i] > max)
      max = _data[i];
  }
  return max;
}
//------------------------------------------------------------------------------
/// Returns the maximal value of the interior points
real_t Grid::InteriorMax() const {
  InteriorIterator intit(_geom);
  real_t max = 0;
  while (intit.Valid()) {
    max = fmax(max, _data[intit]);
    intit.Next();
  }
  return max;
}
//------------------------------------------------------------------------------
/// Returns the minimal value of the grid
real_t Grid::Min() const {
  real_t min = _data[0];
  for (index_t i = 1; i < _geom->Size()[0]*_geom->Size()[1]; i++) {
    if (_data[i] < min)
      min = _data[i];
  }
  return min;
}
//------------------------------------------------------------------------------
/// Returns the minimal value of the interior points
real_t Grid::InteriorMin() const {
  InteriorIterator intit(_geom);
  real_t min = 0;
  while (intit.Valid()) {
    min = fmin(min, _data[intit]);
    intit.Next();
  }
  return min;
}
//------------------------------------------------------------------------------
/// Returns the absolute maximal value
real_t Grid::AbsMax() const {
  real_t max = this->Max();
  real_t min = this->Min();
  if ((max + min) > 0)
    return max;
  else
    return -min;
}
//------------------------------------------------------------------------------
/// Returns a pointer to the raw data
real_t* Grid::Data() {return _data;}
//------------------------------------------------------------------------------
/// Getter for the grid offset
const multi_real_t& Grid::getOffset() const {return _offset;}
//------------------------------------------------------------------------------
/// Returns a pointer to the Geometry
const Geometry* Grid::getGeometry() const {return _geom;}
//------------------------------------------------------------------------------