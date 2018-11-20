#include "geometry.hpp"
#include "iterator.hpp"
#include "grid.hpp"
#include "communicator.hpp"
#include <stdio.h>
#include <string.h>
//------------------------------------------------------------------------------
/// Constructs a default geometry
Geometry::Geometry() : _comm(NULL) {
  _length[0] = 1.0;
  _length[1] = 1.0;
  _size[0] = 128;
  _size[1] = 128;
  _h[0] = _length[0]/_size[0];
  _h[1] = _length[1]/_size[1];
  _pressure = 0.0;
  _velocity[0] = 1.0;
  _velocity[1] = 0.0;
  // Create boundary halo
  _size[0] += 2;
  _size[1] += 2;
  // Block is identical to whole domain
  _bsize = _size;
  _blength = _length;
}
//------------------------------------------------------------------------------
/// Constructs a default geometry with partition set up using the communicator object
Geometry::Geometry(const Communicator* comm) : _comm(comm) {
  _length[0] = 1.0;
  _length[1] = 1.0;
  _size[0] = 128;
  _size[1] = 128;
  _h[0] = _length[0]/_size[0];
  _h[1] = _length[1]/_size[1];
  _pressure = 0.0;
  _velocity[0] = 1.0;
  _velocity[1] = 0.0;
  // Make blocks for each communicator object including halo
  _bsize[0] = _size[0]/_comm->ThreadDim()[0] + 2;
  _bsize[1] = _size[1]/_comm->ThreadDim()[1] + 2;
  // The last communicator objects right/top take the rest of the physical geometry
  if (_comm->ThreadIdx()[0] == _comm->ThreadDim()[0] - 1)
    _bsize[0] += _size[0] % _comm->ThreadDim()[0];
  if (_comm->ThreadIdx()[1] == _comm->ThreadDim()[1] - 1)
    _bsize[1] += _size[1] % _comm->ThreadDim()[1];
  // Compute block length using block size
  _blength[0] = _h[0]*(_bsize[0] - 2);
  _blength[1] = _h[1]*(_bsize[1] - 2);
  // Create global boundary halo
  _size[0] += 2;
  _size[1] += 2;
}
//------------------------------------------------------------------------------
/// Loads a geometry from a file
void Geometry::Load(const char* file) {
  FILE *handle = fopen(file, "r");
  double inval[2];
  char name[20];
  while (!feof(handle)) {
    if (!fscanf(handle, "%s =", name))
      continue;
    if (strcmp(name, "size") == 0) {
      if (fscanf(handle, " %lf %lf\n", &inval[0], &inval[1])) {
        _size[0] = inval[0];
        _size[1] = inval[1];
      }
      continue;
    }
    if (strcmp(name, "length") == 0) {
      if (fscanf(handle, " %lf %lf\n", &inval[0], &inval[1])) {
        _length[0] = inval[0];
        _length[1] = inval[1];
      }
      continue;
    }
    if (strcmp(name, "velocity") == 0) {
      if (fscanf(handle, " %lf %lf\n", &inval[0], &inval[1])) {
        _velocity[0] = inval[0];
        _velocity[1] = inval[1];
      }
      continue;
    }
    if (strcmp(name, "pressure") == 0) {
      if (fscanf(handle, " %lf\n", &inval[0]))
        _pressure = inval[0];
      continue;
    }
  }
  fclose(handle);
  // Set width values of the grid
  _h[0] = _length[0]/_size[0];
  _h[1] = _length[1]/_size[1];
  // If !_comm --> default, see first constructor
  _blength = _length;
  // If communicator objects already exist
  if (_comm) {
    // Make blocks for each communicator object including halo
    _bsize[0] = _size[0]/_comm->ThreadDim()[0] + 2;
    _bsize[1] = _size[1]/_comm->ThreadDim()[1] + 2;
    // The last communicator objects right/top take the rest of the physical geometry
    if (_comm->ThreadIdx()[0] == _comm->ThreadDim()[0] - 1)
      _bsize[0] += _size[0] % _comm->ThreadDim()[0];
    if (_comm->ThreadIdx()[1] == _comm->ThreadDim()[1] - 1)
      _bsize[1] += _size[1] % _comm->ThreadDim()[0];
    // Compute block length using block size
    _blength[0] = _h[0]*(_bsize[0] - 2);
    _blength[1] = _h[1]*(_bsize[1] - 2);
  }
  // Create global boundary halo
  _size[0] += 2;
  _size[1] += 2;
  // If !_comm --> default, see first constructor
  if (!_comm)
    _bsize = _size;
}
//------------------------------------------------------------------------------
/* Getter functions for some parameters
*/
const multi_index_t& Geometry::Size()        const {return _bsize;}
//------------------------------------------------------------------------------
const multi_index_t& Geometry::TotalSize()   const {return _size;}
//------------------------------------------------------------------------------
const multi_real_t&  Geometry::Length()      const {return _blength;}
//------------------------------------------------------------------------------
const multi_real_t&  Geometry::TotalLength() const {return _length;}
//------------------------------------------------------------------------------
const multi_real_t&  Geometry::Mesh()        const {return _h;}
//------------------------------------------------------------------------------
/// Updates the velocity field u on the boundary
// @param u grid for the velocity in x-direction
void Geometry::Update_U(Grid* u) const {
  // Copy all boundaries except of the physical boundary conditions
  _comm->copyBoundary(u);
  // Check if there are physical boundary conditions to be set
  BoundaryIterator bit(this);
  real_t vel_x = _velocity[0];
  // Update on bottom boundary
  if (_comm->isBottom()) {
    bit.SetBoundary(0);
    while (bit.Valid()) {
      u->Cell(bit) = -u->Cell(bit.Top());
      bit.Next();
    }
  }
  // Update on right boundary
  if (_comm->isRight()) {
    bit.SetBoundary(1);
    while (bit.Valid()) {
      u->Cell(bit)        = 0;
      u->Cell(bit.Left()) = 0;
      bit.Next();
    }
  }
  // Update on top boundary
  if (_comm->isTop()) {
    bit.SetBoundary(2);
    while (bit.Valid()) {
      u->Cell(bit) = 2*vel_x - u->Cell(bit.Down());
      bit.Next();
    }
  }
  // Update on left boundary
  if (_comm->isLeft()) {
    bit.SetBoundary(3);
    while (bit.Valid()) {
      u->Cell(bit) = 0;
      bit.Next();
    }
  }
}
//------------------------------------------------------------------------------
/// Updates the velocity field v on the boundary
// @param v grid for the velocity in y-direction
void Geometry::Update_V(Grid* v) const {
  // Copy all boundaries except of the physical boundary conditions
  _comm->copyBoundary(v);
  // Check if there are physical boundary conditions to be set
  BoundaryIterator bit(this);
  // Update on bottom boundary
  if (_comm->isBottom()) {
    bit.SetBoundary(0);
    while (bit.Valid()) {
      v->Cell(bit) = 0;
      bit.Next();
    }
  }
  // Update on right boundary
  if (_comm->isRight()) {
    bit.SetBoundary(1);
    while (bit.Valid()) {
      v->Cell(bit) = -v->Cell(bit.Left());
      bit.Next();
    }
  }
  // Update on top boundary
  if (_comm->isTop()) {
    bit.SetBoundary(2);
    while (bit.Valid()) {
      v->Cell(bit)        = 0;
      v->Cell(bit.Down()) = 0; 
      bit.Next();
    }
  }
  // Update on left boundary
  if (_comm->isLeft()) {
    bit.SetBoundary(3);
    while (bit.Valid()) {
      v->Cell(bit) = -v->Cell(bit.Right());
      bit.Next();
    }
  }
}
//------------------------------------------------------------------------------
/// Updates the pressure field p on the boundary
// @param p grid for the pressure
void Geometry::Update_P(Grid* p) const {
  // Copy all boundaries except of the physical boundary conditions
  _comm->copyBoundary(p);
  // Check if there are physical boundary conditions to be set
  BoundaryIterator bit(this);
  // Update on bottom boundary
  if (_comm->isBottom()) {
    bit.SetBoundary(0);
    while (bit.Valid()) {
      p->Cell(bit) = p->Cell(bit.Top());
      bit.Next();
    }
  }
  // Update on right boundary
  if (_comm->isRight()) {
    bit.SetBoundary(1);
    while (bit.Valid()) {
      p->Cell(bit) = p->Cell(bit.Left());
      bit.Next();
    }
  }
  // Update on top boundary
  if (_comm->isTop()) {
    bit.SetBoundary(2);
    while (bit.Valid()) {
      p->Cell(bit) = p->Cell(bit.Down());
      bit.Next();
    }
  }
  // Update on left boundary
  if (_comm->isLeft()) {
    bit.SetBoundary(3);
    while (bit.Valid()) {
      p->Cell(bit) = p->Cell(bit.Right());
      bit.Next();
    }
  }
}
//------------------------------------------------------------------------------
