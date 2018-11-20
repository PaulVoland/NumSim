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

#include "geometry.hpp"
#include "iterator.hpp"
#include "grid.hpp"
#include "communicator.hpp"
#include <stdio.h>
#include <string.h>
//------------------------------------------------------------------------------
Geometry::Geometry() : _comm(NULL) {
  _length[0] = 1.0;
  _length[1] = 1.0;
  _size[0] = 128;
  _size[1] = 128;
  _h[0] = _length[0] / _size[0];
  _h[1] = _length[1] / _size[1];
  _pressure = 0.0;
  _velocity[0] = 1.0;
  _velocity[1] = 0.0;

  // create boundary halo
  _size[0] += 2;
  _size[1] += 2;

  _bsize = _size;
  _blength = _length;
}
//------------------------------------------------------------------------------
Geometry::Geometry(const Communicator *comm) : _comm(comm) {
  _length[0] = 1.0;
  _length[1] = 1.0;
  _size[0] = 128;
  _size[1] = 128;
  _h[0] = _length[0] / _size[0];
  _h[1] = _length[1] / _size[1];
  _pressure = 0.0;
  _velocity[0] = 1.0;
  _velocity[1] = 0.0;

  _bsize[0] = _size[0] / _comm->ThreadDim()[0] + 2;
  _bsize[1] = _size[1] / _comm->ThreadDim()[1] + 2;
  if (_comm->ThreadIdx()[0] == _comm->ThreadDim()[0] - 1)
    _bsize[0] += _size[0] % _comm->ThreadDim()[0];
  if (_comm->ThreadIdx()[1] == _comm->ThreadDim()[1] - 1)
    _bsize[1] += _size[1] % _comm->ThreadDim()[0];

  _blength[0] = _h[0] * (_bsize[0] - 2);
  _blength[1] = _h[1] * (_bsize[1] - 2);

  // create boundary halo
  _size[0] += 2;
  _size[1] += 2;
}
//------------------------------------------------------------------------------
void Geometry::Load(const char *file) {
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
  _h[0] = _length[0] / _size[0];
  _h[1] = _length[1] / _size[1];

  _blength = _length;

  if (_comm) {
    _bsize[0] = _size[0] / _comm->ThreadDim()[0] + 2;
    _bsize[1] = _size[1] / _comm->ThreadDim()[1] + 2;
    if (_comm->ThreadIdx()[0] == _comm->ThreadDim()[0] - 1)
      _bsize[0] += _size[0] % _comm->ThreadDim()[0];
    if (_comm->ThreadIdx()[1] == _comm->ThreadDim()[1] - 1)
      _bsize[1] += _size[1] % _comm->ThreadDim()[0];

    _blength[0] = _h[0] * (_bsize[0] - 2);
    _blength[1] = _h[1] * (_bsize[1] - 2);
  }

  _size[0] += 2;
  _size[1] += 2;

  if (!_comm)
    _bsize = _size;
}
//------------------------------------------------------------------------------
const multi_index_t &Geometry::Size() const { return _bsize; }
//------------------------------------------------------------------------------
const multi_index_t &Geometry::TotalSize() const { return _size; }
//------------------------------------------------------------------------------
const multi_real_t &Geometry::Length() const { return _blength; }
//------------------------------------------------------------------------------
const multi_real_t &Geometry::TotalLength() const { return _length; }
//------------------------------------------------------------------------------
const multi_real_t &Geometry::Mesh() const { return _h; }

//------------------------------------------------------------------------------
/// Updates the velocity field u on the boundary
// @param u grid for the velocity in x-direction
void Geometry::Update_U(Grid *u) const {
  //copy all boundary except the physical BC
  _comm->copyBoundary(u);

  //check if there are physical BC to be set
  BoundaryIterator bit = BoundaryIterator(this);
  //falls mehrere Ränder betroffen sind, dann stimmt bit nicht mehr oder?

  if(_comm->isBottom()){
    bit.SetBoundary(0);
    while ( bit.Valid() ) {
      u->Cell(bit) = -1*u->Cell(bit.Top());
      bit.Next();
    }
  }

  if(_comm->isRight()){
    bit.SetBoundary(1);
    while (bit.Valid()) {
        u->Cell(bit)        = 0;
        u->Cell(bit.Left()) = 0; // necessary?!
        bit.Next();
    }
  }

  if(_comm->isTop()){
    bit.SetBoundary(2);
    while (bit.Valid()) {
        u->Cell(bit) = 2*vel_x - u->Cell(bit.Down());
        bit.Next();
    }
  }

  if(_comm->isLeft()){
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
void Geometry::Update_V(Grid *v) const {
  //copy all boundary except the physical BC
  _comm->copyBoundary(v);

  //check if there are physical BC to be set
  BoundaryIterator bit = BoundaryIterator(this);

  if(_comm->isBottom()){
    bit.SetBoundary(0);
    while (bit.Valid()) {
        v->Cell(bit) = 0;
        bit.Next();
    }
  }

  if(_comm->isRight()){
    bit.SetBoundary(1);
    while (bit.Valid()) {
        v->Cell(bit) = -v->Cell(bit.Left());
        bit.Next();
    }
  }

  if(_comm->isTop()){
    bit.SetBoundary(2);
    while (bit.Valid()) {
        v->Cell(bit) =        0;
        v->Cell(bit.Down()) = 0; // notwendig?!
        bit.Next();
    }
  }

  if(_comm->isLeft()){
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
void Geometry::Update_P(Grid *p) const {
  //copy all boundary except the physical BC
  _comm->copyBoundary(p);

  // Preparation of boundary iteration
  BoundaryIterator bit = BoundaryIterator(this);

  if(_comm->isBottom()){
    bit.SetBoundary(0);
    while (bit.Valid()) {
        p->Cell(bit) = p->Cell(bit.Top());
        bit.Next();
    }
  }

  if(_comm->isRight()){
    bit.SetBoundary(1);
    while (bit.Valid()) {
        p->Cell(bit) = p->Cell(bit.Left());
        bit.Next();
    }
  }

  if(_comm->isTop()){
    bit.SetBoundary(2);
    while (bit.Valid()) {
        p->Cell(bit) = p->Cell(bit.Down());
        bit.Next();
    }
  }

  if(_comm->isLeft()){
    bit.SetBoundary(3);
    while (bit.Valid()) {
        p->Cell(bit) = p->Cell(bit.Right());
        bit.Next();
    }
  }
}
//------------------------------------------------------------------------------
