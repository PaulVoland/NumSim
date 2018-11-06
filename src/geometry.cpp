#include "geometry.hpp"
#include "grid.hpp"
#include "iterator.hpp"
#include <iostream>
#include <cstdio>
#include <cstring>
#include <cstdlib>

using namespace std;
//------------------------------------------------------------------------------
/// Constructs a default geometry for the lid driven cavity problem
Geometry::Geometry() {
  // Set the initial velocities
  _velocity[0] = 1;
  _velocity[1] = 0;
  // Set the initial pressure
  _pressure    = 0;
  // Set the number of cells per row/column
  _size[0] = 128;
  _size[1] = 128;
  // Physical length of lid driven cavity
  _length[0] = 1;
  _length[1] = 1;
  // Resulting cell lengths
  _h[0] = _length[0]/_size[0];
  _h[1] = _length[1]/_size[1];
  // Message of success
  cout << "Loaded default geometry for lid driven cavity." << endl;
}

/// Loads a geometry setting from a specified file
// @param file filename including path information
void Geometry::Load(const char* file) {
  FILE* handle = fopen(file,"r");
  double inval[2];
  char name[20];
  while (!feof(handle)) {
    if (!fscanf(handle, "%s =", name)) continue;
    if (strcmp(name,"size") == 0) {
     if (fscanf(handle," %lf %lf\n",&inval[0],&inval[1])) {
        _size[0] = inval[0];
        _size[1] = inval[1];
      }
      continue;
    }
    if (strcmp(name,"length") == 0) {
      if (fscanf(handle," %lf %lf\n",&inval[0],&inval[1])) {
        _length[0] = inval[0];
        _length[1] = inval[1];
      }
      continue;
    }
    if (strcmp(name,"velocity") == 0) {
      if (fscanf(handle," %lf %lf\n",&inval[0],&inval[1])) {
        _velocity[0] = inval[0];
        _velocity[1] = inval[1];
      }
      continue;
    }
    if (strcmp(name,"pressure") == 0) {
      if (fscanf(handle," %lf\n",&inval[0]))
        _pressure = inval[0];
      continue;
    }
  }
  fclose(handle);
  // Set width values of the grid
  _h[0] = _length[0]/_size[0];
  _h[1] = _length[1]/_size[1];
  cout << "Loaded geometry data from file " << file << "." << endl;
}

/// Prints the geometry data
void Geometry::PrintData() {
  cout << "Showing geometry parameters..."                << endl;
  cout << "vel_x = "                      << _velocity[0] << endl;
  cout << "vel_y = "                      << _velocity[1] << endl;
  cout << "p = "                          << _pressure    << endl;
  cout << "cells_x = "                    << _size[0]     << endl;
  cout << "cells_y = "                    << _size[1]     << endl;
  cout << "l_x = "                        << _length[0]   << endl;
  cout << "l_y = "                        << _length[1]   << endl;
  cout << "h_x = "                        << _h[0]        << endl;
  cout << "h_y = "                        << _h[1]        << endl;
}

/// Getter functions for all class attributes
// Number of cells per dimension
const multi_index_t& Geometry::Size()  const {return _size;}
// Physical length of the domain
const multi_real_t& Geometry::Length() const {return _length;}
// Width of the resulting mesh
const multi_real_t& Geometry::Mesh()    const {return _h;}

/// Updates the velocity field u on the boundary
// @param u grid for the velocity in x-direction
void Geometry::Update_U(Grid* u) const {
  // Preparation of boundary iteration
  real_t vel_x = _velocity[0];
  BoundaryIterator bit = BoundaryIterator(this);

  // Update on bottom boundary
  bit.SetBoundary(0);
  while (bit.Valid()) {
    u->Cell(bit) = -u->Cell(bit.Top());
    bit.Next();
  }

  // Update on right boundary
  bit.SetBoundary(1);
  while (bit.Valid()) {
    u->Cell(bit)        = 0;
    u->Cell(bit.Left()) = 0; // necessary?!
    bit.Next();
  }

  // Update on top boundary
  bit.SetBoundary(2);
  while (bit.Valid()) {
    u->Cell(bit) = 2*vel_x - u->Cell(bit.Down());
    // cout << "Value of u on upper boundary: " << u->Cell(bit) << "\n" << endl;
    bit.Next();
  }

  // cout << "----------\n" << endl; 

  // Update on left boundary
  bit.SetBoundary(3);
  while (bit.Valid()) {
    u->Cell(bit) = 0;
    bit.Next();
  }
}

/// Updates the velocity field v on the boundary
// @param v grid for the velocity in y-direction
void Geometry::Update_V(Grid* v) const {
  // Preparation of boundary iteration
  BoundaryIterator bit = BoundaryIterator(this);

  // Update on bottom boundary
  bit.SetBoundary(0);
  while (bit.Valid()) {
    v->Cell(bit) = 0;
    bit.Next();
  }

  // Update on right boundary
  bit.SetBoundary(1);
  while (bit.Valid()) {
    v->Cell(bit) = -v->Cell(bit.Left());
    bit.Next();
  }

  // Update on top boundary
  bit.SetBoundary(2);
  while (bit.Valid()) {
    v->Cell(bit) =        0;
    v->Cell(bit.Down()) = 0; // notwendig?!
    bit.Next();
  }

  // Update on left boundary
  bit.SetBoundary(3);
  while (bit.Valid()) {
    v->Cell(bit) = -v->Cell(bit.Right());
    bit.Next();
  }
}

/// Updates the pressure field p on the boundary
// @param p grid for the pressure
void Geometry::Update_P(Grid* p) const {
  // Preparation of boundary iteration
  BoundaryIterator bit = BoundaryIterator(this);

  // Update on bottom boundary
  bit.SetBoundary(0);
  while (bit.Valid()) {
    p->Cell(bit) = p->Cell(bit.Top());
    bit.Next();
  }

 // Update on right boundary
  bit.SetBoundary(1);
  while (bit.Valid()) {
    p->Cell(bit) = p->Cell(bit.Left());
    bit.Next();
  }

  // Update on top boundary
  bit.SetBoundary(2);
  while (bit.Valid()) {
    p->Cell(bit) = p->Cell(bit.Down());
    bit.Next();
  }

  // Update on left boundary
  bit.SetBoundary(3);
  while (bit.Valid()) {
    p->Cell(bit) = p->Cell(bit.Right());
    bit.Next();
  }
}
//------------------------------------------------------------------------------