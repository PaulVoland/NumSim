#include "geometry.hpp"
#include "iterator.hpp"
#include "grid.hpp"
#include <iostream>
#include <stdio.h>
#include <string.h>

using namespace std;
//------------------------------------------------------------------------------
void Geometry::UpdateCellDirichlet_U(Grid* u, const real_t& value,
    const Iterator& it) const {
  switch (_cell[it.Pos()[0] + it.Pos()[1]*_size[0]].fluid) {
  case cellN:
    u->Cell(it) = 2.0*value - u->Cell(it.Top());
    break;
  case cellW:
    u->Cell(it.Left()) = value;
    u->Cell(it) = value;
    break;
  case cellNW:
    u->Cell(it.Left()) = value;
    u->Cell(it) = 2.0*value - u->Cell(it.Top());
    break;
  case cellS:
    u->Cell(it) = 2.0*value - u->Cell(it.Down());
    break;
  case cellSW:
    u->Cell(it.Left()) = value;
    u->Cell(it) = 2.0*value - u->Cell(it.Down());
    break;
  case cellE:
    u->Cell(it) = value;
    break;
  case cellNE:
    u->Cell(it) = value;
    // u->Cell(it) = value - (u->Cell(it.Top()) - u->Cell(it.Right()))/2.0;
    // Alternative version Larissa
    break;
  case cellSE:
    u->Cell(it) = value;
    // u->Cell(it) = value - (u->Cell(it.Down()) - u->Cell(it.Right()))/2.0;
    // Alternative version Larissa
    break;
  default:
    u->Cell(it) = value;
    break;
  };
}
//------------------------------------------------------------------------------
void Geometry::UpdateCellDirichlet_V(Grid* v, const real_t& value,
    const Iterator& it) const {
  switch (_cell[it.Pos()[0] + it.Pos()[1]*_size[0]].fluid) {
  case cellN:
    v->Cell(it) = value;
    break;
  case cellW:
    v->Cell(it) = 2.0*value - v->Cell(it.Left());
    break;
  case cellNW:
    v->Cell(it) = value;
    // v->Cell(it) = value - (v->Cell(it.Left()) - v->Cell(it.Top()))/2.0;
    // Alternative version Larissa
    break;
  case cellS:
    v->Cell(it.Down()) = value;
    v->Cell(it) = value;
    break;
  case cellSW:
    v->Cell(it.Down()) = value;
    v->Cell(it) = 2.0*value - v->Cell(it.Left());
    break;
  case cellE:
    v->Cell(it) = 2.0*value - v->Cell(it.Right());
    break;
  case cellNE:
    v->Cell(it) = value;
    // v->Cell(it) = value - (v->Cell(it.Right()) - v->Cell(it.Top()))/2.0;
    // Alternative version Larissa
    break;
  case cellSE:
    v->Cell(it.Down()) = value;
    v->Cell(it) = 2.0*value - v->Cell(it.Right());
    break;
  default:
    v->Cell(it) = value;
    break;
  };
}
//------------------------------------------------------------------------------
void Geometry::UpdateCellNeumann(Grid* grid, const Iterator& it) const {
  switch (_cell[it.Pos()[0] + it.Pos()[1]*_size[0]].fluid) {
  case cellN:
    grid->Cell(it) = grid->Cell(it.Top());
    break;
  case cellW:
    grid->Cell(it) = grid->Cell(it.Left());
    break;
  case cellNW:
    grid->Cell(it) = (grid->Cell(it.Left()) + grid->Cell(it.Top()))/2.0;
    break;
  case cellS:
    grid->Cell(it) = grid->Cell(it.Down());
    break;
  case cellSW:
    grid->Cell(it) = (grid->Cell(it.Left()) + grid->Cell(it.Down()))/2.0;
    break;
  case cellE:
    grid->Cell(it) = grid->Cell(it.Right());
    break;
  case cellNE:
    grid->Cell(it) = (grid->Cell(it.Right()) + grid->Cell(it.Top()))/2.0;
    break;
  case cellSE:
    grid->Cell(it) = (grid->Cell(it.Right()) + grid->Cell(it.Down()))/2.0;
    break;
  default:
    break;
  };
}
//------------------------------------------------------------------------------
void Geometry::UpdateCellNeumann_P(Grid* grid, const Iterator& it) const {
  UpdateCellNeumann(grid, it);
}
//------------------------------------------------------------------------------
// Not used in Larissas version
void Geometry::UpdateCellDirichlet_T(Grid* T, const real_t& value,
    const Iterator& it) const {
  switch (_cell[it.Pos()[0] + it.Pos()[1]*_size[0]].fluid) {
  case cellN:
    T->Cell(it) = 2.0*value - T->Cell(it.Top());
    break;
  case cellW:
    T->Cell(it) = 2.0*value - T->Cell(it.Left());
    break;
  case cellNW:
    T->Cell(it) = 2.0*value - (T->Cell(it.Top()) + T->Cell(it.Left()))/2.0;
    break;
  case cellS:
    T->Cell(it) = 2.0*value - T->Cell(it.Down());
    break;
  case cellSW:
    T->Cell(it) = 2.0*value - (T->Cell(it.Down()) + T->Cell(it.Left()))/2.0;
    break;
  case cellE:
    T->Cell(it) = 2.0*value - T->Cell(it.Right());
    break;
  case cellNE:
    T->Cell(it) = 2.0*value - (T->Cell(it.Top()) + T->Cell(it.Right()))/2.0;
    break;
  case cellSE:
    T->Cell(it) = 2.0*value - (T->Cell(it.Down()) + T->Cell(it.Right()))/2.0;
    break;
  default:
    T->Cell(it) = value;
    break;
  };
}
//------------------------------------------------------------------------------
/// Constructs a default geometry
Geometry::Geometry() {
  _length[0] = 1.0;
  _length[1] = 1.0;
  _size[0] = 50;
  _size[1] = 50;
  _h[0] = _length[0]/_size[0];
  _h[1] = _length[1]/_size[1];
  _velocity[0] = 0.0;
  _velocity[1] = 0.0;
  _pressure    = 0.0;
  _temperature = 0.0;
  // Create boundary halo
  _size[0] += 2;
  _size[1] += 2;
  // No cell field set yet
  _cell    = NULL;
  _num_coupling = 0;
  // Message of success
  cout << "Loaded default geometry for lid driven cavity." << endl;
}
//------------------------------------------------------------------------------
/// Destructor (deletes the cell field information)
Geometry::~Geometry() {
  if (_cell)
    delete[] _cell;
}
//------------------------------------------------------------------------------
/// Loads a geometry from a file
void Geometry::Load(const char* file) {
  FILE *handle = fopen(file, "r");
  double inval[2];
  char name[200000];
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
    if (strcmp(name, "temperature") == 0) {
      if (fscanf(handle, " %lf\n", &inval[0]))
        _temperature = inval[0];
      continue;
    }
    if (strcmp(name, "geometry") == 0) {
      if (!fscanf(handle, "%s\n", name))
        continue;
      if (strcmp(name, "free") == 0) {
        if (_cell)
          delete[] _cell;
        _cell = new Cell_t[_size[0]*_size[1]];
        bool parabolic = false;
        // Read stuff from file
        for (int y = _size[1]; y-- > 0;) {
          if (feof(handle)) {
            delete[] _cell;
            _cell = NULL;
            break;
          }
          for (int x = 0; x < _size[0]; ++x) {
            _cell[x + y*_size[0]].fluid = cellNone;
            _cell[x + y*_size[0]].factor = 1.0;
            switch (getc(handle)) {
            case '#':
              _cell[x + y*_size[0]].type = typeSolid;
              break;
            case 'I':
              _cell[x + y*_size[0]].type = typeIn;
              break;
            case 'H':
              _cell[x + y*_size[0]].type = typeInH;
              parabolic = true;
              break;
            case 'V':
              _cell[x + y*_size[0]].type = typeInV;
              parabolic = true;
              break;
            case '-':
              _cell[x + y*_size[0]].type = typeSlipH;
              break;
            case '|':
              _cell[x + y*_size[0]].type = typeSlipV;
              break;
            case 'O':
              _cell[x + y*_size[0]].type = typeOut;
              break;
            case 'h':
              _cell[x + y*_size[0]].type = typeTDir_h;
              break;
            case 'c':
              _cell[x + y*_size[0]].type = typeTDir_c;
              break;
            case 'C':
              _cell[x + y*_size[0]].type = typeCoupling;
              _num_coupling++;
              break;
            default: // All other cases, box for bottom/top/left/right layer
              if (x == 0 || x == _size[0] - 1 || y == 0 || y == _size[1] - 1)
                _cell[x + y*_size[0]].type = typeSolid;
              else
                _cell[x + y*_size[0]].type = typeFluid;
              break;
            };
          }
          if (!fscanf(handle, "\n"))
            continue;
        }
        if (!_cell)
          break;
        // Process it
        for (int y = 0; y < _size[1]; ++y) {
          for (int x = 0; x < _size[0]; ++x) {
            int check = 0;
            if (_cell[x + y*_size[0]].type == typeFluid)
              continue;
            if (x < _size[0] - 1 && _cell[x + 1 + y*_size[0]].type == typeFluid)
              check |= 8;
            if (x > 0 && _cell[x - 1 + y*_size[0]].type == typeFluid)
              check |= 2;
            if (y < _size[1] - 1 && _cell[x + (y + 1)*_size[0]].type == typeFluid)
              check |= 1;
            if (y > 0 && _cell[x + (y - 1)*_size[0]].type == typeFluid)
              check |= 4;
            switch (check) {
            case 5:
            case 7:
            case 10:
            case 11:
            case 13:
            case 14:
            case 15:
              _cell[x + y*_size[0]].type = typeFluid;
              _cell[x + y*_size[0]].fluid = cellNone;
              // Remove single 'obstacle' cell; reset iteration to x--, y--
              if (x > 0)
                x--;
              if (y > 0)
                y--;
              break;
            case 1:
              _cell[x + y*_size[0]].fluid = cellN;
              break;
            case 2:
              _cell[x + y*_size[0]].fluid = cellW;
              break;
            case 3:
              _cell[x + y*_size[0]].fluid = cellNW;
              break;
            case 4:
              _cell[x + y*_size[0]].fluid = cellS;
              break;
            case 6:
              _cell[x + y*_size[0]].fluid = cellSW;
              break;
            case 8:
              _cell[x + y*_size[0]].fluid = cellE;
              break;
            case 9:
              _cell[x + y*_size[0]].fluid = cellNE;
              break;
            case 12:
              _cell[x + y*_size[0]].fluid = cellSE;
              break;
            };
          }
        }
        // Parabolic stuff
        if (parabolic) {
          for (int y = 0; y < _size[1]; ++y) {
            for (int x = 0; x < _size[0]; ++x) {
              int32_t dist1 = 0;
              int32_t dist2 = 0;
              switch (_cell[x + y*_size[0]].type) {
              case typeInH:
                while (x - dist1 >= 0 && _cell[x - dist1 + y*_size[0]].type == typeInH)
                  ++dist1;
                while (x + dist2 < _size[0] && _cell[x + dist2 + y*_size[0]].type == typeInH)
                  ++dist2;
                _cell[x + y*_size[0]].factor = 4.0*((real_t)(dist1) - 0.5)*((real_t)(dist2) - 0.5)/
                  (real_t)((dist1 + dist2 - 1)*(dist1 + dist2 - 1));
                break;
              case typeInV:
                while (y - dist1 >= 0 && _cell[x + (y - dist1)*_size[0]].type == typeInV)
                  ++dist1;
                while (y + dist2 < _size[1] && _cell[x + (dist2 + y)*_size[0]].type == typeInV)
                  ++dist2;
                _cell[x + y*_size[0]].factor = 4.0*((real_t)(dist1) - 0.5)*((real_t)(dist2) - 0.5)/
                  (real_t)((dist1 + dist2 - 1)*(dist1 + dist2 - 1));
                break;
              default:
                break;
              };
            }
          }
        }
        _size[0] -= 2;
        _size[1] -= 2;
      }
    }
  }
  fclose(handle);
  // Set width values of the grid
  _h[0] = _length[0]/_size[0];
  _h[1] = _length[1]/_size[1];
  // Create global boundary halo
  _size[0] += 2;
  _size[1] += 2;
  cout << "Loaded geometry data from file " << file << "." << endl;
}
//------------------------------------------------------------------------------
/* Getter functions for some parameters
*/
const multi_index_t& Geometry::TotalSize()   const {return _size;}
//------------------------------------------------------------------------------
const multi_real_t&  Geometry::TotalLength() const {return _length;}
//------------------------------------------------------------------------------
const multi_real_t&  Geometry::Mesh()        const {return _h;}
//------------------------------------------------------------------------------
const Cell_t& Geometry::Cell(const Iterator& it) const {
  return _cell[it]; // Uses cast command via Iterator::operator
}
//------------------------------------------------------------------------------
const multi_real_t&  Geometry::Velocity()    const {return _velocity;}
//------------------------------------------------------------------------------
const real_t&        Geometry::Pressure()    const {return _pressure;}
//------------------------------------------------------------------------------
const real_t&        Geometry::Temperature() const {return _temperature;}

const index_t& Geometry::Num_Coupling() const{return _num_coupling;}
//------------------------------------------------------------------------------
/// Updates the velocity field u on the boundary
// @param u     grid for the velocity in x-direction
// @param u_Dir Dirichlet value for u from .param
void Geometry::Update_U(Grid* u, const real_t& u_Dir) const {
  if (_cell) {
    /// Cell_t is used for free geometries
    Iterator it(this);
    while (it.Valid()) {
      switch (_cell[it.Pos()[0] + it.Pos()[1]*_size[0]].type) {
      case typeSolid:
        UpdateCellDirichlet_U(u, 0.0, it);
        break;
      case typeIn:
        UpdateCellDirichlet_U(u, u_Dir, it);
        break;
      case typeInH:
        UpdateCellDirichlet_U(u, 0.0, it);
        break;
      case typeInV:
        UpdateCellDirichlet_U(u, u_Dir*
          _cell[it.Pos()[0] + it.Pos()[1]*_size[0]].factor, it);
        break;
      case typeSlipH:
        UpdateCellNeumann(u, it);
        break;
      case typeSlipV:
        UpdateCellDirichlet_U(u, 0.0, it);
        break;
      case typeOut:
        UpdateCellNeumann(u, it);
        break;
      case typeTDir_h:
        UpdateCellDirichlet_U(u, 0.0, it);
        break;
      case typeTDir_c:
        UpdateCellDirichlet_U(u, 0.0, it);
        break;
      default:
        break;
      };
      it.Next();
    }
  } else {
    /// Default lid driven cavity example
    // Check if there are physical boundary conditions to be set
    BoundaryIterator bit(this);
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
      u->Cell(bit.Left()) = 0;
      bit.Next();
    }
    // Update on top boundary
    bit.SetBoundary(2);
    while (bit.Valid()) {
      u->Cell(bit) = 2*u_Dir - u->Cell(bit.Down());
      bit.Next();
    }
    // Update on left boundary
    bit.SetBoundary(3);
    while (bit.Valid()) {
      u->Cell(bit) = 0;
      bit.Next();
    }
  }
}
//------------------------------------------------------------------------------
/// Updates the velocity field v on the boundary
// @param v     grid for the velocity in y-direction
// @param v_Dir Dirichlet value for v from .param
void Geometry::Update_V(Grid* v, const real_t& v_Dir) const {
  if (_cell) {
    /// Cell_t is used for free geometries
    Iterator it(this);
    while (it.Valid()) {
      switch (_cell[it.Pos()[0] + it.Pos()[1]*_size[0]].type) {
      case typeSolid:
        UpdateCellDirichlet_V(v, 0.0, it);
        break;
      case typeIn:
        UpdateCellDirichlet_V(v, v_Dir, it);
        break;
      case typeInH:
        UpdateCellDirichlet_V(v, v_Dir*
          _cell[it.Pos()[0] + it.Pos()[1]*_size[0]].factor, it);
        break;
      case typeInV:
        UpdateCellDirichlet_V(v, 0.0, it);
        break;
      case typeSlipH:
        UpdateCellDirichlet_V(v, 0.0, it);
        break;
      case typeSlipV:
        UpdateCellNeumann(v, it);
        break;
      case typeOut:
        UpdateCellNeumann(v, it);
        break;
      case typeTDir_h:
        UpdateCellDirichlet_V(v, 0.0, it);
        break;
      case typeTDir_c:
        UpdateCellDirichlet_V(v, 0.0, it);
        break;
      default:
        break;
      };
      it.Next();
    }
  } else {
    /// Default lid driven cavity example
    // Check if there are physical boundary conditions to be set
    BoundaryIterator bit(this);
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
      v->Cell(bit)        = 0;
      v->Cell(bit.Down()) = 0;
      bit.Next();
    }
    // Update on left boundary
    bit.SetBoundary(3);
    while (bit.Valid()) {
      v->Cell(bit) = -v->Cell(bit.Right());
      bit.Next();
    }
  }
}
//------------------------------------------------------------------------------
/// Updates the pressure field p on the boundary
// @param p     grid for the pressure
// @param p_Dir Dirichlet value for p from .param
void Geometry::Update_P(Grid* p, const real_t& p_Dir) const {
  if (_cell) {
    /// Cell_t is used for free geometries
    Iterator it(this);
    while (it.Valid()) {
      switch (_cell[it.Pos()[0] + it.Pos()[1]*_size[0]].type) {
      case typeSolid:
        UpdateCellNeumann_P(p, it);
        break;
      case typeIn:
        UpdateCellNeumann_P(p, it);
        break;
      case typeInH:
        UpdateCellNeumann_P(p, it);
        break;
      case typeInV:
        UpdateCellNeumann_P(p, it);
        break;
      case typeSlipH:
        UpdateCellNeumann_P(p, it);
        // p->Cell(it) = p_Dir; // Alternative version 1
        break;
      case typeSlipV:
        UpdateCellNeumann_P(p, it);
        // p->Cell(it) = p_Dir; // Alternative version 1
        break;
      case typeOut:
        UpdateCellNeumann_P(p, it);
        // p->Cell(it) = 0.0; // Alternative version 1
        break;
      case typeTDir_h:
        UpdateCellNeumann_P(p, it);
        break;
      case typeTDir_c:
        UpdateCellNeumann_P(p, it);
        break;
      default:
        break;
      };
      it.Next();
    }
  } else {
    /// Default lid driven cavity example
    // Check if there are physical boundary conditions to be set
    BoundaryIterator bit(this);
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
}
//------------------------------------------------------------------------------
/// Updates the temperature field T on the boundary
// @param T   grid for the temperature
// @param T_h higher temeparture value from .param
// @param T_c lower temperature vlaue from .param
void Geometry::Update_T(Grid* T, const real_t& T_h, const real_t& T_c) const {
  if (_cell) {
    /// Cell_t is used for free geometries
    Iterator it(this);
    while (it.Valid()) {
      switch (_cell[it.Pos()[0] + it.Pos()[1]*_size[0]].type) {
      case typeSolid:
        UpdateCellNeumann(T, it);
        break;
      case typeIn:
        UpdateCellNeumann(T, it);
        break;
      case typeInH:
        UpdateCellNeumann(T, it);
        break;
      case typeInV:
        UpdateCellNeumann(T, it);
        break;
      case typeSlipH:
        UpdateCellNeumann(T, it);
        // T->Cell(it) = T_h; // Alternative version 1
        // T->Cell(it) = T_c; // Alternative version 2
        break;
       case typeSlipV:
        UpdateCellNeumann(T, it);
        // T->Cell(it) = T_h; // Alternative version 1
        // T->Cell(it) = T_c; // Alternative version 2
        break;
      case typeOut:
        UpdateCellNeumann(T, it);
        // T->Cell(it) = 0.0; // Alternative version 1
        break;
      case typeTDir_h: // possibly use T_h as positive offset to TI (which is null niveau)
        UpdateCellDirichlet_T(T, T_h, it);
        break;
      case typeTDir_c: // possibly use T_c as negative offset to TI (which is null niveau)
        UpdateCellDirichlet_T(T, T_c, it);
        break;
      default:
        break;
      };
      it.Next();
    }
  } else {
    /// Default lid driven cavity example
    // Check if there are physical boundary conditions to be set
    BoundaryIterator bit(this);
    // Update on bottom boundary
    bit.SetBoundary(0);
    while (bit.Valid()) {
      T->Cell(bit) = T->Cell(bit.Top());
      bit.Next();
    }
    // Update on right boundary
    bit.SetBoundary(1);
    while (bit.Valid()) {
      T->Cell(bit) = T->Cell(bit.Left());
      bit.Next();
    }
    // Update on top boundary
    bit.SetBoundary(2);
    while (bit.Valid()) {
      T->Cell(bit) = T->Cell(bit.Down());
      bit.Next();
    }
    // Update on left boundary
    bit.SetBoundary(3);
    while (bit.Valid()) {
      T->Cell(bit) = T->Cell(bit.Right());
      bit.Next();
    }
  }
}
//------------------------------------------------------------------------------
