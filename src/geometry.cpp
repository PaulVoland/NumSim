#include "geometry.hpp"
#include "iterator.hpp"
#include "grid.hpp"
#include "communicator.hpp"
#include <iostream>
#include <stdio.h>
#include <string.h>

using namespace std;
//------------------------------------------------------------------------------
void Geometry::UpdateCellDirichlet_U(Grid* u, const real_t& value,
    const Iterator& it) const {
  switch (_cell[_boffset + it.Pos()[0] + it.Pos()[1]*_size[0]].fluid) {
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
    break;
  case cellSE:
    u->Cell(it) = value;
    break;
  default:
    u->Cell(it) = value;
    break;
  };
}
//------------------------------------------------------------------------------
void Geometry::UpdateCellDirichlet_V(Grid* v, const real_t& value,
    const Iterator& it) const {
  switch (_cell[_boffset + it.Pos()[0] + it.Pos()[1]*_size[0]].fluid) {
  case cellN:
    v->Cell(it) = value;
    break;
  case cellW:
    v->Cell(it) = 2.0*value - v->Cell(it.Left());
    break;
  case cellNW:
    v->Cell(it) = value;
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
  switch (_cell[_boffset + it.Pos()[0] + it.Pos()[1]*_size[0]].fluid) {
  case cellN:
    grid->Cell(it) = grid->Cell(it.Top());
    break;
  case cellW:
    grid->Cell(it) = grid->Cell(it.Left());
    break;
  case cellNW:
    grid->Cell(it) = 0.5*(grid->Cell(it.Left()) + grid->Cell(it.Top()));
    break;
  case cellS:
    grid->Cell(it) = grid->Cell(it.Down());
    break;
  case cellSW:
    grid->Cell(it) = 0.5*(grid->Cell(it.Left()) + grid->Cell(it.Down()));
    break;
  case cellE:
    grid->Cell(it) = grid->Cell(it.Right());
    break;
  case cellNE:
    grid->Cell(it) = 0.5*(grid->Cell(it.Right()) + grid->Cell(it.Top()));
    break;
  case cellSE:
    grid->Cell(it) = 0.5*(grid->Cell(it.Right()) + grid->Cell(it.Down()));
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
void Geometry::UpdateCellDirichlet_T(Grid* T, const real_t& value,
    const Iterator& it) const {
  switch (_cell[_boffset + it.Pos()[0] + it.Pos()[1]*_size[0]].fluid) {
  case cellN:
    T->Cell(it) = 2.0*value - T->Cell(it.Top());
    break;
  case cellW:
    T->Cell(it) = 2.0*value - T->Cell(it.Left());
    break;
  case cellNW:
    T->Cell(it) = 2.0*value - (T->Cell(it.Top()) + T_>Cell(it.Left()))/2.0;
    break;
  case cellS:
    T->Cell(it) = 2.0*value - T->Cell(it.Down());
    break;
  case cellSW:
    T->Cell(it) = 2.0*value - (T->Cell(it.Down()) + T_>Cell(it.Left()))/2.0;
    break;
  case cellE:
    T->Cell(it) = 2.0*value - T->Cell(it.Right());
    break;
  case cellNE:
    T->Cell(it) = 2.0*value - (T->Cell(it.Top()) + T_>Cell(it.Right()))/2.0;
    break;
  case cellSE:
    T->Cell(it) = 2.0*value - (T->Cell(it.Down()) + T_>Cell(it.Right()))/2.0;
    break;
  default:
    u->Cell(it) = value;
    break;
  };
}
//------------------------------------------------------------------------------
/// Constructs a default geometry
Geometry::Geometry() : _comm(NULL) {
  _length[0] = 1.0;
  _length[1] = 1.0;
  _size[0] = 128;
  _size[1] = 128;
  _h[0] = _length[0]/_size[0];
  _h[1] = _length[1]/_size[1];
  _pressure    = 0.0;
  _temperature = 0.0;
  _velocity[0] = 1.0;
  _velocity[1] = 0.0;
  // Create boundary halo
  _size[0] += 2;
  _size[1] += 2;
  // Block is identical to whole domain
  _bsize   = _size;
  _blength = _length;
  _cell    = NULL;
  _boffset = 0;
  // Message of success
  cout << "Loaded default geometry for lid driven cavity without communicator." << endl;
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
  _pressure    = 0.0;
  _temperature = 0.0;
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
  _cell    =  NULL;
  _boffset =  0;
  // Message of success
  cout << "Loaded default geometry for lid driven cavity with a communicator object." << endl;
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
        for (uint32_t y = _size[1]; y-- > 0;) {
          if (feof(handle)) {
            delete[] _cell;
            _cell = NULL;
            break;
          }
          for (uint32_t x = 0; x < _size[0]; ++x) {
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
        for (uint32_t y = 0; y < _size[1]; ++y) {
          for (uint32_t x = 0; x < _size[0]; ++x) {
            uint32_t check = 0;
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
  // If !_comm --> default, see first constructor
  _blength = _length;
  _boffset = 0;
  // If communicator objects already exist
  if (_comm) {
    // Make blocks for each communicator object including halo
    _bsize[0] = _size[0]/_comm->ThreadDim()[0] + 2;
    _bsize[1] = _size[1]/_comm->ThreadDim()[1] + 2;
    // Golbal starting index of cell field for each thread
    _boffset = (_bsize[0] - 2)*_comm->ThreadIdx()[0] +
      (_size[0] + 2)*(_bsize[1] - 2)*_comm->ThreadIdx()[1];
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
  cout << "Loaded geometry data from file " << file << "." << endl;
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
const Cell_t& Geometry::Cell(const Iterator& it) const {
  return _cell[it]; // Uses cast command via Iterator::operator
}
//------------------------------------------------------------------------------
const real_t&        Geometry::Temperature() const {return _temperature;}
//------------------------------------------------------------------------------
/// Updates the velocity field u on the boundary
// @param u grid for the velocity in x-direction
void Geometry::Update_U(Grid* u) const {
  if (_cell) {
    /// Cell_t is used for free geometries
    Iterator it(this);
    while (it.Valid()) {
      switch (_cell[_boffset + it.Pos()[0] + it.Pos()[1]*_size[0]].type) {
      case typeSolid:
        UpdateCellDirichlet_U(u, 0.0, it);
        break;
      case typeIn:
        UpdateCellDirichlet_U(u, _velocity[0], it);
        break;
      case typeInH:
        UpdateCellDirichlet_U(u, 0.0, it);
        break;
      case typeInV:
        UpdateCellDirichlet_U(u, _velocity[0]*
          _cell[_boffset + it.Pos()[0] + it.Pos()[1]*_size[0]].factor, it);
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
}
//------------------------------------------------------------------------------
/// Updates the velocity field v on the boundary
// @param v grid for the velocity in y-direction
void Geometry::Update_V(Grid* v) const {
  if (_cell) {
    /// Cell_t is used for free geometries
    Iterator it(this);
    while (it.Valid()) {
      switch (_cell[_boffset + it.Pos()[0] + it.Pos()[1]*_size[0]].type) {
      case typeSolid:
        UpdateCellDirichlet_V(v, 0.0, it);
        break;
      case typeIn:
        UpdateCellDirichlet_V(v, _velocity[1], it);
        break;
      case typeInH:
        UpdateCellDirichlet_V(v, _velocity[1]*
          _cell[_boffset + it.Pos()[0] + it.Pos()[1]*_size[0]].factor, it);
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
}
//------------------------------------------------------------------------------
/// Updates the pressure field p on the boundary
// @param p grid for the pressure
void Geometry::Update_P(Grid* p) const {
  if (_cell) {
    /// Cell_t is used for free geometries
    Iterator it(this);
    while (it.Valid()) {
      switch (_cell[_boffset + it.Pos()[0] + it.Pos()[1]*_size[0]].type) {
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
        p->Cell(it) = _pressure; // why?
        // UpdateCellNeumann_P(p, it); // alternative version
      case typeSlipV:
        p->Cell(it) = _pressure; // why?
        // UpdateCellNeumann_P(p, it); // alternative version
        break;
      case typeOut:
        p->Cell(it) = 0.0;
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
      switch (_cell[_boffset + it.Pos()[0] + it.Pos()[1]*_size[0]].type) {
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
        T->Cell(it) = _temperature; // why?
        // UpdateCellNeumann(T, it); // alternative version
      case typeSlipV:
        T->Cell(it) = _temperature; // why?
        // UpdateCellNeumann(T, it); // alternative version
        break;
      case typeOut:
        T->Cell(it) = 0.0;
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
    if (_comm->isBottom()) {
      bit.SetBoundary(0);
      while (bit.Valid()) {
        T->Cell(bit) = T->Cell(bit.Top());
        bit.Next();
      }
    }
    // Update on right boundary
    if (_comm->isRight()) {
      bit.SetBoundary(1);
      while (bit.Valid()) {
      T->Cell(bit) = T->Cell(bit.Left());
      bit.Next();
      }
    }
    // Update on top boundary
    if (_comm->isTop()) {
      bit.SetBoundary(2);
      while (bit.Valid()) {
        T->Cell(bit) = T->Cell(bit.Down());
        bit.Next();
      }
    }
    // Update on left boundary
    if (_comm->isLeft()) {
      bit.SetBoundary(3);
      while (bit.Valid()) {
      T->Cell(bit) = T->Cell(bit.Right());
      bit.Next();
      }
    }
  }
}
//------------------------------------------------------------------------------