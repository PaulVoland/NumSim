#include "compute.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
#include "grid.hpp"
#include "solver.hpp"
#include "iterator.hpp"

#include <iostream>
#include <cmath>

using namespace std;
//------------------------------------------------------------------------------
/// Creates a compute instance with given geometry and parameter
//  @param geom  given geometry
//  @param param given parameter data
Compute::Compute(const Geometry* geom, const Parameter* param) {
  // Initialize
  _geom   = geom;
  _param  = param;
  _solver = new SOR(geom, param->Omega());

  // Set timestep data
  real_t _t = 0.0;
  _dtlimit  = param->Dt();
  
  // Tolerance for Poisson solver
  _epslimit = param->Eps();

  // Offsets for the grid variables
  multi_real_t off_u;
  multi_real_t off_v;
  multi_real_t off_p;

  off_u[0] = geom->Mesh()[0];
  off_u[1] = geom->Mesh()[1]/2;
  off_v[0] = geom->Mesh()[0]/2;
  off_v[1] = geom->Mesh()[1];
  off_p[0] = geom->Mesh()[0]/2;
  off_p[1] = geom->Mesh()[1]/2;

  // Instantiate grids with offsets
  _u   = new Grid(geom, off_u);
  _v   = new Grid(geom, off_v);
  _p   = new Grid(geom, off_p);
  _F   = new Grid(geom);
  _G   = new Grid(geom);
  _rhs = new Grid(geom);
  _tmp = new Grid(geom);

  // Set initial values for physical variables
  _u->Initialize(0.0);
  _v->Initialize(0.0);
  _p->Initialize(0.0);
}

/// Deletes all grids
Compute::~Compute() {
  delete[] _u;
  delete[] _v;
  delete[] _p;

  delete[] _F;
  delete[] _G;
  delete[] _rhs;

  delete[] _tmp;
  delete _solver;
}

/// Execute one time step of the fluid simulation (with or without debug info)
// @ param printInfo print information about current solver state (residual etc.)
void Compute::TimeStep(bool printInfo) {
  // Refresh boundary values
  _geom->Update_U(_u);
  _geom->Update_V(_v);
  _geom->Update_P(_p);

  // Find timestep as a minimum of different criteria
  real_t dt_1 = _param->Tau()*fmin(_geom->Mesh()[0], _geom->Mesh()[1])/
    fmax(_u->AbsMax(), _v->AbsMax());
  real_t dt_2 = _param->Tau()*_param->Re()*
    (_geom->Mesh()[0]*_geom->Mesh()[0]*_geom->Mesh()[1]*_geom->Mesh()[1])/
    (2*(_geom->Mesh()[0]*_geom->Mesh()[0] + _geom->Mesh()[1]*_geom->Mesh()[1]));
  real_t dt = min(dt_1, dt_2);
  // If explicit timestep > 0 is given in the paramater file, use this instead
  if (_param->Dt() > 0) {
    dt = _param->Dt();
  }

  // Compute temporary velocities F, G using difference schemes
  MomentumEqu(dt);

  // Compute RHS of the Poisson equation
  RHS(dt);

  // Find solution of the Poisson equation using a SOR solver
  real_t res = 1000000.0;
  index_t i  = 0;
  while (res > _param->Eps() && i < _param->IterMax()) {
    res = _solver->Cycle(_p, _rhs);
    
    i++;
  }

  // Compute 'new' velocities using the pressure
  NewVelocities(dt);

  // (optionally) printing informations
  if (printInfo) {
    cout << "_t = " << fixed << _t << "\tdt = " << scientific << dt << " \tres = " << res
      << " \tprogress: " << fixed << (uint32_t) 100*_t/_param->Tend() << "%\n" << endl;
  }

  // Next timestep
  _t += dt;
}
/* Getter functions
*/
/// Returns the simulated time in total
const real_t& Compute::GetTime() const {return _t;}

/// Returns the pointer to u
const Grid* Compute::GetU() const {return _u;}

// Returns the pointer to v
const Grid* Compute::GetV() const {return _v;}

// Returns the pointer to p
const Grid* Compute::GetP() const {return _p;}

// Returns the pointer to the RHS
const Grid* Compute::GetRHS() const {return _rhs;}

/// Computes and returns the absolute velocity
const Grid* Compute::GetVelocity() {
  // Initialize full Iterator
  Iterator it(_geom);

  // Go through all cells
  while (it.Valid()) {
    real_t _u_m = (_u->Cell(it.Left()) + _u->Cell(it))/2;
    real_t _v_m = (_v->Cell(it.Down()) + _v->Cell(it))/2;
    _tmp->Cell(it) = sqrt(_u_m*_u_m + _v_m*_v_m);
    it.Next();
  }

  return _tmp;
}

/// Computes and returns the vorticity
const Grid* Compute::GetVorticity() {
  // Not to be done now
  return _tmp;
}

/// Computes and returns the stream line values
const Grid* Compute::GetStream() {
  // Not to be done now
  return _tmp;
}

/* Private functions
*/
/// Compute the 'new' velocites u, v
// @param dt timestep width
void Compute::NewVelocities(const real_t& dt) {
  // Initialize interior Iterator
  InteriorIterator intit(_geom);

  // Cycle through all inner cells
  while (intit.Valid()) {
    // Read access to temporary velocities F, G
    const real_t F = _F->Cell(intit);
    const real_t G = _G->Cell(intit);

    // Calculate 'new' velocities
    _u->Cell(intit) = F - dt*_p->dx_r(intit);
    _v->Cell(intit) = G - dt*_p->dy_t(intit);

    // Next cell
    intit.Next();
  }
}

/// Compute the temporary velocites F, G
// @param dt timestep width
void Compute::MomentumEqu(const real_t& dt) {
  // Initialize interior Iterator
  InteriorIterator intit(_geom);

  // Cycle through all inner cells
  while (intit.Valid()) {
    // read access to u, v
    const real_t u = _u->Cell(intit);
    const real_t v = _v->Cell(intit);

    // Parameter for Donor-Cell
    const real_t Re_inv = 1.0/_param->Re();
    const real_t alpha  = _param->Alpha();

    // Update correlation, see lecture
    _F->Cell(intit) = u + dt*( Re_inv*(_u->dxx(intit) + _u->dyy(intit)) -
      _u->DC_udu_x(intit, alpha) - _u->DC_vdu_y(intit, alpha, _v));
    _G->Cell(intit) = v + dt*( Re_inv*(_v->dxx(intit) + _v->dyy(intit)) -
      _v->DC_vdv_y(intit, alpha) - _v->DC_udv_x(intit, alpha, _u));

    // Next cell
    intit.Next();
  }  
  // Boundary update for new values of F, G
  _geom->Update_U(_F);
  _geom->Update_V(_G);
}

/// Compute the RHS of the Poisson equation
// @param dt timestep width
void Compute::RHS(const real_t& dt) {
  // Initialize interior Iterator
  InteriorIterator intit(_geom);

  // Cycle through all inner cells
  while (intit.Valid()) {
    // Update correlation for RHS, see lecture
    _rhs->Cell(intit) = (_F->dx_l(intit) + _G->dy_d(intit))/dt;

    // Next cell
    intit.Next();
  }
}
//------------------------------------------------------------------------------