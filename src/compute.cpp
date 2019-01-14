#include "compute.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
#include "grid.hpp"
#include "solver.hpp"
#include "iterator.hpp"
#include "zeitgeist.hpp"

#include <iostream>
#include <cmath>

using namespace std;
//------------------------------------------------------------------------------
/// Creates a compute instance with given geometry and parameter
//  @param geom  given geometry
//  @param param given parameter data
Compute::Compute(const Geometry* geom, const Parameter* param)
  : _geom(geom), _param(param) {
  // Initialize solver
  _solver = new SOR(geom, param->Omega());
  // Set timestep data
  _t = 0.0;
  // _dtlimit  = param->Dt_Val(); // usage?
  // Tolerance for Poisson solver
  _epslimit = param->Eps(); //evtl. anpassen Musterlösung 01
  // Offsets for the grid variables
  multi_real_t off_u;
  multi_real_t off_v;
  multi_real_t off_p;
  off_u[0] = geom->Mesh()[0];
  off_u[1] = geom->Mesh()[1]/2.0;
  off_v[0] = geom->Mesh()[0]/2.0;
  off_v[1] = geom->Mesh()[1];
  off_p[0] = geom->Mesh()[0]/2.0;
  off_p[1] = geom->Mesh()[1]/2.0;
  // Instantiate grids with offsets
  _u   = new Grid(geom, off_u);
  _v   = new Grid(geom, off_v);
  _p   = new Grid(geom, off_p);
  _T   = new Grid(geom, off_p);
  _F   = new Grid(geom, off_u);
  _G   = new Grid(geom, off_v);
  _rhs = new Grid(geom, off_p);
  _tmp = new Grid(geom, off_p);
  // Set initial values for physical variables
  _u->Initialize(geom->Velocity()[0]);
  _v->Initialize(geom->Velocity()[1]);
  _p->Initialize(geom->Pressure());
  _T->Initialize(geom->Temperature());
  _F->Initialize(geom->Velocity()[0]);
  _G->Initialize(geom->Velocity()[1]);
  _rhs->Initialize(0.0);
  _tmp->Initialize(0.0);
  _geom->Update_U(_u, _param->u_D());
  _geom->Update_V(_v, _param->v_D());
  _geom->Update_T(_T, _param->T_H(), _param->T_C(), _param->K_S());
  this->Comp_TimeStep(0.0);
}
//------------------------------------------------------------------------------
/// Deletes all grids
Compute::~Compute() {
  delete[] _u;
  delete[] _v;
  delete[] _p;
  delete[] _T;
  delete[] _F;
  delete[] _G;
  delete[] _rhs;
  delete[] _tmp;
  delete _solver;
}

//-------------------------------------------------------------------

void Compute::Comp_TimeStep(const real_t &dt_param) {
  // If explicit timestep > 0 is given, use this instead
  if (dt_param > 0) {
    _dt = dt_param;
  } else {
      cout << "Hier." << endl;
  // Find timestep as a minimum of different criteria --> first timestep without tau
  real_t dt_1 = _param->Re()*(_geom->Mesh()[0]*_geom->Mesh()[0]*_geom->Mesh()[1]*_geom->Mesh()[1])/
    (2.0*(_geom->Mesh()[0]*_geom->Mesh()[0] + _geom->Mesh()[1]*_geom->Mesh()[1]));
  // Second timestep
  real_t dt_2 = _param->Tau()*fmin(_geom->Mesh()[0], _geom->Mesh()[1])/
    fmax(_u->AbsMax(), _v->AbsMax());
  // Use minimum of dt_1 and dt_2
  _dt = min(dt_1, dt_2);
  if (_param->Pr() != 0 && _param->Pr() < 1) {
    // Third timestep
    real_t dt_3 = dt_1*_param->Pr(); // if pr > 1, this criterium is never used?!
    // Use minimum of dt and dt_3
    _dt = min(_dt, dt_3);
    }
    cout << _dt << endl;
  }
}
//------------------------------------------------------------------------------
/// Execute one time step of the fluid simulation (with or without debug info)
// @ param printInfo print information about current solver state (residual etc.)
void Compute::TimeStep(bool printInfo) {
  // Measuring of computational times
  // ZeitGeist zg;
  // Refresh boundary values
  _geom->Update_U(_u, _param->u_D());
  _geom->Update_V(_v, _param->v_D());
  _geom->Update_T(_T, _param->T_H(), _param->T_C(), _param->K_S());
  // _geom->Update_P(_p); // not necessary here
  // Measuring of computational times
  // zg.Start();



  /* // 2.4.2) measuring computation of the time step
  if (_comm->getRank() == 0 && printInfo) {
      cout << "Computational time of the time step = "
        << zg.Step() << " µs\n" << endl;
  } */
  // Compute the new temperature values explicitly (only influence of old velocity values)
  HeatTransport(_dt);
  _geom->Update_T(_T, _param->T_H(), _param->T_C(), _param->K_S()); // possibly necessary?
  // Compute temporary velocities F, G using difference schemes
  MomentumEqu(_dt);
  // Boundary update for new values of F, G
  _geom->Update_U(_F, _param->u_D());
  _geom->Update_V(_G, _param->v_D());
  // Compute RHS of the Poisson equation
  RHS(_dt);
  // Boundary update for new values of rhs
  // _geom->Update_P(_rhs); // not necessary here
  // Find solution of the Poisson equation using a SOR solver
  real_t res = 1000000.0;
  index_t i  = 0;
  /* // 2.4.1B) harmonic mean for an average step of the iteration
  real_t time_comp_inv = 0.0;
  real_t time_comm_inv = 0.0; */
  /* // 2.4.1A) one iteration only, no mean usage
  unsigned long time_comp, time_comm; */
  /* // 2.4.2) measuring the whole iteration
  unsigned long time_comp = 0;
  unsigned long time_comm = 0; */
  while (res > _param->Eps() && i < _param->IterMax()) {
    // Measuring of computational times
    /* // only 2.4.1A)
    time_comp = 0; time_comm = 0; */
    // zg.Start();
    /* if (_comm->getSize() > 1) {
      real_t res_red = ((RedOrBlackSOR*)_solver)->RedCycle(_p, _rhs);
      // 2.4.1A), 2.4.2)    time_comp += zg.Step();
      // 2.4.1B)    time_comp_inv += 1.0/zg.Step();
      // Update boundary values for pressure
      _geom->Update_P(_p);
      // 2.4.1A), 2.4.2)    time_comm += zg.Step();
      // 2.4.1B)    time_comm_inv += 1.0/zg.Step();
      real_t res_black = ((RedOrBlackSOR*)_solver)->BlackCycle(_p, _rhs);
      // 2.4.1A), 2.4.2)    time_comp += zg.Step();
      // 2.4.1B)    time_comp_inv += 1.0/zg.Step();
      // Update boundary values for pressure
      _geom->Update_P(_p);
      real_t res_comm = res_red*res_red + res_black*res_black;
      res = sqrt(_comm->gatherSum(res_comm)/_comm->getSize());
      // 2.4.1A), 2.4.2)    time_comm += zg.Step();
      // 2.4.1B)    time_comm_inv += 1.0/zg.Step();
      /* // only 2.4.1A)
      if (_comm->getRank() == 0 && printInfo) {
        cout << "Computational time for one step = "
          << time_comp << " µs\n" << endl;
        cout << "Communication time for one step = "
          << time_comm << " µs\n" << endl;
      } */
    // }
      res = _solver->Cycle(_p, _rhs);
      // 2.4.1A), 2.4.2)    time_comp += zg.Step();
      // 2.4.1B)    time_comp_inv += 1.0/zg.Step();
      // Update boundary values for pressure
      _geom->Update_P(_p, _param->p_D());
      // 2.4.1A), 2.4.2)    time_comm += zg.Step();
      // 2.4.1B)    time_comm_inv += 1.0/zg.Step();
      /* // only 2.4.1A)
      if (_comm->getRank() == 0 && printInfo) {
        cout << "Computational time for one step = "
          << time_comp << " µs\n" << endl;
        cout << "Communication time for one step = "
          << time_comm << " µs\n" << endl;
      } */
    i++;
  }
  /* // only 2.4.1B)
  if (_comm->getRank() == 0 && printInfo) {
    cout << "Harmonic mean of computational time over used steps = "
      << (index_t) (i/time_comp_inv) << " µs\n" << endl;
    cout << "Harmonic mean of communication time over used steps = "
      << (index_t) (i/time_comm_inv) << " µs\n" << endl;
  } */
  /* // only 2.4.2)
  if (_comm->getRank() == 0 && printInfo) {
    cout << "Total computational time over used steps = "
      << time_comp << " µs\n" << endl;
    cout << "Total communication time over used steps = "
      << time_comm << " µs\n" << endl;
  } */
  // Compute 'new' velocities using the pressure
  NewVelocities(_dt);
  // (optionally) printing informations
  if (printInfo) {
    cout << "_t = " << fixed << _t << "\tdt = " << scientific << _dt << " \tres = " << res
      << " \tprogress: " << fixed << (uint32_t) 100*_t/_param->Tend() << "%\n" << endl;
  }

  // Next timestep
  _t += _dt;
}
/* Getter functions
*/
/// Returns the simulated time in total
const real_t& Compute::GetTime() const {return _t;}


/// Returns the timestep
const real_t& Compute::GetTimeStep() const{return _dt;}


/// Returns the pointer to u
const Grid* Compute::GetU() const {return _u;}

// Returns the pointer to v
const Grid* Compute::GetV() const {return _v;}

// Returns the pointer to p
const Grid* Compute::GetP() const {return _p;}

// Returns the pointer to the RHS
const Grid* Compute::GetRHS() const {return _rhs;}

// Returns the pointer to T
const Grid* Compute::GetT() const {return _T;}

/// Computes and returns the absolute velocity
const Grid* Compute::GetVelocity() {
  // Initialize full Iterator
  Iterator it(_geom);



  // Go through all cells
  while (it.Valid()) {
    real_t _u_m = (_u->Cell(it.Left()) + _u->Cell(it))/2.0;
    real_t _v_m = (_v->Cell(it.Down()) + _v->Cell(it))/2.0;
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
    if (_geom->Cell(intit).type == typeFluid) {
      if (_geom->Cell(intit.Right()).type == typeFluid) {
        // Read access to temporary velocity F
        const real_t F = _F->Cell(intit);
        // Calculate 'new' velocity
        _u->Cell(intit) = F - dt*_p->dx_r(intit);
      }
      if (_geom->Cell(intit.Top()).type == typeFluid) {
        // Read access to temporary velocity G
        const real_t G = _G->Cell(intit);
        // Calculate 'new' velocity
        _v->Cell(intit) = G - dt*_p->dy_t(intit);
      }
    }

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
    if (_geom->Cell(intit).type == typeFluid) {
      // read access to u, v
      const real_t u = _u->Cell(intit);
      const real_t v = _v->Cell(intit);

      // Parameter for Donor-Cell
      const real_t Re_inv = 1.0/_param->Re();
      const real_t alpha  = _param->Alpha();

      // Update correlation, see lecture
      if (_geom->Cell(intit.Right()).type == typeFluid) {
        // Additional term through temperature inclusion (temperature value at u position)
        real_t add_u = -_param->Gx()*_param->Beta()*
          ((_T->Cell(intit) + _T->Cell(intit.Right())/2.0));
        /* real_t add_u = -_param->Gx()*_param->Beta()*
          ((_T->Cell(intit) + _T->Cell(intit.Right())/2.0)); // Larissa */
        _F->Cell(intit) = u + dt*(Re_inv*(_u->dxx(intit) + _u->dyy(intit)) -
          _u->DC_udu_x(intit, alpha) - _u->DC_vdu_y(intit, alpha, _v) + add_u);
      }
      if (_geom->Cell(intit.Top()).type == typeFluid) {
        // Additional term through temperature inclusion (temperature value at v position)
        real_t add_v = -_param->Gy()*_param->Beta()*
          ((_T->Cell(intit) + _T->Cell(intit.Top())/2.0));
        /* real_t add_v = -_param->Gy()*_param->Beta()*
          ((_T->Cell(intit) + _T->Cell(intit.Top())/2.0)); // Larissa */
        _G->Cell(intit) = v + dt*(Re_inv*(_v->dxx(intit) + _v->dyy(intit)) -
          _v->DC_vdv_y(intit, alpha) - _v->DC_udv_x(intit, alpha, _u) + add_v);
      }
    }
    // Next cell
    intit.Next();
  }
}

/// Compute the RHS of the Poisson equation
// @param dt timestep width
void Compute::RHS(const real_t& dt) {
  // Initialize interior Iterator
  InteriorIterator intit(_geom);

  // Cycle through all inner cells
  while (intit.Valid()) {
    if (_geom->Cell(intit).type == typeFluid) {
      // Update correlation for RHS, see lecture
      _rhs->Cell(intit) = (_F->dx_l(intit) + _G->dy_d(intit))/dt;
    }

    // Next cell
    intit.Next();
  }
}

/// Compute the new temperature values with the heat transport equation (using old velocity values)
// @param dt timestep width
void Compute::HeatTransport(const real_t& dt) {
  // Initialize interior Iterator
  InteriorIterator intit(_geom);

  // Cycle through all inner cells
  while (intit.Valid()) {
    if (_geom->Cell(intit).type == typeFluid && (_param->Pr() != 0)) {
    // if (_param->Pr() != 0) {
      // read access to T
      const real_t T = _T->Cell(intit);

      // Parameter for Donor-Cell
      const real_t Re_Pr_inv = 1.0/(_param->Re()*_param->Pr());
      const real_t alpha  = _param->Alpha();

      // Update correlation, see lecture
      _T->Cell(intit) = T + dt*(Re_Pr_inv*(_T->dxx(intit) + _T->dyy(intit)) -
        _T->DC_udT_x(intit, alpha, _u) - _T->DC_vdT_y(intit, alpha, _v));
    }
    // Next cell
    intit.Next();
  }
}
//------------------------------------------------------------------------------
