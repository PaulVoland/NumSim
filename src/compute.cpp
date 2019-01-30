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
Compute::Compute(Geometry* geom, const Parameter* param)
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
  _u     = new Grid(geom, off_u);
  _u_alt = new Grid(geom, off_u);
  _v     = new Grid(geom, off_v);
  _v_alt = new Grid(geom, off_u);
  _p     = new Grid(geom, off_p);
  _T     = new Grid(geom, off_p);
  _F     = new Grid(geom, off_u);
  _G     = new Grid(geom, off_v);
  _rhs   = new Grid(geom, off_p);
  _tmp   = new Grid(geom, off_p);
  // Set initial values for physical variables
  _u->Initialize(geom->Velocity()[0]);
  _u_alt->Initialize(geom->Velocity()[0]);
  _v->Initialize(geom->Velocity()[1]);
  _v_alt->Initialize(geom->Velocity()[1]);
  _p->Initialize(geom->Pressure());
  _T->Initialize(geom->Temperature());
  _F->Initialize(geom->Velocity()[0]);
  _G->Initialize(geom->Velocity()[1]);
  _rhs->Initialize(0.0);
  _tmp->Initialize(0.0);

  // Instantiate particle counter field and set to zero
  _ppc = new index_t[geom->TotalSize()[0]*geom->TotalSize()[1]];
  for (int i = 0; i < geom->TotalSize()[0]*geom->TotalSize()[1]; i++) {
    _ppc[i] = 0;
  }

  // Set particles initially to the _part_trace vector
  SetParticles();

  // Find timestep as a minimum of different criteria --> first timestep without tau
  real_t dt_1 = _param->Re()*(_geom->Mesh()[0]*_geom->Mesh()[0]*_geom->Mesh()[1]*_geom->Mesh()[1])/
    (2.0*(_geom->Mesh()[0]*_geom->Mesh()[0] + _geom->Mesh()[1]*_geom->Mesh()[1]));
  // Second timestep
  real_t dt_2 = _param->Tau()*fmin(_geom->Mesh()[0], _geom->Mesh()[1])/
    fmax(_u->AbsMax(), _v->AbsMax());
  // Use minimum of dt_1 and dt_2
  real_t _dt = min(dt_1, dt_2);
  if (_param->Pr() != 0 && _param->Pr() < 1) {
    // Third timestep
    real_t dt_3 = dt_1*_param->Pr(); // if pr > 1, this criterium is never used?!
    // Use minimum of dt and dt_3
    _dt = min(_dt, dt_3);
  }
  // If explicit timestep > 0 is given in the paramater file, use this instead
  if (_param->Dt() > 0) {
    _dt = _param->Dt();
  }
}
//------------------------------------------------------------------------------
/// Deletes all grids
Compute::~Compute() {
  delete[] _u;
  delete[] _u_alt;
  delete[] _v;
  delete[] _v_alt;
  delete[] _p;
  delete[] _T;
  delete[] _F;
  delete[] _G;
  delete[] _rhs;
  delete[] _tmp;
  delete _solver;
  //delete[] _ppc;
}
//------------------------------------------------------------------------------
/// Execute one time step of the fluid simulation (with or without debug info)
// @ param printInfo print information about current solver state (residual etc.)
void Compute::TimeStep(bool printInfo) {
  // Measuring of computational times
  // ZeitGeist zg;
  // Refresh boundary values
  _geom->Update_U(_u, _v, _param->u_D(), _dt, _param->Gx());
  _geom->Update_V(_v, _u, _param->v_D(), _dt, _param->Gy());
  _geom->Update_T(_T, _param->T_H(), _param->T_C());

  // copy velocities of the old timestep to _._alt
  CopyVelocities();

  // _geom->Update_P(_p); // not necessary here
  // Measuring of computational times
  // zg.Start();
  // Find timestep as a minimum of different criteria --> first timestep without tau
  real_t dt_1 = _param->Re()*(_geom->Mesh()[0]*_geom->Mesh()[0]*_geom->Mesh()[1]*_geom->Mesh()[1])/
    (2.0*(_geom->Mesh()[0]*_geom->Mesh()[0] + _geom->Mesh()[1]*_geom->Mesh()[1]));
  // Second timestep
  real_t dt_2 = _param->Tau()*fmin(_geom->Mesh()[0], _geom->Mesh()[1])/
    fmax(_u->AbsMax(), _v->AbsMax());
  // Use minimum of dt_1 and dt_2
  real_t _dt = min(dt_1, dt_2);
  if (_param->Pr() != 0 && _param->Pr() < 1) {
    // Third timestep
    real_t dt_3 = dt_1*_param->Pr(); // if pr > 1, this criterium is never used?!
    // Use minimum of dt and dt_3
    _dt = min(_dt, dt_3);
  }
  // If explicit timestep > 0 is given in the paramater file, use this instead
  if (_param->Dt() > 0) {
    _dt = _param->Dt();
  }
  /* // 2.4.2) measuring computation of the time step
  if (_comm->getRank() == 0 && printInfo) {
      cout << "Computational time of the time step = "
        << zg.Step() << " µs\n" << endl;
  } */
  // Compute the new temperature values explicitly (only influence of old velocity values)
  HeatTransport(_dt);
  _geom->Update_T(_T, _param->T_H(), _param->T_C()); // possibly necessary?
  // Compute temporary velocities F, G using difference schemes
  MomentumEqu(_dt);
  // Boundary update for new values of F, G
  _geom->Update_U(_F, _G, _param->u_D(), _dt, _param->Gx());
  _geom->Update_V(_G, _F, _param->v_D(), _dt, _param->Gy());
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
      _geom->Update_P(_p, _u, _v, _param->p_D(), _param->Re());
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

  //real_t text1 = PhysToVelocity(1.0 , 2.0, 'U');
  //real_t text2 = PhysToVelocity(1.0 , 2.0, 'V');
  //cout << text1 << " | " << text2  << endl;

  // Particle trace per timestep
  ParticleTrace(_dt);


  // Update cell type list for the new particle constellation
  _geom->DynamicNeighbourhood();

  // (optionally) printing informations
  if (printInfo) {
    cout << "_t = " << fixed << _t << "\tdt = " << scientific << _dt << " \tres = " << res
      << " \tprogress: " << fixed << (uint32_t) 100*_t/_param->Tend() << "%\n" << endl;
  }

  // Next timestep
  _t += _dt;
}
//------------------------------------------------------------------------------
/* Getter functions
*/
//------------------------------------------------------------------------------
/// Returns the simulated time in total
const real_t& Compute::GetTime() const {return _t;}
//------------------------------------------------------------------------------
/// Returns the pointer to u
const Grid* Compute::GetU() const {return _u;}
//------------------------------------------------------------------------------
// Returns the pointer to v
const Grid* Compute::GetV() const {return _v;}
//------------------------------------------------------------------------------
// Returns the pointer to p
const Grid* Compute::GetP() const {return _p;}
//------------------------------------------------------------------------------
// Returns the pointer to the RHS
const Grid* Compute::GetRHS() const {return _rhs;}
//------------------------------------------------------------------------------
// Returns the pointer to T
const Grid* Compute::GetT() const {return _T;}
//------------------------------------------------------------------------------
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
//------------------------------------------------------------------------------
/// Computes and returns the vorticity
const Grid* Compute::GetVorticity() {
  // Not to be done now
  return _tmp;
}
//------------------------------------------------------------------------------
/// Computes and returns the stream line values
const Grid* Compute::GetStream() {
  // Not to be done now
  return _tmp;
}
//------------------------------------------------------------------------------
/* Private functions
*/
//------------------------------------------------------------------------------
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
//------------------------------------------------------------------------------
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
        real_t add_u = _param->Gx()*(1.0 - _param->Beta()*
          ((_T->Cell(intit) + _T->Cell(intit.Right()))/2.0));
        /* real_t add_u = -_param->Gx()*_param->Beta()*
          ((_T->Cell(intit) + _T->Cell(intit.Right())/2.0)); // Larissa */
        _F->Cell(intit) = u + dt*(Re_inv*(_u->dxx(intit) + _u->dyy(intit)) -
          _u->DC_udu_x(intit, alpha) - _u->DC_vdu_y(intit, alpha, _v) + add_u);
      }
      if (_geom->Cell(intit.Top()).type == typeFluid) {
        // Additional term through temperature inclusion (temperature value at v position)
        real_t add_v = _param->Gy()*(1.0 - _param->Beta()*
          ((_T->Cell(intit) + _T->Cell(intit.Top()))/2.0));
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
//------------------------------------------------------------------------------
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
//------------------------------------------------------------------------------
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
// Set particles initially to the _part_trace vector
void Compute::SetParticles(){
  Iterator it_pc(_geom);
  index_t s = 0;
  //index_t an_inflow = 4; // should be even?!
  //index_t an_inflow = 3; // not used here
  index_t an_cell = 9;
  //it_pc.First();
  while (it_pc.Valid()) {
    if (_geom->Cell(it_pc).type == typeFluid) {
      //cout << "Ref. ##### x: " << (it_pc.Pos()[0]-1)*_geom->TotalLength()[0]/(_geom->TotalSize()[0]-2) << " y: " << (it_pc.Pos()[1]-1)*_geom->TotalLength()[1]/(_geom->TotalSize()[1]-2) << endl;
      for (int i = 0; i < an_cell; i++) {
          real_t *foo;
          foo = new real_t[2];
          _part_trace.push_back(foo);
          //cout << "hiervor x: " << _part_trace[s][0] << " y: " << _part_trace[s][1] << " s="<< s << endl;
          _part_trace[s][0] = RandNumb(it_pc.Pos()[0],it_pc.Pos()[0] - 1)*_geom->Mesh()[0];
          _part_trace[s][1] = RandNumb(it_pc.Pos()[1],it_pc.Pos()[1] - 1)*_geom->Mesh()[1];
          //cout << "Hiernach x: " << _part_trace[s][0] << " y: " << _part_trace[s][1] << " s="<< s << endl;
          //cout << "x: " << RandNumb(it_pc.Pos()[0],it_pc.Pos()[0]-1)*_geom->TotalLength()[0]/(_geom->TotalSize()[0]-2) << " y: " << RandNumb(it_pc.Pos()[1],it_pc.Pos()[1]-1)*_geom->TotalLength()[1]/(_geom->TotalSize()[1]-2) << " s="<< s << endl;
          s++;
        }
      }
    it_pc.Next();
  }
  // Set new particles in each timestep at the inflow boundaries
  SetNewInflowParticles();
  //int i = 0;
  /*for(std::vector<double*>::iterator it = _part_trace.begin(); it != _part_trace.end(); ++it) {
      cout << "Integer :" << i << " with Value 1: " <<  _part_trace[i][0] << "with Value 2:"<< _part_trace[i][1] << "\n ";
      i++;
    }*/
  //cout << " x=  " <<  _part_trace[15777][0] << " y: "<< _part_trace[15777][1];
  //cout << " x=  " <<  _part_trace[15778][0] << " y: "<< _part_trace[15778][1];
  //cout << " x=  " <<  _part_trace[15777][0] << " y: "<< _part_trace[15777][1];
  // Inflow
  // zusätzliche frage nach nord süd und west
}
//------------------------------------------------------------------------------
// Random number between max and min
real_t Compute::RandNumb(const real_t& max, const real_t& min) const {
  return ((real_t)(rand()/RAND_MAX))*(max - min) + min;
}
//------------------------------------------------------------------------------
// set new particles in each timestep at the inflow boundaries
void Compute::SetNewInflowParticles() {
  Iterator it_pc(_geom);
  index_t s = 0;
  //index_t an_inflow = 4; // should be even?!
  index_t an_inflow = 3;
  //index_t an_cell = 9; // not used here
  //it_pc.First();
    while (it_pc.Valid()) {
      if (_geom->Cell(it_pc).type == typeIn || _geom->Cell(it_pc).type == typeInH || _geom->Cell(it_pc).type == typeInV) {
        switch (_geom->Cell(it_pc).neighbour) {
        case cellN:
          for (int i = 0; i < an_inflow; i++) {
            real_t *foo;
            foo = new real_t[2];
            _part_trace.push_back(foo);
            _part_trace[s][0] = RandNumb(it_pc.Pos()[0], it_pc.Pos()[0] - 1)*_geom->Mesh()[0];
            _part_trace[s][1] = it_pc.Pos()[1]*_geom->Mesh()[1];
            s++;
          }
          break;
        case cellW:
          for (int i = 0; i < an_inflow; i++) {
            real_t *foo;
            foo = new real_t[2];
            _part_trace.push_back(foo);
            _part_trace[s][0] = (it_pc.Pos()[0] - 1)*_geom->Mesh()[0];
            _part_trace[s][1] = RandNumb(it_pc.Pos()[1], it_pc.Pos()[1] - 1)*_geom->Mesh()[1];
            s++;
          }
          break;
        case cellNW:
          // for (int i = 0; i < an_inflow; i++) {
          for (int i = 0; i < (index_t)(an_inflow/2.0); i++) {
            real_t * foo;
            foo = new real_t[2];
            _part_trace.push_back(foo);
            _part_trace[s][0] = RandNumb(it_pc.Pos()[0], it_pc.Pos()[0] - 1)*_geom->Mesh()[0];
            _part_trace[s][1] = it_pc.Pos()[1]*_geom->Mesh()[1];
            s++;
          }
          // for (int i = 0; i < an_inflow; i++) {
          for (int i = 0; i < (index_t)(an_inflow/2.0); i++) {
            real_t * foo;
            foo = new real_t[2];
            _part_trace.push_back(foo);
            _part_trace[s][0] = (it_pc.Pos()[0] - 1)*_geom->Mesh()[0];
            _part_trace[s][1] = RandNumb(it_pc.Pos()[1], it_pc.Pos()[1] - 1)*_geom->Mesh()[1];
            s++;
          }
          break;
        case cellS:
          for (int i = 0; i < an_inflow; i++) {
            real_t *foo;
            foo = new real_t[2];
            _part_trace.push_back(foo);
            _part_trace[s][0] = RandNumb(it_pc.Pos()[0], it_pc.Pos()[0] - 1)*_geom->Mesh()[0];
            _part_trace[s][1] = (it_pc.Pos()[1] - 1)*_geom->Mesh()[1];
            s++;
          }
          break;
        case cellSW:
          // for (int i = 0; i < an_inflow; i++) {
          for (int i = 0; i < (index_t)(an_inflow/2.0); i++) {
            real_t *foo;
            foo = new real_t[2];
            _part_trace.push_back(foo);
            _part_trace[s][0] = RandNumb(it_pc.Pos()[0], it_pc.Pos()[0] - 1)*_geom->Mesh()[0];
            _part_trace[s][1] = (it_pc.Pos()[1] - 1)*_geom->Mesh()[1];
            s++;
          }
          // for (int i = 0; i < an_inflow; i++) {
          for (int i = 0; i < (index_t)(an_inflow/2.0); i++) {
            real_t *foo;
            foo = new real_t[2];
            _part_trace.push_back(foo);
            _part_trace[s][0] = (it_pc.Pos()[0] - 1)*_geom->Mesh()[0];
            _part_trace[s][1] = RandNumb(it_pc.Pos()[1], it_pc.Pos()[1] - 1)*_geom->Mesh()[1];
            s++;
          }
          break;
        case cellE:
          for (int i = 0; i < an_inflow; i++) {
            real_t *foo;
            foo = new real_t[2];
            _part_trace.push_back(foo);
            _part_trace[s][0] = it_pc.Pos()[0]*_geom->Mesh()[0];
            _part_trace[s][1] = RandNumb(it_pc.Pos()[1], it_pc.Pos()[1] - 1)*_geom->Mesh()[1];
            s++;
          }
          break;
        case cellNE:
          // for (int i = 0; i < an_inflow; i++) {
          for (int i = 0; i < (index_t)(an_inflow/2.0); i++) {
            real_t * foo;
            foo = new real_t[2];
            _part_trace.push_back(foo);
            _part_trace[s][0] = RandNumb(it_pc.Pos()[0], it_pc.Pos()[0] - 1)*_geom->Mesh()[0];
            _part_trace[s][1] = it_pc.Pos()[1]*_geom->Mesh()[1];
            s++;
          }
          // for (int i = 0; i < an_inflow; i++) {
          for (int i = 0; i < (index_t)(an_inflow/2.0); i++) {
            real_t *foo;
            foo = new real_t[2];
            _part_trace.push_back(foo);
            _part_trace[s][0] = it_pc.Pos()[0]*_geom->Mesh()[0];
            _part_trace[s][1] = RandNumb(it_pc.Pos()[1], it_pc.Pos()[1] - 1)*_geom->Mesh()[1];
            s++;
          }
          break;
        case cellSE:
          // for (int i = 0; i < an_inflow; i++) {
          for (int i = 0; i < (index_t)(an_inflow/2.0); i++) {
            real_t *foo;
            foo = new real_t[2];
            _part_trace.push_back(foo);
            _part_trace[s][0] = RandNumb(it_pc.Pos()[0], it_pc.Pos()[0] - 1)*_geom->Mesh()[0];
            _part_trace[s][1] = (it_pc.Pos()[1] - 1)*_geom->Mesh()[1];
            s++;
          }
          // for (int i = 0; i < an_inflow; i++) {
          for (int i = 0; i < (index_t)(an_inflow/2.0); i++) {
            real_t *foo;
            foo = new real_t[2];
            _part_trace.push_back(foo);
            _part_trace[s][0] = it_pc.Pos()[0]*_geom->Mesh()[0];
            _part_trace[s][1] = RandNumb(it_pc.Pos()[1], it_pc.Pos()[1] - 1)*_geom->Mesh()[1];
            s++;
          }
          break;
        default:
          break;
        };
      }
    it_pc.Next();
  }
}
//------------------------------------------------------------------------------
// Particle trace per timestep
void Compute::ParticleTrace(const real_t &dt) {
  index_t _increm_x = _geom->TotalSize()[1];
  index_t _increm_y = _geom->TotalSize()[0];
  index_t _num_cell = _increm_x*_increm_y;
  real_t vel_u_old = 0;
  real_t vel_u_new = 0;
  real_t vel_v_old = 0;
  real_t vel_v_new = 0;
  multi_index_t index_pos;
  index_t cell_number;
  index_t part_crit = 1; // particle criterion to set a cell to fluid cell

  for (int i=0; i < _num_cell; i++) {
    _ppc[i] = 0;
  }

  // Leap Frog for prediction of particle positions
  int i= 0;
  for (vec_arr::iterator it = _part_trace.begin(); it != _part_trace.end(); it++) {
      // Get all necessary velocities to propagate the recent particle
      vel_u_old = PhysToVelocity(_part_trace[i][0], _part_trace[i][1], 'u');
      vel_u_new = PhysToVelocity(_part_trace[i][0], _part_trace[i][1], 'U');
      vel_v_old = PhysToVelocity(_part_trace[i][0], _part_trace[i][1], 'v');
      vel_v_new = PhysToVelocity(_part_trace[i][0], _part_trace[i][1], 'V');

      // Calculate next position via Leap Frog with mean of old and new velocity
      _part_trace[i][0] += dt*(vel_u_old + vel_u_new)/2.0;
      _part_trace[i][1] += dt*(vel_v_old + vel_v_new)/2.0;
      // Calculate the corresponding index of the physical coordinates
      index_pos = PhysToIndex(_part_trace[i][0], _part_trace[i][1]);
      //cout << "New x=" << index_pos[0] << " y=" << index_pos[1] << endl;
      //cout << "Size x=" << _increm_x << " y=" << _increm_y << endl;
      //cout << "New x=" << _part_trace[i][0] << " y=" << _part_trace[i][1] <<" Cell Number " << cell_number << endl;
      // sort Partical in to the right cellnumber or delete
      //cout << _num_cell << endl;
      // if (_part_trace[i][0] < 0 || _part_trace[i][0] > _geom->TotalLength()[0]
      //   || _part_trace[i][1] < 0 || _part_trace[i][1] > _geom->TotalLength()[1]) {
      if ( 1 > index_pos[0] || index_pos[0] >= _increm_y || 1 > index_pos[1] || index_pos[1] >= _increm_x) {
        //cout << "Raus x=" << index_pos[0] << " y=" << index_pos[1] << endl;
        _part_trace.erase(it);
        it--;
      } else {
        // Calculate the fitting cell number
        cell_number = IndexToCell(index_pos);
        //cout << cell_number << endl;
        // #####_ppc[cell_number] = _ppc[cell_number]+1;
        // New Interator at position cell_number
        Iterator it_cell(_geom, cell_number);
        // Set the type of the cell
        if (_geom->Cell(it_cell).type == typeFluid || _geom->Cell(it_cell).type == typeEmpty
          || _geom->Cell(it_cell).type == typeSurf) { // particle in inner cell where at least there is no obstacle
          _ppc[cell_number] += 1;
          if (part_crit <= _ppc[cell_number]) {
            _geom->SetCell(it_cell).type = typeFluid;
          } else {
            _geom->SetCell(it_cell).type = typeEmpty;
          }
        } else { // particle in inner obstacle cell
          _part_trace.erase(it);
          it--;
          i--;
        }
        i++;
      }
    }
  // ############################ Hier zu Debugzwecken #############################
  string zeile;
  string spalte;
  for (int i=0; i < _num_cell; i++) {
    if (i % _increm_y == 0) {
      spalte = zeile + "\n" + spalte;
      zeile = "";
      // only formatting
      if (_ppc[i] > 99) {
        zeile = zeile + "|" + to_string(_ppc[i]);
      } else if (_ppc[i] > 9) {
        zeile = zeile + "|" + " " + to_string(_ppc[i]);
      } else {
        zeile = zeile + "|" + "  " + to_string(_ppc[i]);
      }
    } else {
      // only formatting
      if (_ppc[i] > 99) {
        zeile = zeile + "|" + to_string(_ppc[i]);
      } else if (_ppc[i] > 9) {
        zeile = zeile + "|" + " " + to_string(_ppc[i]);
      } else {
        zeile = zeile + "|" + "  " + to_string(_ppc[i]);
      }
    }
  }
  spalte = zeile + "\n" + spalte;
  cout << spalte << endl;
  // ###############################################################################
  // Set new particles in each timestep at the inflow boundaries
  SetNewInflowParticles();
}
//------------------------------------------------------------------------------
// Calculate the index positions from physical coordinates
multi_index_t Compute::PhysToIndex(const real_t& x , const real_t& y) const {
  multi_real_t h    = _geom->Mesh();
  // Instantiate indices and distances
  index_t i, j;
  multi_index_t value;
  //cout << "Variablen in PhysToIndex x= " << x << " y= " << y <<endl;
  //cout << "Variablen in PhysToIndex i= " << i << " j= " << j <<endl;
  // find inner cell index for anchor cell in format
  // h(0)*[i,i+1) x h(1)*[j,j+1) (i,j = 0,...,_geom->TotalSize()[0,1])-1) (inner numbering)
  // if x,y >= 0
  if (x < 0) {
    i = 0; // is in outer index format
  } else {
    i = (index_t)(x/h[0]); // is in inner index format
    i++; // convert to outer index format
  }
  if (y < 0) {
    j = 0; // is in outer index format
  } else {
    j = (index_t)(y/h[1]); // is in inner index format
    j++; // convert to outer index format
  }
  value[0] = i;
  value[1] = j;
  //cout << "Variablen in PhysToIndex i= " << i << " j= " << j <<endl;
  //cout << "Variablen in PhysToIndex i= " << value[0] << " j= " << value[1] <<endl;
  return value;
}
//------------------------------------------------------------------------------
// Return the index position of a multi_index-type
index_t Compute::IndexToCell(const multi_index_t& value) const {
  index_t x = value[0];
  index_t y = value[1];
  index_t _increm_y = _geom->TotalSize()[0];
  //benutze version aus interpolate
  //cout << "x = " << x << " y=" << y<< endl;
  //cout << "TotalSize x=" << _geom->TotalSize()[0]-2 << " y=" << _geom->TotalSize()[1]-2 <<" TotalLength x=" << _geom->TotalLength()[0]<< " y=" << _geom->TotalLength()[1] << endl;
  return x + y*_increm_y;
}
//------------------------------------------------------------------------------
// Get special velocity defined by char f in {u, U, v, V} at physical position
real_t Compute::PhysToVelocity(const real_t& x , const real_t& y ,const char& f) const {
  multi_real_t pos;
  pos[0] = x;
  pos[1] = y;
  real_t value;
  if (f == 'u' || f == 'U' || f == 'v' || f == 'V' ) {
    if (f == 'u') {
      value = _u_alt->Interpolate(pos);
      //cout <<"value" << value << endl;
    } else if (f == 'U') {
      value = _u->Interpolate(pos);
      //cout <<"value" << value << endl;
    } else if (f == 'v') {
      value = _v_alt->Interpolate(pos);
      //cout <<"value" << value << endl;
    } else if (f == 'V') {
      value = _v->Interpolate(pos);
      //cout <<"value" << value << endl;
    }
  }
  return value;
}
//------------------------------------------------------------------------------
// Copy velocities of the old timestep to _._alt
void Compute::CopyVelocities() {
  Iterator it(_geom);
  while (it.Valid()) {
    _u_alt->Cell(it) = _u->Cell(it);
    _v_alt->Cell(it) = _v->Cell(it);
    it.Next();
  }
}
//------------------------------------------------------------------------------
