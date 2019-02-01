#include "compute.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
#include "grid.hpp"
#include "solver.hpp"
#include "iterator.hpp"
#include "zeitgeist.hpp"
#include "typedef.hpp"


#include <iostream>
#include <cmath>
#include <string>

using namespace std;

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"


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
  for (index_t i = 0; i < geom->TotalSize()[0]*geom->TotalSize()[1]; i++) {
    _ppc[i] = 0;
  }

  // set partical trace array with fluid cells
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
  delete[] _v;
  delete[] _u_alt;
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
  //ShowVelocitysPressures('U');
  //ShowVelocitysPressures('V');
  // Refresh boundary values
  _geom->Update_U(_u, _v, _param->u_D(), _dt, _param->Gx());
  //ShowVelocitysPressures('U');
  //ShowVelocitysPressures('V');
  _geom->Update_V(_v, _u, _param->v_D(), _dt, _param->Gy());
  //ShowVelocitysPressures('U');
  _geom->Update_P(_p, _u, _v, _param->p_D(), _param->Re());
  //ShowVelocitysPressures('U');
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
  //ShowVelocitysPressures('F');
  // Boundary update for new values of F, G
  _geom->Update_U(_F, _G, _param->u_D(), _dt, _param->Gx());
  //ShowVelocitysPressures('F');
  _geom->Update_V(_G, _F, _param->v_D(), _dt, _param->Gy());
  //ShowVelocitysPressures('F');
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
      //ShowVelocitysPressures('U');
      // Update boundary values for pressure
      _geom->Update_P(_p, _u, _v, _param->p_D(), _param->Re());
      //ShowVelocitysPressures('U');
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
  //ShowVelocitysPressures('U');
  _geom->Update_U(_u, _v, _param->u_D(), _dt, _param->Gx());
  //ShowVelocitysPressures('U');
  _geom->Update_V(_v, _u, _param->v_D(), _dt, _param->Gy());
  //ShowVelocitysPressures('U');

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
    //if (_geom->Cell(intit).type == typeFluid || _geom->Cell(intit).type == typeSurf) {
    if (_geom->Cell(intit).type == typeFluid) {
      if (_geom->Cell(intit.Right()).type == typeFluid || _geom->Cell(intit.Right()).type == typeSurf) {
        // Read access to temporary velocity F
        const real_t F = _F->Cell(intit);
        // Calculate 'new' velocity
        _u->Cell(intit) = F - dt*_p->dx_r(intit);
      }
      if (_geom->Cell(intit.Top()).type == typeFluid || _geom->Cell(intit.Top()).type == typeSurf) {
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
    //if (_geom->Cell(intit).type == typeFluid || _geom->Cell(intit).type == typeSurf) {
    if (_geom->Cell(intit).type == typeFluid) {
      // read access to u, v
      const real_t u = _u->Cell(intit);
      const real_t v = _v->Cell(intit);

      // Parameter for Donor-Cell
      const real_t Re_inv = 1.0/_param->Re();
      const real_t alpha  = _param->Alpha();

      // Update correlation, see lecture
      if (_geom->Cell(intit.Right()).type == typeFluid || _geom->Cell(intit.Right()).type == typeSurf) {
        // Additional term through temperature inclusion (temperature value at u position)
        real_t add_u = _param->Gx()*(1.0 - _param->Beta()*
          ((_T->Cell(intit) + _T->Cell(intit.Right()))/2.0));
        /* real_t add_u = -_param->Gx()*_param->Beta()*
          ((_T->Cell(intit) + _T->Cell(intit.Right())/2.0)); // Larissa */
        _F->Cell(intit) = u + dt*(Re_inv*(_u->dxx(intit) + _u->dyy(intit)) -
          _u->DC_udu_x(intit, alpha) - _u->DC_vdu_y(intit, alpha, _v) + add_u);
      }
      if (_geom->Cell(intit.Top()).type == typeFluid || _geom->Cell(intit.Top()).type == typeSurf) {
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
    //if (_geom->Cell(intit).type == typeFluid || _geom->Cell(intit).type == typeSurf) {
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
    //if ((_geom->Cell(intit).type == typeFluid || _geom->Cell(intit).type == typeSurf) && (_param->Pr() != 0)) {
    if (_geom->Cell(intit).type == typeFluid && (_param->Pr() != 0)) {
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
void Compute::SetParticles() {
  Iterator it_pc(_geom);
  index_t s = 0;
  //index_t an_inflow = 4; // should be even?!
  //index_t an_inflow = 3; // not used here
  index_t an_cell = 9;
  //it_pc.First();
  while (it_pc.Valid()) {
    if (_geom->Cell(it_pc).type == typeFluid || _geom->Cell(it_pc).type == typeSurf) {
        //cout << "Ref. ##### x: " << (it_pc.Pos()[0]-1)*_geom->TotalLength()[0]/(_geom->TotalSize()[0]-2) << " y: " << (it_pc.Pos()[1]-1)*_geom->TotalLength()[1]/(_geom->TotalSize()[1]-2) << endl;
      for (index_t i = 0; i < an_cell; ++i)
        {
          real_t * foo;
          foo = new real_t[2];
          _part_trace.push_back(foo);
          //cout << "hiervor x: " << _part_trace[s][0] << " y: " << _part_trace[s][1] << " s="<< s << endl;
          _part_trace[s][0] = RandNumb(it_pc.Pos()[0],it_pc.Pos()[0] - 1)*_geom->Mesh()[0];
          _part_trace[s][1] = RandNumb(it_pc.Pos()[1],it_pc.Pos()[1] - 1)*_geom->Mesh()[1];
          //ShowParticleToCellDebug(_part_trace[s][0],_part_trace[s][1]);
          //cout << "x: " << RandNumb(it_pc.Pos()[0],it_pc.Pos()[0]-1)*_geom->TotalLength()[0]/(_geom->TotalSize()[0]-2) << " y: " << RandNumb(it_pc.Pos()[1],it_pc.Pos()[1]-1)*_geom->TotalLength()[1]/(_geom->TotalSize()[1]-2) << " s="<< s << endl;
          s++;
        }
      }
  it_pc.Next();
  }
  SetNewInflowParticles();
  //index_t i = 0;
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
real_t Compute::RandNumb(const real_t& max, const real_t& min) const {
  return ((real_t)(rand())/RAND_MAX)*(max - min) + min;
}
//------------------------------------------------------------------------------
void Compute::SetNewInflowParticles() {
  Iterator it_pc(_geom);
  index_t s = _part_trace.size();
  //index_t an_inflow = 4; // should be even?!
  index_t an_inflow = 3;
  //index_t an_cell = 9; // not used here
  //it_pc.First();
    while (it_pc.Valid()) {
      if (_geom->Cell(it_pc).type == typeIn || _geom->Cell(it_pc).type == typeInH || _geom->Cell(it_pc).type == typeInV) {
        switch (_geom->Cell(it_pc).neighbour) {
        case cellN:
          for (index_t i = 0; i < an_inflow; i++) {
            real_t *foo;
            foo = new real_t[2];
            _part_trace.push_back(foo);
            _part_trace[s][0] = RandNumb(it_pc.Pos()[0], it_pc.Pos()[0] - 1)*_geom->Mesh()[0];
            _part_trace[s][1] = it_pc.Pos()[1]*_geom->Mesh()[1];
            s++;
          }
          break;
        case cellW:
          for (index_t i = 0; i < an_inflow; i++) {
            real_t *foo;
            foo = new real_t[2];
            _part_trace.push_back(foo);
            _part_trace[s][0] = (it_pc.Pos()[0] - 1)*_geom->Mesh()[0];
            _part_trace[s][1] = RandNumb(it_pc.Pos()[1], it_pc.Pos()[1] - 1)*_geom->Mesh()[1];
            s++;
          }
          break;
        case cellNW:
          // for (index_t i = 0; i < an_inflow; i++) {
          for (index_t i = 0; i < (index_t)(an_inflow/2.0); i++) {
            real_t * foo;
            foo = new real_t[2];
            _part_trace.push_back(foo);
            _part_trace[s][0] = RandNumb(it_pc.Pos()[0], it_pc.Pos()[0] - 1)*_geom->Mesh()[0];
            _part_trace[s][1] = it_pc.Pos()[1]*_geom->Mesh()[1];
            s++;
          }
          // for (index_t i = 0; i < an_inflow; i++) {
          for (index_t i = 0; i < (index_t)(an_inflow/2.0); i++) {
            real_t * foo;
            foo = new real_t[2];
            _part_trace.push_back(foo);
            _part_trace[s][0] = (it_pc.Pos()[0] - 1)*_geom->Mesh()[0];
            _part_trace[s][1] = RandNumb(it_pc.Pos()[1], it_pc.Pos()[1] - 1)*_geom->Mesh()[1];
            s++;
          }
          break;
        case cellS:
          for (index_t i = 0; i < an_inflow; i++) {
            real_t *foo;
            foo = new real_t[2];
            _part_trace.push_back(foo);
            _part_trace[s][0] = RandNumb(it_pc.Pos()[0], it_pc.Pos()[0] - 1)*_geom->Mesh()[0];
            _part_trace[s][1] = (it_pc.Pos()[1] - 1)*_geom->Mesh()[1];
            s++;
          }
          break;
        case cellSW:
          // for (index_t i = 0; i < an_inflow; i++) {
          for (index_t i = 0; i < (index_t)(an_inflow/2.0); i++) {
            real_t *foo;
            foo = new real_t[2];
            _part_trace.push_back(foo);
            _part_trace[s][0] = RandNumb(it_pc.Pos()[0], it_pc.Pos()[0] - 1)*_geom->Mesh()[0];
            _part_trace[s][1] = (it_pc.Pos()[1] - 1)*_geom->Mesh()[1];
            s++;
          }
          // for (index_t i = 0; i < an_inflow; i++) {
          for (index_t i = 0; i < (index_t)(an_inflow/2.0); i++) {
            real_t *foo;
            foo = new real_t[2];
            _part_trace.push_back(foo);
            _part_trace[s][0] = (it_pc.Pos()[0] - 1)*_geom->Mesh()[0];
            _part_trace[s][1] = RandNumb(it_pc.Pos()[1], it_pc.Pos()[1] - 1)*_geom->Mesh()[1];
            s++;
          }
          break;
        case cellE:
          for (index_t i = 0; i < an_inflow; i++) {
            real_t *foo;
            foo = new real_t[2];
            _part_trace.push_back(foo);
            _part_trace[s][0] = it_pc.Pos()[0]*_geom->Mesh()[0];
            _part_trace[s][1] = RandNumb(it_pc.Pos()[1], it_pc.Pos()[1] - 1)*_geom->Mesh()[1];
            s++;
          }
          break;
        case cellNE:
          // for (index_t i = 0; i < an_inflow; i++) {
          for (index_t i = 0; i < (index_t)(an_inflow/2.0); i++) {
            real_t * foo;
            foo = new real_t[2];
            _part_trace.push_back(foo);
            _part_trace[s][0] = RandNumb(it_pc.Pos()[0], it_pc.Pos()[0] - 1)*_geom->Mesh()[0];
            _part_trace[s][1] = it_pc.Pos()[1]*_geom->Mesh()[1];
            s++;
          }
          // for (index_t i = 0; i < an_inflow; i++) {
          for (index_t i = 0; i < (index_t)(an_inflow/2.0); i++) {
            real_t *foo;
            foo = new real_t[2];
            _part_trace.push_back(foo);
            _part_trace[s][0] = it_pc.Pos()[0]*_geom->Mesh()[0];
            _part_trace[s][1] = RandNumb(it_pc.Pos()[1], it_pc.Pos()[1] - 1)*_geom->Mesh()[1];
            s++;
          }
          break;
        case cellSE:
          // for (index_t i = 0; i < an_inflow; i++) {
          for (index_t i = 0; i < (index_t)(an_inflow/2.0); i++) {
            real_t *foo;
            foo = new real_t[2];
            _part_trace.push_back(foo);
            _part_trace[s][0] = RandNumb(it_pc.Pos()[0], it_pc.Pos()[0] - 1)*_geom->Mesh()[0];
            _part_trace[s][1] = (it_pc.Pos()[1] - 1)*_geom->Mesh()[1];
            s++;
          }
          // for (index_t i = 0; i < an_inflow; i++) {
          for (index_t i = 0; i < (index_t)(an_inflow/2.0); i++) {
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
void Compute::ParticleTrace(const real_t &dt){
  index_t _increm_x = _geom->TotalSize()[1];
  index_t _increm_y = _geom->TotalSize()[0];
  index_t _num_cell = _increm_x*_increm_y;
  real_t vel_v_old = 0;
  real_t vel_v_new = 0;
  real_t vel_u_old = 0;
  real_t vel_u_new = 0;
  index_t cell_number = 0;
  index_t part_crit = 1; // criterion to set fluid cell
  multi_index_t index_pos;
  for (index_t i=0;i<_num_cell;i++)
    _ppc[i] = 0;

  
// Leap Frog
  index_t i= 0;
  for(vec_arr::iterator it = _part_trace.begin(); it != _part_trace.end(); ++it) {
      // Calculate Velocity V_i+1/2
      vel_v_old =  PhysToVelocity(_part_trace[i][0],_part_trace[i][1] , 'v');
      vel_v_new =  PhysToVelocity(_part_trace[i][0],_part_trace[i][1] , 'V');
      vel_u_old =  PhysToVelocity(_part_trace[i][0],_part_trace[i][1] , 'u');
      vel_u_new =  PhysToVelocity(_part_trace[i][0],_part_trace[i][1] , 'U');
      // Calculate next Position
      //cout << "########## Before ##################" << endl;
      //ShowParticleToCellDebug(_part_trace[i][0],_part_trace[i][1]);
      _part_trace[i][0] = _part_trace[i][0] +  dt*(vel_u_old + vel_u_new)/2.0;
      _part_trace[i][1] = _part_trace[i][1] +  dt*(vel_v_old + vel_v_new)/2.0;
      //cout << "########## After ###################" << endl;
      //ShowParticleToCellDebug(_part_trace[i][0],_part_trace[i][1]);
      // calculate the index from the phys coord.
      index_pos = PhysToIndex(_part_trace[i][0],_part_trace[i][1]);
      //cout << "New x=" << index_pos[0] << " y=" << index_pos[1] << endl;
      //cout << "Size x=" << _increm_x << " y=" << _increm_y << endl;
      //cout << "New x=" << _part_trace[i][0] << " y=" << _part_trace[i][1] <<" Cell Number " << cell_number << endl;
      // sort Partical in to the right cellnumber or delete
      //cout << _num_cell << endl;
      if ( 1 > index_pos[0] || index_pos[0] >= _increm_y || 1 > index_pos[1] || index_pos[1] >= _increm_x)// abfrage nach physikalische positionen
      {
        //cout << "Raus x=" << index_pos[0] << " y=" << index_pos[1] << endl;
        _part_trace.erase(it);
        it--;

      }else{
        // calculate the Cell Number
        cell_number = IndexToCell(index_pos);
        //cout << cell_number << endl;
        // #####_ppc[cell_number] = _ppc[cell_number]+1;
        // New Interator on cell_number
        Iterator it_cell = Iterator(_geom,cell_number);
        //  set the type of the cell
        if (_geom->Cell(it_cell).type == typeFluid || _geom->Cell(it_cell).type == typeEmpty
          || _geom->Cell(it_cell).type == typeSurf ) // hier sollte noch E und g abgefragt werden
        {
          _ppc[cell_number] = _ppc[cell_number]+1;
        } else {
          _part_trace.erase(it);
          it--;
          i--;
        }
        i++;
      }
    }
      for (index_t i = 0; i < _num_cell; i++) {
      Iterator it_cell = Iterator(_geom,i); 
        if (_geom->Cell(it_cell).type == typeFluid || _geom->Cell(it_cell).type == typeEmpty
          || _geom->Cell(it_cell).type == typeSurf ) // hier sollte noch E und g abgefragt werden
        {
          if ( part_crit <= _ppc[i])
          {
            _geom->setCell(it_cell).type = typeFluid;
            _geom->setCell(it_cell).neighbour = cellNone;
          } else{
            _geom->setCell(it_cell).type = typeEmpty;
            //cout << "Zelle " << i << " ist Leer" << endl;
            _geom->setCell(it_cell).neighbour = cellNone;
            _u->Cell(it_cell) = 0.0;
            _v->Cell(it_cell) = 0.0;
            _p->Cell(it_cell) = 0.0;
          }
        }
      }



  //################ Debug #################################
  //ShowParticle();
  //ShowNeighbour();
  //ShowType();
  //########################################################
  // new Partikel from Inflow
  SetNewInflowParticles();
}
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


index_t Compute::IndexToCell(const multi_index_t& value) const {
  index_t x = value[0];
  index_t y = value[1];
  index_t _increm_y = _geom->TotalSize()[0];
  //benutze version aus interpolate
  //cout << "x = " << x << " y=" << y<< endl;
  //cout << "TotalSize x=" << _geom->TotalSize()[0]-2 << " y=" << _geom->TotalSize()[1]-2 <<" TotalLength x=" << _geom->TotalLength()[0]<< " y=" << _geom->TotalLength()[1] << endl;
  return x + y*_increm_y;
}
real_t Compute::PhysToVelocity(const real_t& x , const real_t& y ,const char& f) const {
  multi_real_t velo;
  velo[0] = x ;
  velo[1] = y;
  real_t value = 0;
  if (f=='v' || f=='V' || f=='u' || f=='U' ) // klein ist alt groß ist neu
  {
    if (f=='u')
    {
      value = _u_alt->Interpolate(velo);
      //cout <<"value" << value << endl;
    } else if (f=='U')
    {
      value = _u->Interpolate(velo);
      //cout <<"value" << value << endl;
    }
    else if (f=='v')
    {
      value = _v_alt->Interpolate(velo);
      //cout <<"value" << value << endl;
    }
    else if (f=='V')
    {
      value = _v->Interpolate(velo);
      //cout <<"value" << value << endl;
    }
  }

  return value;
}
void Compute::CopyVelocities(){
  Iterator it = Iterator(_geom);
  while (it.Valid()){
    _u_alt->Cell(it) = _u->Cell(it);
    _v_alt->Cell(it) = _v->Cell(it);
    it.Next();
  }
}

  void Compute::ShowParticle(){
  cout << "\n";
  green("Das ist die Ausgabe für Anzahl Partikel pro Zelle \n ");
  cout << "\n";
  index_t _increm_x = _geom->TotalSize()[1];
  index_t _increm_y = _geom->TotalSize()[0];
  index_t _num_cell = _increm_x*_increm_y;
  string zeile;
  string spalte;
  //cout << _num_cell << endl;
  //cout << _increm_y << endl;
  for (index_t i=0;i<_num_cell;i++){
    if (i%(_increm_y) ==0 )
    {

      spalte = zeile + "\n" + spalte;
      zeile = "";

      if (_ppc[i] > 99)
      {
        zeile = zeile + "|" + to_string(_ppc[i]);
      } else if (_ppc[i] > 9){
        zeile = zeile + "|" + " " + to_string(_ppc[i]);
      } else {
        zeile = zeile + "|" + "  " + to_string(_ppc[i]);
      }

    } else{
      if (_ppc[i] > 99)
      {
        zeile = zeile + "|" + to_string(_ppc[i]);
      } else if (_ppc[i] > 9){
        zeile = zeile + "|" + " " + to_string(_ppc[i]);
      } else {
        zeile = zeile + "|" + "  " + to_string(_ppc[i]);
      }
    }
  }
  spalte = zeile + "\n" + spalte;
  cout << spalte;


  }
  // Show Neighbour
  void Compute::ShowNeighbour(){
  cout << "\n";
  green("Das ist die Ausgabe für den Zellnachbarn \n");
  blue("\n\
  cellN = 1,    // Cell N\n\
  cellW = 2,    // Cell W\n\
  cellNW = 3,   // Cells N & W\n\
  cellS = 4,    // Cell S\n\
  cellNS = 5,   // Cells N & S\n\
  cellSW = 6,   // Cells S & W\n\
  cellNWS = 7,  // Cells N & W & S\n\
  cellE = 8,    // Cell E\n\
  cellNE = 9,   // Cells N & E\n\
  cellWE = 10,  // Cells W & E\n\
  cellNWE = 11, // Cells N & W & E\n\
  cellSE = 12,  // Cells S & E\n\
  cellNSE = 13, // Cells N & S & E\n\
  cellWSE = 14, // Cells W & E & S\n\
  cellAll = 15  // Cells N & W & S & E\n\
  ");
  cout << "\n";
  index_t _increm_x = _geom->TotalSize()[1];
  index_t _increm_y = _geom->TotalSize()[0];
  index_t _num_cell = _increm_x*_increm_y;
  string zeile;
  string spalte;
  //cout << _num_cell << endl;
  //cout << _increm_y << endl;
  for (index_t i=0;i<_num_cell;i++){
    Iterator it(_geom, i);
    if (i%(_increm_y) ==0 )
    {

      spalte = zeile + "\n" + spalte;
      zeile = "";

      if (_geom->Cell(it).neighbour > 99)
      {
        zeile = zeile + "|" + to_string(_geom->Cell(it).neighbour);
      } else if (_geom->Cell(it).neighbour > 9){
        zeile = zeile + "|" + " " + to_string(_geom->Cell(it).neighbour);
      } else {
        zeile = zeile + "|" + "  " + to_string(_geom->Cell(it).neighbour);
      }

    } else{
      if (_geom->Cell(it).neighbour > 99)
      {
        zeile = zeile + "|" + to_string(_geom->Cell(it).neighbour);
      } else if (_geom->Cell(it).neighbour > 9){
        zeile = zeile + "|" + " " + to_string(_geom->Cell(it).neighbour);
      } else {
        zeile = zeile + "|" + "  " + to_string(_geom->Cell(it).neighbour);
      }
    }
  }
  spalte = zeile + "\n" + spalte;
  cout << spalte;

  }
  // Show Type
  void Compute::ShowType(){
  cout << "\n";
  green("Das ist die Ausgabe für den Zelltypen \n");
  blue("\n\
  0  typeFluid,   // Standard fluid cell\n\
  1  typeSolid,   // Simple wall, no slip, isolated against heat transfer\n\
  2  typeIn,      // Simple inflow (forced velocity)\n\
  3  typeInH,     // Horizontal inflow (parabolic)\n\
  4  typeInV,     // Vertical inflow (parabolic)\n\
  5  typeSlipH,   // Horizontal slip boundary\n\
  6  typeSlipV,   // Vertical slip boundary\n\
  7  typeOut,     // Outflow\n\
  8  typeTDir_h,  // Dirichlet value for higher temperature (u,v,p treated as no slip)\n\
  9  typeTDir_c,  // Dirichlet value for lower temperature (u,v,p treated as no slip)\n\
  10 typeEmpty,   // Empty cell (air)\n\
  11 typeSurf     // Inner surface cell between fluid and empty space\n\
  ");
  cout << "\n";

  index_t _increm_x = _geom->TotalSize()[1];
  index_t _increm_y = _geom->TotalSize()[0];
  index_t _num_cell = _increm_x*_increm_y;
  string zeile;
  string spalte;

  //cout << _num_cell << endl;
  //cout << _increm_y << endl;
  for (index_t i=0;i<_num_cell;i++){
    Iterator it(_geom, i);
    if (i%(_increm_y) ==0 )
    {

      spalte = zeile + "\n" + spalte;
      zeile = "";

      if (_geom->Cell(it).type > 99)
      {
        zeile = zeile + "|" + to_string_color(_geom->Cell(it).type);
      } else if (_geom->Cell(it).type > 9){
        zeile = zeile + "|" + " " + to_string_color(_geom->Cell(it).type);
      } else {
        zeile = zeile + "|" + "  " + to_string_color(_geom->Cell(it).type);
      }

    } else{
      if (_geom->Cell(it).type > 99)
      {
        zeile = zeile + "|" + to_string_color(_geom->Cell(it).type);
      } else if (_geom->Cell(it).type > 9){
        zeile = zeile + "|" + " " + to_string_color(_geom->Cell(it).type);
      } else {
        zeile = zeile + "|" + "  " + to_string_color(_geom->Cell(it).type);
      }
    }
  }
  spalte = zeile + "\n" + spalte;
  cout << spalte;

  }

  void Compute::ShowParticleToCellDebug(const real_t &x , const real_t &y){
    multi_index_t versuch;
    versuch = PhysToIndex(x,y);
    cout << "Phys: x= " << x << " y= " << y << " Index x: " << versuch[0] << " y: " << versuch[1] <<" Cell: "<<IndexToCell(versuch) << endl;
  }
  // Options 'v' 'u' for old velocities, 'V' 'U' new Velocitys  
  // 'F' 'G' 'p' 'T' for the other dynamic parameter
  void Compute::ShowVelocitysPressures(char f){
    cout << "\n";
    Iterator it(_geom);
    // Create Boundary Iterator instances for debug purposes
    BoundaryIterator bit0(_geom);
    BoundaryIterator bit1(_geom);
    BoundaryIterator bit2(_geom);
    BoundaryIterator bit3(_geom);
    // Create Interior Iterator instance for debug purposes
    //InteriorIterator intit(_geom);
    InteriorIterator init(_geom);

    if (f=='v' || f=='V' || f=='u' || f=='U' || f=='p' || f=='F' || f=='G' || f=='T' ) // klein ist alt groß ist neu
  {
    if (f=='u')
    {
    green("Das ist die Ausgabe für _u_alt \n"); 
    cout << "\n"; 
    string writelocation = "";
    string writevalue = "";
    bit0.SetBoundary(2);
    bit0.First();
    while (bit0.Valid()){
      writelocation = writelocation + "(" + to_string(bit0.Pos()[0])+ "," ; 
      writelocation = writelocation + to_string(bit0.Pos()[1]) + ")   |  ";
      writevalue = writevalue + plusminus_to_string(_u_alt->Cell(bit0))+ "  ";
    bit0.Next();
    }
    yellow("       " +writelocation + " \n ");
    blue("     " + writevalue + "\n");
  
    string writelocation1 = "";
    string writevalue1 = "";
    int counter = 0;
    int counter1 = 0;
    int counter2 = 0;
    
    bit0.SetBoundary(1);
    bit2.SetBoundary(3);
    bit3.SetBoundary(0);
    bit1.SetBoundary(3);
    bit0.First();
    bit2.First();
    bit3.First();
    while (bit0.Valid()){
        counter ++;
        counter1 ++;
        counter2 ++;
        //red(to_string(bit3.Pos()[1]));
        bit0.Next();
    }
    bit0.First();
    
    while (bit2.Valid()){

      bit1.First();

      while (bit1.Valid()){
        if (counter1 == bit1.Pos()[1])
        {
          writevalue1 = plusminus_to_string(_u_alt->Cell(bit1));
          writelocation1 = "(" + to_string(bit1.Pos()[0])  ;
          writelocation1 = writelocation1 +","+ to_string(bit1.Pos()[1]) +")";
          yellow(writelocation1 +" ");
          blue(writevalue1 + "  ");
          counter1 --;    
        }
        bit1.Next();
      }

      init.First();
      while (init.Valid()){
        if (counter2 == init.Pos()[1])
        {

          writevalue1 = plusminus_to_string(_u_alt->Cell(init));
          red(writevalue1 + "  "); 
        }
        init.Next();
      }      
      if (counter2 > 1)
      {
        counter2 --;
      }


      bit0.First();

      while (bit0.Valid()){
        if (counter == bit0.Pos()[1])
        {
          writevalue1 = plusminus_to_string(_u_alt->Cell(bit0));
          writelocation1 = "(" + to_string(bit0.Pos()[0])  ;
          writelocation1 = writelocation1 +","+ to_string(bit0.Pos()[1]) +")";
          blue(writevalue1 +" ");
          yellow(writelocation1 + "\n");
          counter --;    
        }
        bit0.Next();
      }
    
    bit2.Next();
    }
  
    writelocation = "      ";
    writevalue = "     ";
    bit0.SetBoundary(0);
    bit0.First();
    while (bit0.Valid()){
      writelocation = writelocation +"(" + to_string(bit0.Pos()[0]) ;
      writelocation = writelocation +"," + to_string(bit0.Pos()[1]) + ")   |  " ;
      writevalue =  writevalue + plusminus_to_string(_u_alt->Cell(bit0))+ "  ";
    bit0.Next();
    }
    blue(" " + writevalue + "\n");
    yellow(" " +writelocation + " \n ");
  
    } else if (f=='U')
    {
    green("Das ist die Ausgabe für _u \n");
    cout << "\n";  
    string writelocation = "";
    string writevalue = "";
    bit0.SetBoundary(2);
    bit0.First();
    while (bit0.Valid()){
      writelocation = writelocation + "(" + to_string(bit0.Pos()[0])+ "," ; 
      writelocation = writelocation + to_string(bit0.Pos()[1]) + ")   |  ";
      writevalue = writevalue + plusminus_to_string(_u->Cell(bit0))+ "  ";
    bit0.Next();
    }
    yellow("       " +writelocation + " \n ");
    blue("     " + writevalue + "\n");
  
    string writelocation1 = "";
    string writevalue1 = "";
    int counter = 0;
    int counter1 = 0;
    int counter2 = 0;
    
    bit0.SetBoundary(1);
    bit2.SetBoundary(3);
    bit3.SetBoundary(0);
    bit1.SetBoundary(3);
    bit0.First();
    bit2.First();
    bit3.First();
    while (bit0.Valid()){
        counter ++;
        counter1 ++;
        counter2 ++;
        //red(to_string(bit3.Pos()[1]));
        bit0.Next();
    }
    bit0.First();
    
    while (bit2.Valid()){

      bit1.First();

      while (bit1.Valid()){
        if (counter1 == bit1.Pos()[1])
        {
          writevalue1 = plusminus_to_string(_u->Cell(bit1));
          writelocation1 = "(" + to_string(bit1.Pos()[0])  ;
          writelocation1 = writelocation1 +","+ to_string(bit1.Pos()[1]) +")";
          yellow(writelocation1 +" ");
          blue(writevalue1 + "  ");
          counter1 --;    
        }
        bit1.Next();
      }

      init.First();
      while (init.Valid()){
        if (counter2 == init.Pos()[1])
        {

          writevalue1 = plusminus_to_string(_u->Cell(init));
          red(writevalue1 + "  "); 
        }
        init.Next();
      }      
      if (counter2 > 1)
      {
        counter2 --;
      }


      bit0.First();

      while (bit0.Valid()){
        if (counter == bit0.Pos()[1])
        {
          writevalue1 = plusminus_to_string(_u->Cell(bit0));
          writelocation1 = "(" + to_string(bit0.Pos()[0])  ;
          writelocation1 = writelocation1 +","+ to_string(bit0.Pos()[1]) +")";
          blue(writevalue1 +" ");
          yellow(writelocation1 + "\n");
          counter --;    
        }
        bit0.Next();
      }
    
    bit2.Next();
    }
  
    writelocation = "      ";
    writevalue = "     ";
    bit0.SetBoundary(0);
    bit0.First();
    while (bit0.Valid()){
      writelocation = writelocation +"(" + to_string(bit0.Pos()[0]) ;
      writelocation = writelocation +"," + to_string(bit0.Pos()[1]) + ")   |  " ;
      writevalue =  writevalue + plusminus_to_string(_u->Cell(bit0))+ "  ";
    bit0.Next();
    }
    blue(" " + writevalue + "\n");
    yellow(" " +writelocation + " \n ");
    }
    else if (f=='v')
    {
    green("Das ist die Ausgabe für _v_alt \n");
    cout << "\n";  
    string writelocation = "";
    string writevalue = "";
    bit0.SetBoundary(2);
    bit0.First();
    while (bit0.Valid()){
      writelocation = writelocation + "(" + to_string(bit0.Pos()[0])+ "," ; 
      writelocation = writelocation + to_string(bit0.Pos()[1]) + ")   |  ";
      writevalue = writevalue + plusminus_to_string(_v_alt->Cell(bit0))+ "  ";
    bit0.Next();
    }
    yellow("       " +writelocation + " \n ");
    blue("     " + writevalue + "\n");
  
    string writelocation1 = "";
    string writevalue1 = "";
    int counter = 0;
    int counter1 = 0;
    int counter2 = 0;
    
    bit0.SetBoundary(1);
    bit2.SetBoundary(3);
    bit3.SetBoundary(0);
    bit1.SetBoundary(3);
    bit0.First();
    bit2.First();
    bit3.First();
    while (bit0.Valid()){
        counter ++;
        counter1 ++;
        counter2 ++;
        //red(to_string(bit3.Pos()[1]));
        bit0.Next();
    }
    bit0.First();
    
    while (bit2.Valid()){

      bit1.First();

      while (bit1.Valid()){
        if (counter1 == bit1.Pos()[1])
        {
          writevalue1 = plusminus_to_string(_v_alt->Cell(bit1));
          writelocation1 = "(" + to_string(bit1.Pos()[0])  ;
          writelocation1 = writelocation1 +","+ to_string(bit1.Pos()[1]) +")";
          yellow(writelocation1 +" ");
          blue(writevalue1 + "  ");
          counter1 --;    
        }
        bit1.Next();
      }

      init.First();
      while (init.Valid()){
        if (counter2 == init.Pos()[1])
        {

          writevalue1 = plusminus_to_string(_v_alt->Cell(init));
          red(writevalue1 + "  "); 
        }
        init.Next();
      }      
      if (counter2 > 1)
      {
        counter2 --;
      }


      bit0.First();

      while (bit0.Valid()){
        if (counter == bit0.Pos()[1])
        {
          writevalue1 = plusminus_to_string(_v_alt->Cell(bit0));
          writelocation1 = "(" + to_string(bit0.Pos()[0])  ;
          writelocation1 = writelocation1 +","+ to_string(bit0.Pos()[1]) +")";
          blue(writevalue1 +" ");
          yellow(writelocation1 + "\n");
          counter --;    
        }
        bit0.Next();
      }
    
    bit2.Next();
    }
  
    writelocation = "      ";
    writevalue = "     ";
    bit0.SetBoundary(0);
    bit0.First();
    while (bit0.Valid()){
      writelocation = writelocation +"(" + to_string(bit0.Pos()[0]) ;
      writelocation = writelocation +"," + to_string(bit0.Pos()[1]) + ")   |  " ;
      writevalue =  writevalue + plusminus_to_string(_v_alt->Cell(bit0))+ "  ";
    bit0.Next();
    }
    blue(" " + writevalue + "\n");
    yellow(" " +writelocation + " \n ");
    }
    else if (f=='V')
    {
    green("Das ist die Ausgabe für _v \n");
    cout << "\n";  
    string writelocation = "";
    string writevalue = "";
    bit0.SetBoundary(2);
    bit0.First();
    while (bit0.Valid()){
      writelocation = writelocation + "(" + to_string(bit0.Pos()[0])+ "," ; 
      writelocation = writelocation + to_string(bit0.Pos()[1]) + ")   |  ";
      writevalue = writevalue + plusminus_to_string(_v->Cell(bit0))+ "  ";
    bit0.Next();
    }
    yellow("       " +writelocation + " \n ");
    blue("     " + writevalue + "\n");
  
    string writelocation1 = "";
    string writevalue1 = "";
    int counter = 0;
    int counter1 = 0;
    int counter2 = 0;
    
    bit0.SetBoundary(1);
    bit2.SetBoundary(3);
    bit3.SetBoundary(0);
    bit1.SetBoundary(3);
    bit0.First();
    bit2.First();
    bit3.First();
    while (bit0.Valid()){
        counter ++;
        counter1 ++;
        counter2 ++;
        //red(to_string(bit3.Pos()[1]));
        bit0.Next();
    }
    bit0.First();
    
    while (bit2.Valid()){

      bit1.First();

      while (bit1.Valid()){
        if (counter1 == bit1.Pos()[1])
        {
          writevalue1 = plusminus_to_string(_v->Cell(bit1));
          writelocation1 = "(" + to_string(bit1.Pos()[0])  ;
          writelocation1 = writelocation1 +","+ to_string(bit1.Pos()[1]) +")";
          yellow(writelocation1 +" ");
          blue(writevalue1 + "  ");
          counter1 --;    
        }
        bit1.Next();
      }

      init.First();
      while (init.Valid()){
        if (counter2 == init.Pos()[1])
        {

          writevalue1 = plusminus_to_string(_v->Cell(init));
          red(writevalue1 + "  "); 
        }
        init.Next();
      }      
      if (counter2 > 1)
      {
        counter2 --;
      }


      bit0.First();

      while (bit0.Valid()){
        if (counter == bit0.Pos()[1])
        {
          writevalue1 = plusminus_to_string(_v->Cell(bit0));
          writelocation1 = "(" + to_string(bit0.Pos()[0])  ;
          writelocation1 = writelocation1 +","+ to_string(bit0.Pos()[1]) +")";
          blue(writevalue1 +" ");
          yellow(writelocation1 + "\n");
          counter --;    
        }
        bit0.Next();
      }
    
    bit2.Next();
    }
  
    writelocation = "      ";
    writevalue = "     ";
    bit0.SetBoundary(0);
    bit0.First();
    while (bit0.Valid()){
      writelocation = writelocation +"(" + to_string(bit0.Pos()[0]) ;
      writelocation = writelocation +"," + to_string(bit0.Pos()[1]) + ")   |  " ;
      writevalue =  writevalue + plusminus_to_string(_v->Cell(bit0))+ "  ";
    bit0.Next();
    }
    blue(" " + writevalue + "\n");
    yellow(" " +writelocation + " \n ");
    }
    else if (f=='p')
    {
    green("Das ist die Ausgabe für _p \n");
    cout << "\n";  
    string writelocation = "";
    string writevalue = "";
    bit0.SetBoundary(2);
    bit0.First();
    while (bit0.Valid()){
      writelocation = writelocation + "(" + to_string(bit0.Pos()[0])+ "," ; 
      writelocation = writelocation + to_string(bit0.Pos()[1]) + ")   |  ";
      writevalue = writevalue + plusminus_to_string(_p->Cell(bit0))+ "  ";
    bit0.Next();
    }
    yellow("       " +writelocation + " \n ");
    blue("     " + writevalue + "\n");
  
    string writelocation1 = "";
    string writevalue1 = "";
    int counter = 0;
    int counter1 = 0;
    int counter2 = 0;
    
    bit0.SetBoundary(1);
    bit2.SetBoundary(3);
    bit3.SetBoundary(0);
    bit1.SetBoundary(3);
    bit0.First();
    bit2.First();
    bit3.First();
    while (bit0.Valid()){
        counter ++;
        counter1 ++;
        counter2 ++;
        //red(to_string(bit3.Pos()[1]));
        bit0.Next();
    }
    bit0.First();
    
    while (bit2.Valid()){

      bit1.First();

      while (bit1.Valid()){
        if (counter1 == bit1.Pos()[1])
        {
          writevalue1 = plusminus_to_string(_p->Cell(bit1));
          writelocation1 = "(" + to_string(bit1.Pos()[0])  ;
          writelocation1 = writelocation1 +","+ to_string(bit1.Pos()[1]) +")";
          yellow(writelocation1 +" ");
          blue(writevalue1 + "  ");
          counter1 --;    
        }
        bit1.Next();
      }

      init.First();
      while (init.Valid()){
        if (counter2 == init.Pos()[1])
        {

          writevalue1 = plusminus_to_string(_p->Cell(init));
          red(writevalue1 + "  "); 
        }
        init.Next();
      }      
      if (counter2 > 1)
      {
        counter2 --;
      }


      bit0.First();

      while (bit0.Valid()){
        if (counter == bit0.Pos()[1])
        {
          writevalue1 = plusminus_to_string(_p->Cell(bit0));
          writelocation1 = "(" + to_string(bit0.Pos()[0])  ;
          writelocation1 = writelocation1 +","+ to_string(bit0.Pos()[1]) +")";
          blue(writevalue1 +" ");
          yellow(writelocation1 + "\n");
          counter --;    
        }
        bit0.Next();
      }
    
    bit2.Next();
    }
  
    writelocation = "      ";
    writevalue = "     ";
    bit0.SetBoundary(0);
    bit0.First();
    while (bit0.Valid()){
      writelocation = writelocation +"(" + to_string(bit0.Pos()[0]) ;
      writelocation = writelocation +"," + to_string(bit0.Pos()[1]) + ")   |  " ;
      writevalue =  writevalue + plusminus_to_string(_p->Cell(bit0))+ "  ";
    bit0.Next();
    }
    blue(" " + writevalue + "\n");
    yellow(" " +writelocation + " \n ");
    }
    else if (f=='F')
    {
      green("Das ist die Ausgabe für _F \n");
      cout << "\n";  
    string writelocation = "";
    string writevalue = "";
    bit0.SetBoundary(2);
    bit0.First();
    while (bit0.Valid()){
      writelocation = writelocation + "(" + to_string(bit0.Pos()[0])+ "," ; 
      writelocation = writelocation + to_string(bit0.Pos()[1]) + ")   |  ";
      writevalue = writevalue + plusminus_to_string(_F->Cell(bit0))+ "  ";
    bit0.Next();
    }
    yellow("       " +writelocation + " \n ");
    blue("     " + writevalue + "\n");
  
    string writelocation1 = "";
    string writevalue1 = "";
    int counter = 0;
    int counter1 = 0;
    int counter2 = 0;
    
    bit0.SetBoundary(1);
    bit2.SetBoundary(3);
    bit3.SetBoundary(0);
    bit1.SetBoundary(3);
    bit0.First();
    bit2.First();
    bit3.First();
    while (bit0.Valid()){
        counter ++;
        counter1 ++;
        counter2 ++;
        //red(to_string(bit3.Pos()[1]));
        bit0.Next();
    }
    bit0.First();
    
    while (bit2.Valid()){

      bit1.First();

      while (bit1.Valid()){
        if (counter1 == bit1.Pos()[1])
        {
          writevalue1 = plusminus_to_string(_F->Cell(bit1));
          writelocation1 = "(" + to_string(bit1.Pos()[0])  ;
          writelocation1 = writelocation1 +","+ to_string(bit1.Pos()[1]) +")";
          yellow(writelocation1 +" ");
          blue(writevalue1 + "  ");
          counter1 --;    
        }
        bit1.Next();
      }

      init.First();
      while (init.Valid()){
        if (counter2 == init.Pos()[1])
        {

          writevalue1 = plusminus_to_string(_F->Cell(init));
          red(writevalue1 + "  "); 
        }
        init.Next();
      }      
      if (counter2 > 1)
      {
        counter2 --;
      }


      bit0.First();

      while (bit0.Valid()){
        if (counter == bit0.Pos()[1])
        {
          writevalue1 = plusminus_to_string(_F->Cell(bit0));
          writelocation1 = "(" + to_string(bit0.Pos()[0])  ;
          writelocation1 = writelocation1 +","+ to_string(bit0.Pos()[1]) +")";
          blue(writevalue1 +" ");
          yellow(writelocation1 + "\n");
          counter --;    
        }
        bit0.Next();
      }
    
    bit2.Next();
    }
  
    writelocation = "      ";
    writevalue = "     ";
    bit0.SetBoundary(0);
    bit0.First();
    while (bit0.Valid()){
      writelocation = writelocation +"(" + to_string(bit0.Pos()[0]) ;
      writelocation = writelocation +"," + to_string(bit0.Pos()[1]) + ")   |  " ;
      writevalue =  writevalue + plusminus_to_string(_F->Cell(bit0))+ "  ";
    bit0.Next();
    }
    blue(" " + writevalue + "\n");
    yellow(" " +writelocation + " \n ");
    }
    else if (f=='G')
    {
      green("Das ist die Ausgabe für _G \n");
      cout << "\n";  
    string writelocation = "";
    string writevalue = "";
    bit0.SetBoundary(2);
    bit0.First();
    while (bit0.Valid()){
      writelocation = writelocation + "(" + to_string(bit0.Pos()[0])+ "," ; 
      writelocation = writelocation + to_string(bit0.Pos()[1]) + ")   |  ";
      writevalue = writevalue + plusminus_to_string(_G->Cell(bit0))+ "  ";
    bit0.Next();
    }
    yellow("       " +writelocation + " \n ");
    blue("     " + writevalue + "\n");
  
    string writelocation1 = "";
    string writevalue1 = "";
    int counter = 0;
    int counter1 = 0;
    int counter2 = 0;
    
    bit0.SetBoundary(1);
    bit2.SetBoundary(3);
    bit3.SetBoundary(0);
    bit1.SetBoundary(3);
    bit0.First();
    bit2.First();
    bit3.First();
    while (bit0.Valid()){
        counter ++;
        counter1 ++;
        counter2 ++;
        //red(to_string(bit3.Pos()[1]));
        bit0.Next();
    }
    bit0.First();
    
    while (bit2.Valid()){

      bit1.First();

      while (bit1.Valid()){
        if (counter1 == bit1.Pos()[1])
        {
          writevalue1 = plusminus_to_string(_G->Cell(bit1));
          writelocation1 = "(" + to_string(bit1.Pos()[0])  ;
          writelocation1 = writelocation1 +","+ to_string(bit1.Pos()[1]) +")";
          yellow(writelocation1 +" ");
          blue(writevalue1 + "  ");
          counter1 --;    
        }
        bit1.Next();
      }

      init.First();
      while (init.Valid()){
        if (counter2 == init.Pos()[1])
        {

          writevalue1 = plusminus_to_string(_G->Cell(init));
          red(writevalue1 + "  "); 
        }
        init.Next();
      }      
      if (counter2 > 1)
      {
        counter2 --;
      }


      bit0.First();

      while (bit0.Valid()){
        if (counter == bit0.Pos()[1])
        {
          writevalue1 = plusminus_to_string(_G->Cell(bit0));
          writelocation1 = "(" + to_string(bit0.Pos()[0])  ;
          writelocation1 = writelocation1 +","+ to_string(bit0.Pos()[1]) +")";
          blue(writevalue1 +" ");
          yellow(writelocation1 + "\n");
          counter --;    
        }
        bit0.Next();
      }
    
    bit2.Next();
    }
  
    writelocation = "      ";
    writevalue = "     ";
    bit0.SetBoundary(0);
    bit0.First();
    while (bit0.Valid()){
      writelocation = writelocation +"(" + to_string(bit0.Pos()[0]) ;
      writelocation = writelocation +"," + to_string(bit0.Pos()[1]) + ")   |  " ;
      writevalue =  writevalue + plusminus_to_string(_G->Cell(bit0))+ "  ";
    bit0.Next();
    }
    blue(" " + writevalue + "\n");
    yellow(" " +writelocation + " \n ");
    }
    else if (f=='T')
    {
      green("Das ist die Ausgabe für _T \n");
      cout << "\n";  
    string writelocation = "";
    string writevalue = "";
    bit0.SetBoundary(2);
    bit0.First();
    while (bit0.Valid()){
      writelocation = writelocation + "(" + to_string(bit0.Pos()[0])+ "," ; 
      writelocation = writelocation + to_string(bit0.Pos()[1]) + ")   |  ";
      writevalue = writevalue + plusminus_to_string(_T->Cell(bit0))+ "  ";
    bit0.Next();
    }
    yellow("       " +writelocation + " \n ");
    blue("     " + writevalue + "\n");
  
    string writelocation1 = "";
    string writevalue1 = "";
    int counter = 0;
    int counter1 = 0;
    int counter2 = 0;
    
    bit0.SetBoundary(1);
    bit2.SetBoundary(3);
    bit3.SetBoundary(0);
    bit1.SetBoundary(3);
    bit0.First();
    bit2.First();
    bit3.First();
    while (bit0.Valid()){
        counter ++;
        counter1 ++;
        counter2 ++;
        //red(to_string(bit3.Pos()[1]));
        bit0.Next();
    }
    bit0.First();
    
    while (bit2.Valid()){

      bit1.First();

      while (bit1.Valid()){
        if (counter1 == bit1.Pos()[1])
        {
          writevalue1 = plusminus_to_string(_T->Cell(bit1));
          writelocation1 = "(" + to_string(bit1.Pos()[0])  ;
          writelocation1 = writelocation1 +","+ to_string(bit1.Pos()[1]) +")";
          yellow(writelocation1 +" ");
          blue(writevalue1 + "  ");
          counter1 --;    
        }
        bit1.Next();
      }

      init.First();
      while (init.Valid()){
        if (counter2 == init.Pos()[1])
        {

          writevalue1 = plusminus_to_string(_T->Cell(init));
          red(writevalue1 + "  "); 
        }
        init.Next();
      }      
      if (counter2 > 1)
      {
        counter2 --;
      }


      bit0.First();

      while (bit0.Valid()){
        if (counter == bit0.Pos()[1])
        {
          writevalue1 = plusminus_to_string(_T->Cell(bit0));
          writelocation1 = "(" + to_string(bit0.Pos()[0])  ;
          writelocation1 = writelocation1 +","+ to_string(bit0.Pos()[1]) +")";
          blue(writevalue1 +" ");
          yellow(writelocation1 + "\n");
          counter --;    
        }
        bit0.Next();
      }
    
    bit2.Next();
    }
  
    writelocation = "      ";
    writevalue = "     ";
    bit0.SetBoundary(0);
    bit0.First();
    while (bit0.Valid()){
      writelocation = writelocation +"(" + to_string(bit0.Pos()[0]) ;
      writelocation = writelocation +"," + to_string(bit0.Pos()[1]) + ")   |  " ;
      writevalue =  writevalue + plusminus_to_string(_T->Cell(bit0))+ "  ";
    bit0.Next();
    }
    blue(" " + writevalue + "\n");
    yellow(" " +writelocation + " \n ");
    }
  } 
  else{
    red("NOT CORRECT CHAR \n");
  }
}

/* Internal methods for debug prints on console in different colours
*/
void Compute::red(string x) {
  string y = string(ANSI_COLOR_RED);
  string z = string(ANSI_COLOR_RESET);
  cout << y << x << z;
}

void Compute::green(string x) {
  string y = string(ANSI_COLOR_GREEN);
  string z = string(ANSI_COLOR_RESET);
  cout << y << x << z;
}

void Compute::yellow(string x) {
  string y = string(ANSI_COLOR_YELLOW);
  string z = string(ANSI_COLOR_RESET);
  cout << y << x << z;
}

void Compute::blue(string x) {
  string y = string(ANSI_COLOR_BLUE);
  string z = string(ANSI_COLOR_RESET);
  cout << y << x << z;
}

void Compute::magenta(string x) {
  string y = string(ANSI_COLOR_MAGENTA);
  string z = string(ANSI_COLOR_RESET);
  cout << y << x << z;
}

void Compute::cyan(string x) {
  string y = string(ANSI_COLOR_CYAN);
  string z = string(ANSI_COLOR_RESET);
  cout << y << x << z;
}

string Compute::plusminus_to_string(double x) {
  if (x >= 0)
    return "+" + to_string(x);
  else
    return to_string(x);
}
string Compute::to_string_color(int a){
  string value;
  if (a == 0)
  {
    value = string(ANSI_COLOR_BLUE) + to_string(a) + string(ANSI_COLOR_RESET);
  } 
  else if ( a == 10)
  {
    value = to_string(a) ;
  }
  else if ( a == 11)
  {
    value = string(ANSI_COLOR_YELLOW) + to_string(a) + string(ANSI_COLOR_RESET);
  }
  else
  {
    value = string(ANSI_COLOR_GREEN) + to_string(a) + string(ANSI_COLOR_RESET);
  }
  return value;
}


//------------------------------------------------------------------------------
