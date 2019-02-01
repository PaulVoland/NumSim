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

#include "typedef.hpp"
#include <string>
//------------------------------------------------------------------------------
#ifndef __COMPUTE_HPP
#define __COMPUTE_HPP
//------------------------------------------------------------------------------
class Compute {
public:
  /// Creates a compute instance with given geometry and parameter
  Compute(Geometry *geom, const Parameter *param);
  /// Deletes all grids
  ~Compute();

  /// Execute one time step of the fluid simulation (with or without debug info)
  // @ param printInfo print information about current solver state (residual
  // etc.)
  void TimeStep(bool printInfo);

  /// Returns the simulated time in total
  const real_t &GetTime() const;

  /// Returns the pointer to U
  const Grid *GetU() const;
  /// Returns the pointer to V
  const Grid *GetV() const;
  /// Returns the pointer to P
  const Grid *GetP() const;
  /// Returns the pointer to RHS
  const Grid *GetRHS() const;
  /// Returns the pointer to T
  const Grid *GetT() const;

  /// Computes and returns the absolute velocity
  const Grid *GetVelocity();
  /// Computes and returns the vorticity
  const Grid *GetVorticity();
  /// Computes and returns the stream line values
  const Grid *GetStream();
  // Show pressure velocity and temperatur 
  // Options 'v' 'u' for old velocities, 'V' 'U' new Velocitys  
  // 'F' 'G' 'p' 'T' for the other dynamic parameter 
  void ShowVelocitysPressures(char f);
  // Show Particle Trace
  void ShowParticle();
  // Show Neighbour
  void ShowNeighbour();
  // Show Type
  void ShowType();

private:
  // current timestep
  real_t _t;

  // current step width
  real_t _dt;

  // donor-cell diffusion condition (p. 27)
  real_t _dtlimit;

  // limit for residual
  real_t _epslimit;

  // particle trace array
  vec_arr _part_trace = vec_arr();

  // counter for particles per cell
  index_t* _ppc;

  // velocities
  Grid *_u;
  Grid *_v;
  Grid *_u_alt;
  Grid *_v_alt;

  // pressure
  Grid *_p;

  // temperature
  Grid *_T;

  // prel. vel
  Grid *_F;
  Grid *_G;

  // right-hand side
  Grid *_rhs;

  // container for interpolating whichever values
  Grid *_tmp;

  Solver *_solver;

  Geometry *_geom;
  const Parameter *_param;

  /// Compute the new velocites u,v
  void NewVelocities(const real_t &dt);
  /// Compute the temporary velocites F,G
  void MomentumEqu(const real_t &dt);
  /// Compute the RHS of the poisson equation
  void RHS(const real_t &dt);
  /// Compute the new temperature values
  void HeatTransport(const real_t &dt);

  /// Set particles initially to the _part_trace vector
  void SetParticles();
  /// Random number between max and min
  real_t RandNumb(const real_t &max, const real_t &min) const;
  /// Set new particles in each timestep at the inflow boundaries
  void SetNewInflowParticles();
  /// Particle trace per timestep
  void ParticleTrace(const real_t &dt);
  /// Return the index position of a multi_index-type
  index_t IndexToCell(const multi_index_t &value) const;
  /// Calculate the index positions from physical coordinates
  multi_index_t PhysToIndex(const real_t &x , const real_t &y) const;
  /// Get special velocity defined by char f in {u, U, v, V} at physical position
  real_t PhysToVelocity(const real_t &x , const real_t &y, const char &f) const;
  /// Copy velocities of the old timestep to _._alt
  void CopyVelocities();
  //  Particle to Cell Number + Output
  void ShowParticleToCellDebug(const real_t &x , const real_t &y);
  
  // Colors 
  void red(std::string x);
  void green(std::string x);
  void yellow(std::string x);
  void blue(std::string x);
  void magenta(std::string x);
  void cyan(std::string x);
  // formating string 
  std::string plusminus_to_string(double x);
};
//------------------------------------------------------------------------------
#endif // __COMPUTE_HPP
