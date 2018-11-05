#include "solver.hpp"
#include "geometry.hpp"
#include "grid.hpp"

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

//----------------------------------------------------------------------------
/// Constructor of the abstract Solver class
//@param geom information about the geometry
Solver::Solver(const Geometry* geom) {
    _geom = geom;
    cout << "Constructed a solver for the given geometry." << endl;
}

/// Destructor of the Solver Class (siehe Paul)
Solver::~Solver() {}

/// Returns the residual at [it] for the pressure-Poisson equation
//@param it  Iterator instance
//@param p   current pressure values on the grid
//@param rhs right hand side (grid values)
real_t Solver::localRes(const Iterator& it, const Grid* p, const Grid* rhs) const {
    return fabs(p->dxx(it) + p->dyy(it) - rhs->Cell(it));
}
//------------------------------------------------------------------------------
/* concrete SOR solver
*/
/// Constructs an actual SOR solver
//@param geom  information about the geometry
//@param omega factor for the correction
SOR::SOR(const Geometry* geom, const real_t& omega) {
    _geom = geom;
    _omega = omega;
    cout << "Constructed a SOR solver for the given geometry and omega = " << _omega << "." << endl;
}

/// Constructs an actual SOR solver 'overloaded' (without an omega input)
// => compute optimal omega
// ... improvement around factor 2  => try out for the work sheet
//@param geom information about the geometry
SOR::SOR(const Geometry* geom) {
    const real_t PI = M_PI;
    _geom = geom;
    _omega = 2.0/(1.0 + sin(PI*geom->Mesh()[0]));
    cout << "Constructed the SOR solver for the given geometry with a computed optimal omega = " << _omega << "." << endl;
}

/// Destructor (siehe Paul)
SOR::~SOR() {}

/// Returns the total residual and executes a solver cycle
//@param p   current pressure values on the grid
//@param rhs right hand side (grid values)
real_t SOR::Cycle(Grid* p, const Grid* rhs) const {
    InteriorIterator intit(_geom);
    // Counter
    index_t n = 0;
    // Preparations
    real_t dx = _geom->Mesh()[0];
    real_t dy = _geom->Mesh()[1];
    real_t res_tot = 0.0;
    real_t norm = (dx*dx*dy*dy)/(2*(dx*dx + dy*dy));

    while(intit.Valid()) {
        n++;
        // Define corresponding p values
        real_t p_ij  = p->Cell(intit);
        real_t p_ijl = p->Cell(intit.Left());
        real_t p_ijr = p->Cell(intit.Right());
        real_t p_ijd = p->Cell(intit.Down());
        real_t p_ijt = p->Cell(intit.Top());
        // Corrector
        corr = (p_ijl + p_ijr)/(dx*dx) + (p_ijd + p_ijt)/(dy*dy) - rhs->Cell(intit);
        // New inner p values with SOR solver approach
        p->Cell(intit) = (1 - _omega)*p_ij + _omega*norm*corr;
        // Compute the local residual and sum
        real_t res_loc = localRes(intit, p, rhs);
        res_tot += res_loc;
        // Interior Iterator goes to the next cell
        intit.Next();
    }
    /* NOT HERE!
    // Update boundary values for pressure
    _geom->Update_P(p);
    */
    // Compute norm residual (weighted with number of grid cells)
    return res_tot/n;
}
//------------------------------------------------------------------------------