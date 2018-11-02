#include "solver.hpp"
#include <cmath>
#include <iostream>

//----------------------------------------------------------------------------

/// Constructor of the abstract Solver class
//@param geom information about the geometry
Solver::Solver(const Geometry* geom){
    _geom = geom;
    cout << "Constructed a Solver for the given geometry." << endl;
}


/// Destructor of the Solver Class (siehe Paul)
  //??virtual
Solver::~Solver(){}


/// Returns the residual at [it] for the pressure-Poisson equation
//@param it Iterator instance
//@param grid current pressure values
//@param rhs right hand side
real_t Solver::localRes(const Iterator& it, const Grid* grid, const Grid* rhs) const{
    //welche Seite ist sicher größer? Betrag?
    real_t _res = grid->dxx(it) + grid->dyy(it) - rhs->Cell(it);
    return{_res};
}

//KANN KOMPLET WEG? Oder was ist die Bedeutung?

/// This function must be implemented in a child class
// @param [in][out] grid current values
// @param [in]      rhs  right hand side values
// @returns accumulated residual
  // ??VIRTUAL und warum const = 0 ?
//real_t Cycle(Grid* grid, const Grid* rhs) const = 0 { }


//------------------------------------------------------------------------------

/** concrete SOR solver
*/

/// Constructs an actual SOR solver
//@param geom information about the geometry
//@param omega factor for the correction
SOR::SOR(const Geometry* geom, const real_t& omega){
    _geom = geom;
    _omega = omega;
    cout << "Constructed the SOR_solver for the given geometry and omega =" << _omega << endl;
}

/// Constructs an actual SOR solver, 'overloaded' (without an omega input)
//-> compute optimal omega
// ... Verbesserung um Faktor 2  -> ausprobieren für Fragen
//@param geom information about the geometry
SOR::SOR(const Geometry* geom){

    const real_t PI = 3.141592653589793;
    _geom = geom;
    _omega = 2.0/(1.0 + sin(PI*geom->Mesh()[0]));
    cout << "Constructed the SOR_solver for the given geometry with a computed optimal omega =" << _omega << endl;
}


/// Destructor (siehe Paul)
SOR::~SOR(){}

/// Returns the total residual and executes a solver cycle
// @param grid current pressure values
// @param rhs right hand side
real_t SOR::Cycle(Grid* grid, const Grid* rhs) const{

    InteriorIterator it = InteriorIterator(_geom);

    real_t dx = _geom->Mesh()[0];
    real_t dy = _geom->Mesh()[1];
    real_t correction = 0.0;
    real_t res_loc= 0.0;
    real_t res_tot = 0.0;
    real_t res_sum = 0.0;

    while(it.Valide()){

        correction = (grid->dxx(it) + grid->dyy(it) - rhs->Cell(it)) /(2.0/(dx*dx)+2.0/(dy*dy));

        //set p_ij^(n+1) = p_ij^(n) + omega*correction (Vgl VL)
        grid->Cell(it) = grid->Cell(it) + _omega * correction;

        //compute the local residual and sum the squares
        res_loc = localRes(it, grid, rhs);
        res_sum = res_sum + res_loc * res_loc;

        //iterator goes to the next cell
        it.Next();
    }

    //compute total residual with the L^2-norm
    res_tot = sqrt(res_sum/(_geom->Size()[0]*_geom->Size()[1]));
    return {res_tot};

    //Update boundary values for pressure
    _geom->Update_P(grid);
}


//------------------------------------------------------------------------------
