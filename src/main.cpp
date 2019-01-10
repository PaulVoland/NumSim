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

//------------------------------------------------------------------------------
#include "typedef.hpp"
#include "compute.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
#include "precice/SolverInterface.hpp"
#include "%PRECICE_ROOT%"

#ifdef USE_DEBUG_VISU
#include "visu.hpp"
#endif // USE_DEBUG_VISU
#ifdef USE_VTK
#include "vtk.hpp"
#include <sys/stat.h>
#endif // USE_VTK

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#include "zeitgeist.hpp"
#include "argvparser.hpp"

#include <iostream>

using namespace std;

/* Internal methods for debug prints on console in different colours
*/
void red(string x) {
  string y = string(ANSI_COLOR_RED);
  string z = string(ANSI_COLOR_RESET);
  cout << y << x << z;
}

void green(string x) {
  string y = string(ANSI_COLOR_GREEN);
  string z = string(ANSI_COLOR_RESET);
  cout << y << x << z;
}

void yellow(string x) {
  string y = string(ANSI_COLOR_YELLOW);
  string z = string(ANSI_COLOR_RESET);
  cout << y << x << z;
}

void blue(string x) {
  string y = string(ANSI_COLOR_BLUE);
  string z = string(ANSI_COLOR_RESET);
  cout << y << x << z;
}

void magenta(string x) {
  string y = string(ANSI_COLOR_MAGENTA);
  string z = string(ANSI_COLOR_RESET);
  cout << y << x << z;
}

void cyan(string x) {
  string y = string(ANSI_COLOR_CYAN);
  string z = string(ANSI_COLOR_RESET);
  cout << y << x << z;
}

string plusminus_to_string(double x) {
  if (x >= 0)
    return "+" + to_string(x);
  else
    return to_string(x);
}

//-----------------------------------------------------------------------------

// include the namespaces
using namespace precice;
using namespace precice::constants;

//TODO funktionen implementieren

precice.setInterfaceVertices()  // Interface _ coordinaten

precice.writeTemp( Grid, temp, num_coupling_cells, couplingValues, &geom)             //

setCouplingBoundary()           //

//------------------------------------------------------------------------------
int main(int argc, char **argv) {

  /* // Measuring of computational times
  ZeitGeist zg;
  zg.Start(); */

  // Create parameter and geometry instances with default values
  Parameter param;
  // param.Load("../default.param"); // Is done by the parser now.
  Geometry geom;
  // geom.Load("../default.geom"); // Is done by the parser now.

  // Using the ARGVParser.hpp template
  // Works with commands from terminal like: -geom ../default.geom -param ../default.param
  ARGVParser parser;
  parser.bind("-geom", [&geom](int ac, char **av) -> int {
    if (ac != 1) return 0;
    geom.Load(av[0]);
    return 1;
  });
  parser.bind("-param", [&param](int ac, char **av) -> int {
    if (ac != 1) return 0;
    param.Load(av[0]);
    return 1;
  });
  parser.bind("-config") -> int {
    if (ac != 1) return 0;
    string path(av[0]);
    return 1;
  });
  parser.exec(argc, argv);

  // Create the fluid solver
  Compute comp(&geom, &param);


  cout << "Configure preCICE..." << endl;


 // Create preCICE with the solver’s name, the rank, and the total number of processes.
 //string solverName = "Fluid";

  // initialize preCICE
  SolverInterface precice("Fluid", 0, 1);

  precice.configure(path);

  // get domain dimension from precice
  int dim = precice.getDimensions();

  //
  int N = geom.Num_Coupling();

  // ######################################
  multi_real_t temperature(N);
  temperature(geom->Temperature());
 // #########################################

  // get IDs from preCICE
  int meshID        = precice.getMeshID("Fluid-Mesh");
  int tempID        = precice.getDataID("Temperature", meshID);
  int heatFluxID    = precice.getDataID("Heat-Flux", meshID);

  // define coupling mesh and data ids from precice
  //int meshID = precice.getMeshID(mesh_name);

  int* vertexIDs;
  double* grid;
  vertexIDs = new int[N];
  grid = new double[dim * N];

  // array für die Interphasewerte als buffer
  //double* vertices = new double [dim * N]
  //Zuordnung der Koordinaten TODO
  // give mesh to precice
  precice.setMeshVertices(meshID, N, grid, vertexIDs);

  precice.initialize();

   // write initial data if required
  if (precice.isActionRequired(actionWriteInitialData())) {
    precice.writeBlockScalarData(tempID, N, vertexIDs, temperature);
    precice.fulfilledAction(actionWriteInitialData());
  }
  // initial data is sent or received if necessary
  precice.initializeData();
   // read data if available
  if (precice.isReadDataAvailable()) {
    precice.readBlockScalarData(heatFluxID, N, vertexIDs, heatFlux);
  }


 /////// NEW preCICE ////////////////////////////////

 //preCICE parameter
 string precice_config =        //the path to the precice-config.xml file,
 string participant_name        // which should typically be Fluid,
 string mesh_name               // which should typically be Fluid-Mesh,
 string read_data_name          // which should typically be Heat-Flux, and
 string write_data_name

 int* vertexIDs = set_interface_vertices(...);  // get coupling cell ids

// neues Blatt:   int dataID = precice.getDataID("data", meshID); // get your data id from precice,
                                                    // data such as pressure, velocity, temperatur
                                                    // each data should have a seperate id!
double* vertices = new double [dim * num_coupling_cells] // array für die Interphasewerte als buffer

// define Dirichlet part of coupling written by this solver
int temperatureID = precice.getDataID(write_data_name, meshID);
double* temperatureCoupled = double[sizeof(double) * num_coupling_cells];

// define Neumann part of coupling read by this solver
int heatFluxID = precicec.getDataID(read_data_name, meshID);
double* heatfluxCoupled = new double[sizeof(double) * num_coupling_cells];


 //////////////////////////////////////////////////////////////

  // Create Iterator instance for debug purposes
  Iterator it(&geom);
  // Create Boundary Iterator instances for debug purposes
  BoundaryIterator bit0(&geom);
  BoundaryIterator bit1(&geom);
  BoundaryIterator bit2(&geom);
  BoundaryIterator bit3(&geom);
  // Create Interior Iterator instance for debug purposes
  InteriorIterator intit(&geom);

  // To put the picture in a nice order on the screen while execution or to control execution
  bool start = true;

  #ifdef USE_VTK
    // Check if folder "VTK" exists
    struct stat info;

    if (stat("VTK", &info) != 0)
      system("mkdir VTK");
  #endif // USE_VTK

  // Create and initialize the debug visualization
  #ifdef USE_DEBUG_VISU
    Renderer visu(geom.TotalLength(), geom.Mesh());
    double ratio = geom.TotalLength()[1]/geom.TotalLength()[0];
    visu.Init(800, 800*ratio);
  #endif // USE_DEBUG_VISU

  #ifdef USE_VTK
    // Create a VTK generator
    VTK vtk(geom.Mesh(), geom.TotalLength());
  #endif

  #ifdef USE_DEBUG_VISU
    const Grid *visugrid;

    visugrid = comp.GetVelocity();
  #endif // USE_DEBUG_VISU


  // ...
  while (precicec_isCouplingOngoing()) {

    #ifdef USE_DEBUG_VISU
      // Render and check if window is closed
      switch (visu.Render(visugrid)) {
        case -1:
          return -1;
        case 0:
          visugrid = comp.GetVelocity();
          break;
        case 1:
          visugrid = comp.GetU();
          break;
        case 2:
          visugrid = comp.GetV();
          break;
        case 3:
          visugrid = comp.GetP();
          break;
        case 4:
          visugrid = comp.GetT();
        default:
          break;
      };
    #endif // USE_DEBUG_VISU

    // Type in anything to start after positioning the pictures or checking the console prints
    if (start) {
      cin.get();
      #ifndef USE_STEP_BY_STEP
        start = false;
      #endif // NOT USE_STEP_BY_STEP
    }

    #ifdef USE_VTK
      // Create VTK Files in the folder VTK
      vtk.Init("VTK/field");
      vtk.AddField("Velocity", comp.GetU(), comp.GetV());
      vtk.AddScalar("Pressure", comp.GetP());
      vtk.AddScalar("Temperature", comp.GetT());
      vtk.Finish();
    #endif // USE_VTK

///////////////////////////////// NEW preCICE  ///////

    // real_t oder double (precice_dt)
    real_t solver_dt;
    real_t dt;

    //1. calculate time step
      this->Comp_TimeStep(0.0);
      solver_dt = comp.GetTimeStep();

      // call precicec_initialize()
      double precice_dt = precice.initialize();
      //neues Blatt interface.initialize(); //samesame

      dt = min(solver_dt, precice_dt);
     this->Comp_TimeStep(dt);


    //3 - 6. calculate temp, F and G | RHS of P eq. | pressure | new U,V:

    // Run a few steps
    //#ifndef USE_STEP_BY_STEP
    //    for (uint32_t i = 0; i < 199; ++i) {
            comp.TimeStep(false);

//#####################################################
    // temp übergeben schreibe in Array temperatur für jede Zelle aus Array CouplingCellsNumbs(enthält Zell Nummern von allen Coupling Zellen (Temperatur von Nachbar Fluid Zelle )) 
    // entsprechende Werte für die Temperatur
//#####################################################
    //7. coupling
    precice.writeBlockScalarData(...); // write new temperature to preCICE buffers
    precice_dt = precicec_advance(dt); // advance coupling
    precicec.readBlockScalarData(...); // read new heatflux from preCICE buffers 
    //#########################################################################
    // rufe setHeatFlux in geom auf und setze das Array mit HeatFlux neu
    // ########################################################################
    //8. output U, V, P for visualization and update iteration values
    }
    precicec_finalize();


    #endif // NOT USE_STEP_BY_STEP

    // Print coordinates with values for the velocity u TODO
    #ifdef USE_DEBUG_PRINT_U
      /* green("Das ist die Ausgabe für u:\n");
      intit.First();
      string writelocation = "";
      string writevalue = "";
      bit2.SetBoundary(2);
      while (bit2.Valid()) {
        writelocation += "(" + to_string(TestBoundIt.Pos()[0]) + "," + to_string(TestBoundIt.Pos()[1]) + ")   |   ";
        writevalue += plusminus_to_string(comp.GetU()->Cell(bit2)) + "   ";
      TestBoundIt.Next();
      }
      yellow("      " + writelocation + "\n");
      blue("      " + writevalue + "\n");
      */
    #endif // USE_DEBUG_PRINT_U

    // Print coordinates with values for the velocity v TODO
    #ifdef USE_DEBUG_PRINT_V
      /* green("Das ist die Ausgabe für u:\n");
      intit.First();
      string writelocation = "";
      string writevalue = "";
      bit2.SetBoundary(2);
      while (bit2.Valid()) {
        writelocation += "(" + to_string(TestBoundIt.Pos()[0]) + "," + to_string(TestBoundIt.Pos()[1]) + ")   |   ";
        writevalue += plusminus_to_string(comp.GetV()->Cell(bit2)) + "   ";
      TestBoundIt.Next();
      }
      yellow("      " + writelocation + "\n");
      blue("      " + writevalue + "\n");
      */
    #endif // USE_DEBUG_PRINT_V

    // Print coordinates with values for the velocity p TODO
    #ifdef USE_DEBUG_PRINT_P
      /* green("Das ist die Ausgabe für u:\n");
      intit.First();
      string writelocation = "";
      string writevalue = "";
      bit2.SetBoundary(2);
      while (bit2.Valid()) {
        writelocation += "(" + to_string(TestBoundIt.Pos()[0]) + "," + to_string(TestBoundIt.Pos()[1]) + ")   |   ";
        writevalue += plusminus_to_string(comp.GetV()->Cell(bit2)) + "   ";
      TestBoundIt.Next();
      }
      yellow("      " + writelocation + "\n");
      blue("      " + writevalue + "\n");
      */
    #endif // USE_DEBUG_PRINT_P

    // Print field types and neighbourhood from geometry.cpp TODO
    #ifdef USE_DEBUG_PRINT_TYPES
      /* TODO
      */
    #endif // USE_DEBUG_PRINT_TYPES

    comp.TimeStep(true);
  }
  /* if (!comm.getRank()) {
    cout << "Total computational time = "
      << zg.Stop() << " µs\n" << endl;
  }*/
  return 0;
}
