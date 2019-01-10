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
using namespace precice;
using namespace precice::constants;

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
  parser.bind("-config", -> int {
    if (ac != 1) return 0;
    string path(av[0]);
    return 1;
  });
  parser.exec(argc, argv);

  // Create the fluid solver
  Compute comp(&geom, &param);

  cout << "Configuring preCICE..." << endl;

  // Create preCICE with the solver’s name, the rank, and the total number of processes.
  SolverInterface precice("Fluid", 0, 1);
  // Give path to config file
  precice.configure(path);
  // Get dimension of physical point coordinates from precice
  index_t dim = (index_t)(precice.getDimensions());
  // Get number of coupling cells from geom
  index_t N = geom.Num_Coupling();
  // Fields for exchange
  temp = new double[N];
  heatFlux = new double[N];
  // Get temperature values from fluid field
  for (index_t i = 0; i < N; i++) {
    Iterator it_C(&geom, geom.GetCouplingCellNumbs()[i]);
    if (geom.Cell(it_C).type != typeCoupling)
      break;
    switch (geom.Cell(it_C).fluid) {
      case cellN:
        temp[i] = comp.GetT()->Cell(it_C.Top());
        break;
      case cellW:
        temp[i] = comp.GetT()->Cell(it_C.Left());
        break;
      case cellNW:
        temp[i] = (comp.GetT()->Cell(it_C.Top()) + comp.GetT()->Cell(it_C.Left()))/2.0;
        break;
      case cellS:
        temp[i] = comp.GetT()->Cell(it_C.Down());
        break;
      case cellSW:
        temp[i] = (comp.GetT()->Cell(it_C.Down()) + comp.GetT()->Cell(it_C.Left()))/2.0;
        break;
      case cellE:
        temp[i] = comp.GetT()->Cell(it_C.Right());
        break;
      case cellNE:
        temp[i] = (comp.GetT()->Cell(it_C.Top()) + comp.GetT()->Cell(it_C.Right()))/2.0;
        break;
      case cellSE:
        temp[i] = (comp.GetT()->Cell(it_C.Down()) + comp.GetT()->Cell(it_C.Right()))/2.0;
        break;
      default:
        break;
    };
  }

  // Get IDs from preCICE
  int meshID        = precice.getMeshID("Fluid-Mesh");
  int tempID        = precice.getDataID("Temperature", meshID);
  int heatFluxID    = precice.getDataID("Heat-Flux", meshID);

  int* vertexIDs;
  double* grid;
  vertexIDs = new int[N];
  grid = new double[dim*N];
  //------------------------------------------------------------------------------
  // Manipulate manually!!!
  double* offset = new double[3];
  offset[0] = 0.0;
  offset[1] = 0.25;
  //------------------------------------------------------------------------------ 
  for (index_t i = 0; i < N; i++) {
    for (index_t j = 0; j < dim; j++) {
      Iterator it_C(&geom, geom.GetCouplingCellNumbs()[i]);
      if (j == 2) {
        grid[i*j + j] = 0.0;
      } else {
        grid[i*j + j] = (double)(it_C.Pos()[j]*geom.Mesh()[j] + geom.Mesh()[j]/2.0) + offset[j];
      }
    }  
  }

  // Give mesh to precice
  precice.setMeshVertices(meshID, (int)(N), grid, vertexIDs);

  cout << "Initializing preCICE..." << endl;
  real_t dt;
  real_t precice_dt = (real_t)(precice.initialize());

  // Write initial data if required
  if (precice.isActionRequired(actionWriteInitialData())) {
    precice.writeBlockScalarData(tempID, (int)(N), vertexIDs, temp);
    precice.fulfilledAction(actionWriteInitialData());
  }
  // Initial data is sent or received if necessary
  precice.initializeData();
  // Read data if available
  if (precice.isReadDataAvailable()) {
    precice.readBlockScalarData(heatFluxID, (int)(N), vertexIDs, heatFlux);
  }
  geom.setHeatFlux(heatFlux);

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

  while (precice.isCouplingOngoing()) {

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

    comp.Comp_TimeStep(0.0);
    dt = min(comp.GetTimeStep(), precice_dt);
    comp.Comp_TimeStep(dt);
    comp.TimeStep(true);

    // Get temperature values from fluid field
  for (index_t i = 0; i < N; i++) {
    Iterator it_C(&geom, geom.GetCouplingCellNumbs()[i]);
    if (geom.Cell(it_C).type != typeCoupling)
      break;
    switch (geom.Cell(it_C).fluid) {
      case cellN:
        temp[i] = comp.GetT()->Cell(it_C.Top());
        break;
      case cellW:
        temp[i] = comp.GetT()->Cell(it_C.Left());
        break;
      case cellNW:
        temp[i] = (comp.GetT()->Cell(it_C.Top()) + comp.GetT()->Cell(it_C.Left()))/2.0;
        break;
      case cellS:
        temp[i] = comp.GetT()->Cell(it_C.Down());
        break;
      case cellSW:
        temp[i] = (comp.GetT()->Cell(it_C.Down()) + comp.GetT()->Cell(it_C.Left()))/2.0;
        break;
      case cellE:
        temp[i] = comp.GetT()->Cell(it_C.Right());
        break;
      case cellNE:
        temp[i] = (comp.GetT()->Cell(it_C.Top()) + comp.GetT()->Cell(it_C.Right()))/2.0;
        break;
      case cellSE:
        temp[i] = (comp.GetT()->Cell(it_C.Down()) + comp.GetT()->Cell(it_C.Right()))/2.0;
        break;
      default:
        break;
    };
  }

  precice.writeBlockScalarData(tempID, (int)(N), vertexIDs, temp);
  precice_dt = precice.advance(dt); // advance coupling
  precice.readBlockScalarData(heatFluxID, (int)(N), vertexIDs, heatFlux);
  geom.setHeatFlux(heatFlux);
  
  }

  precice.finalize();


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

  /* if (!comm.getRank()) {
    cout << "Total computational time = "
      << zg.Stop() << " µs\n" << endl;
  }*/
  return 0;
}