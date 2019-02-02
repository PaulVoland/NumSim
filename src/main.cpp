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
  parser.exec(argc, argv);

  // Create the fluid solver
  Compute comp(&geom, &param);
  

  // Create Iterator instance for debug purposes
  Iterator it(&geom);
  // Create Boundary Iterator instances for debug purposes
  BoundaryIterator bit0(&geom);
  BoundaryIterator bit1(&geom);
  BoundaryIterator bit2(&geom);
  BoundaryIterator bit3(&geom);
  // Create Interior Iterator instance for debug purposes
  InteriorIterator intit(&geom);
  InteriorIterator init(&geom);
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
    if (geom.TotalLength()[0] > geom.TotalLength()[1]) {
      double ratio = geom.TotalLength()[1]/geom.TotalLength()[0];
      visu.Init(600, 600*ratio);
    } else {
      double ratio = geom.TotalLength()[0]/geom.TotalLength()[1];
      visu.Init(600*ratio, 600);
    }
  #endif // USE_DEBUG_VISU

  #ifdef USE_VTK
    // Create a VTK generator
    VTK vtk(geom.Mesh(), geom.TotalLength());
  #endif

  #ifdef USE_DEBUG_VISU
    const Grid *visugrid;

    visugrid = comp.GetVelocity();
  #endif // USE_DEBUG_VISU

  // Run the time steps until the end is reached
  while (comp.GetTime() < param.Tend()) {
    #ifdef USE_DEBUG_VISU
      // Render and check if window is closed
      switch (visu.Render(visugrid, visugrid->Min(), visugrid->Max())) {
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
      Grid *temp1 = new Grid(&geom);
      Grid *temp2 = new Grid(&geom);
      index_t count = geom.TotalSize()[0]*geom.TotalSize()[1];
      for (index_t i = 0; i < count; ++i)
      {
        Iterator it(&geom,i);
        temp1->Cell(it) = geom.Cell(it).type;
        temp2->Cell(it) = geom.Cell(it).neighbour;
      }

      // Create VTK Files in the folder VTK
      vtk.Init("VTK/field");
      vtk.AddField("Velocity", comp.GetU(), comp.GetV());
      vtk.AddField("Type and Neighbourhood", temp1 , temp2);
      vtk.AddScalar("Pressure", comp.GetP());
      vtk.AddScalar("Temperature", comp.GetT());
      vtk.Finish();
    #endif // USE_VTK

    // Run a few steps
    #ifndef USE_STEP_BY_STEP
      for (uint32_t i = 0; i < 9; ++i) {
        comp.TimeStep(false);
      }
    #endif // NOT USE_STEP_BY_STEP

    // Print coordinates with values for the velocity u TODO
    #ifdef USE_DEBUG_PRINT_U
    comp.ShowVelocitysPressures('U');
    #endif // USE_DEBUG_PRINT_U

    #ifdef USE_DEBUG_PRINT_U_OLD
    comp.ShowVelocitysPressures('u');
    #endif // USE_DEBUG_PRINT_U

    // Print coordinates with values for the velocity v TODO
    #ifdef USE_DEBUG_PRINT_V
    comp.ShowVelocitysPressures('V');
    #endif // USE_DEBUG_PRINT_V

    #ifdef USE_DEBUG_PRINT_V_OLD
    comp.ShowVelocitysPressures('v');
    #endif // USE_DEBUG_PRINT_U

    // Print coordinates with values for the velocity p TODO
    #ifdef USE_DEBUG_PRINT_P
    comp.ShowVelocitysPressures('p');
    #endif // USE_DEBUG_PRINT_P

    #ifdef USE_DEBUG_PRINT_T
    comp.ShowVelocitysPressures('T');
    #endif // USE_DEBUG_PRINT_T

    #ifdef USE_DEBUG_PRINT_F
    comp.ShowVelocitysPressures('F');
    #endif // USE_DEBUG_PRINT_F

    #ifdef USE_DEBUG_PRINT_G
    comp.ShowVelocitysPressures('G');
    #endif // USE_DEBUG_PRINT_F

    // Print field types and neighbourhood from geometry.cpp TODO
    #ifdef USE_DEBUG_PRINT_TYPES
    comp.ShowType();
    #endif // USE_DEBUG_PRINT_TYPES
    
    #ifdef USE_DEBUG_PRINT_PARTICLE
    comp.ShowParticle();
    #endif // USE_DEBUG_PRINT_PARTICLE

    #ifdef USE_DEBUG_PRINT_NEIGHBOUR
    comp.ShowNeighbour();
    #endif // USE_DEBUG_PRINT_NEIGHBOUR



    comp.TimeStep(true);
  }
  /* if (!comm.getRank()) {
    cout << "Total computational time = "
      << zg.Stop() << " µs\n" << endl;
  }*/
  return 0;
}