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
#include "communicator.hpp"
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

#include "zeitgeist.hpp"
#include "argvparser.hpp"

#include <iostream>

using namespace std;
int main(int argc, char **argv) {

  /* // Measuring of computational times
  ZeitGeist zg;
  zg.Start(); */

  // Create parameter and geometry instances with default values, set up a communicator
  Communicator comm(&argc, &argv);
  Parameter param;
  // param.Load("../default.param"); // Is done by the parser now.
  Geometry geom(&comm);
  // geom.Load("../default.geom"); // Is done by the parser now.
  // Create the fluid solver
  Compute comp(&geom, &param, &comm);

  // Using the ARGVParser.hpp template
  // Works with commands from terminal like: -geom ../example_1.geom -param ../example_1.param
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

  // To put the single pictures per thread in a nice order on the screen while execution
  bool start = true;

  #ifdef USE_VTK
  if (comm.getRank() == 0) {
    // check if folder "VTK" exists
    struct stat info;

    if (stat("VTK", &info) != 0) {
      system("mkdir VTK");
    }
  }
  #endif

  // Create and initialize the visualization
  #ifdef USE_DEBUG_VISU
    Renderer visu(geom.Length(), geom.Mesh());
    double ratio = geom.Length()[1]/geom.Length()[0];
    visu.Init(800/comm.ThreadDim()[0], 800*ratio/comm.ThreadDim()[1],
      comm.getRank() + 1);
  #endif // USE_DEBUG_VISU

  #ifdef USE_VTK
  // Create a VTK generator;
  // use offset as the domain shift
  multi_real_t offset;
  offset[0] = comm.ThreadIdx()[0]*(geom.Mesh()[0]*(double)(geom.Size()[0] - 2));
  offset[1] = comm.ThreadIdx()[1]*(geom.Mesh()[1]*(double)(geom.Size()[1] - 2));
  VTK vtk(geom.Mesh(), geom.Length(), geom.TotalLength(), offset, comm.getRank(),
    comm.getSize(), comm.ThreadDim());
  #endif

  #ifdef USE_DEBUG_VISU
    const Grid *visugrid;

    visugrid = comp.GetVelocity();
  #endif // USE_DEBUG_VISU

  // Run the time steps until the end is reached
  while (comp.GetTime() < param.Tend()) {
    #ifdef USE_DEBUG_VISU
      // Render and check if window is closed
      switch (visu.Render(visugrid, comm.gatherMin(visugrid->Min()),
        comm.gatherMax(visugrid->Max()))) {
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
        default:
          break;
      };
    #endif // USE_DEBUG_VISU

    // Type in anything to start after positioning the pictures
    if (start) {
      std::cin.get();
      start = false;
    }

    #ifdef USE_VTK
    // Create VTK Files in the folder VTK
    // Note that when using VTK module as it is you first have to write cell
    // information, then call SwitchToPointData(), and then write point data.
    vtk.Init("VTK/field");
    vtk.AddRank();
    vtk.AddCellField("Cell Velocity", comp.GetU(), comp.GetV());
    vtk.SwitchToPointData();
    vtk.AddPointField("Velocity", comp.GetU(), comp.GetV());
    vtk.AddPointScalar("Pressure", comp.GetP());
    vtk.Finish();
    #endif

    // Run a few steps
    for (uint32_t i = 0; i < 9; ++i)
      comp.TimeStep(false);

    bool printOnlyOnMaster = !comm.getRank();
    comp.TimeStep(printOnlyOnMaster);
  }

  /* if (!comm.getRank()) {
    cout << "Total computational time = "
      << zg.Stop() << " µs\n" << endl;
  }*/
  
  return 0;
}
