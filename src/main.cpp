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
#define USE_VTK
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


void blue(string x)
{
  string y = string(ANSI_COLOR_BLUE);
  string z = string(ANSI_COLOR_RESET);
  cout << y << x << z;
}
void red(string x)
{
  string y = string(ANSI_COLOR_RED);
  string z = string(ANSI_COLOR_RESET);
  cout << y << x << z;
}
void green(string x)
{
  string y = string(ANSI_COLOR_GREEN);
  string z = string(ANSI_COLOR_RESET);
  cout << y << x << z;
}
void yellow(string x)
{
  string y = string(ANSI_COLOR_YELLOW);
  string z = string(ANSI_COLOR_RESET);
  cout << y << x << z;
}
void magenta(string x)
{
  string y = string(ANSI_COLOR_MAGENTA);
  string z = string(ANSI_COLOR_RESET);
  cout << y << x << z;
}
void cyan(string x)
{
  string y = string(ANSI_COLOR_CYAN);
  string z = string(ANSI_COLOR_RESET);
  cout << y << x << z;
}

string plusminus_to_string(double x)
{
  if (x >= 0)
  {
    return "+"+to_string(x);
  } else {
    return to_string(x);
  }

}

int main(int argc, char **argv) {

  /* // Measuring of computational times
  ZeitGeist zg;
  zg.Start(); */
  string writelocation1 = "";
  int county = 0;
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

  Iterator TestFullIt(&geom);
  // Create Boundary Iterator
  BoundaryIterator TestBoundIt(&geom);
  BoundaryIterator TestBoundIt2(&geom);
  BoundaryIterator TestBoundIt3(&geom);
  BoundaryIterator TestBoundIt22(&geom);
  // Create Interior Iterator
  InteriorIterator TestInterIt(&geom);

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
        case 4:
          visugrid = comp.GetT();
        default:
          break;
      };
    #endif // USE_DEBUG_VISU

    // Type in anything to start after positioning the pictures
    if (start) {
      std::cin.get();
      //start = false;
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
    vtk.AddPointScalar("Temperature", comp.GetT());
    vtk.Finish();
    #endif

    // Run a few steps
    //for (uint32_t i = 0; i < 9; ++i)
     // comp.TimeStep(false);

//print Coordinates with values for the value u
#ifdef USE_DEBUG_PRINT_U 
  if(true){
    green("Das ist die Ausgabe für u \n");
    TestInterIt.First();
    string writelocation = "";
    string writevalue = "";

    TestBoundIt.SetBoundary(2);
    TestBoundIt.First();
    while (TestBoundIt.Valid()){
      writelocation = writelocation + "(" + to_string(TestBoundIt.Pos()[0])+ "," ;
      writelocation = writelocation + to_string(TestBoundIt.Pos()[1]) + ")   |  ";
      writevalue = writevalue +  plusminus_to_string(comp.GetU()->Cell(TestBoundIt)) + "  ";
    TestBoundIt.Next();
    }
    yellow("        " +writelocation + " \n ");
    blue("     " + writevalue + "\n");
  
    string writelocation1 = "";
    string writevalue1 = "";
    int counter = 0;
    int counter1 = 0;
    int counterfirst = 0;
    
    
    TestBoundIt.SetBoundary(1);
    TestBoundIt2.SetBoundary(3);
    TestBoundIt3.SetBoundary(0);
    TestBoundIt.First();
    TestBoundIt2.First();
    TestBoundIt3.First();
    TestBoundIt22.SetBoundary(3);
    TestBoundIt22.First();
    TestFullIt.First();
    while (TestBoundIt.Valid()){
        counter ++;
        counter1 ++;
        counterfirst ++;
        //red(to_string(TestBoundIt3.Pos()[1]));
        TestBoundIt.Next();
    }

    TestBoundIt.First();
    while (TestBoundIt2.Valid()){

      TestBoundIt22.First();
      //blue(to_string(counterfirst));
      while (TestBoundIt22.Valid()){
        //blue(to_string(counterfirst)); 
        //green(to_string(TestBoundIt22.Pos()[1]));
        if (counterfirst == TestBoundIt22.Pos()[1])
        {
          //blue(to_string(counterfirst));
          writelocation1 = "(" + to_string(TestBoundIt22.Pos()[0])  ;
          writelocation1 = writelocation1 +","+ to_string(TestBoundIt22.Pos()[1]) +")";
          writevalue1 = plusminus_to_string(comp.GetU()->Cell(TestBoundIt22));
          yellow(writelocation1 + " ");
          blue(writevalue1 +" ");
          counterfirst --;    
        }
        TestBoundIt22.Next();
      }

      TestBoundIt.First();
      TestInterIt.First();
      TestFullIt.First();
      TestBoundIt3.First();
        cout << " ";
      while(TestBoundIt3.Valid()){  
            TestBoundIt3.Next();
            TestFullIt.Next();
          }
          TestBoundIt3.First();
      while (TestBoundIt.Valid()){
        if (counter == TestBoundIt.Pos()[1])
        {
          TestBoundIt3.Next();
          TestBoundIt3.Next();
          TestFullIt.Next();
          while(TestBoundIt3.Valid()){
          red(plusminus_to_string(comp.GetU()->Cell(TestFullIt))+"  ");
          TestBoundIt3.Next();
          TestFullIt.Next();
          }
        TestFullIt.Next();
        counter1 --;
        TestBoundIt3.First();
        }
        else{
          while(TestBoundIt3.Valid()){  
            TestBoundIt3.Next();
            TestFullIt.Next();
          }
        counter1 --;
        TestBoundIt3.First();
        }
        TestBoundIt.Next();
      }
        
      TestBoundIt.First();

      while (TestBoundIt.Valid()){
        if (counter == TestBoundIt.Pos()[1])
        {
          writevalue1 = plusminus_to_string(comp.GetU()->Cell(TestBoundIt));
          writelocation1 = "(" + to_string(TestBoundIt.Pos()[0])  ;
          writelocation1 = writelocation1 +","+ to_string(TestBoundIt.Pos()[1]) +")";
          blue(writevalue1 +" ");
          yellow(writelocation1 + "\n");
          counter --;    
        }
        TestBoundIt.Next();
      }
    
    TestBoundIt2.Next();
    }
  
    writelocation = "      ";
    writevalue = "     ";
    TestBoundIt.SetBoundary(0);
    TestBoundIt.First();
    while (TestBoundIt.Valid()){
      writelocation = writelocation +"(" + to_string(TestBoundIt.Pos()[0]) ;
      writelocation = writelocation +"," + to_string(TestBoundIt.Pos()[1]) + ")   |  " ;
      writevalue =  writevalue + plusminus_to_string(comp.GetU()->Cell(TestBoundIt))+ "  ";
    TestBoundIt.Next();
    }
    blue(writevalue + "\n");
    yellow("  " +writelocation + " \n ");
  }
  TestFullIt.First();
#endif // #ifdef USE_DEBUG_PRINT_U
  
#ifdef USE_DEBUG_PRINT_V  
//print Coordinates with values for the value v 
  if(true){
    green("Das ist die Ausgabe für v \n");
    TestInterIt.First();
    string writelocation = "";
    string writevalue = "";

    TestBoundIt.SetBoundary(2);
    TestBoundIt.First();
    while (TestBoundIt.Valid()){
      writelocation = writelocation + "(" + to_string(TestBoundIt.Pos()[0])+ "," ;
      writelocation = writelocation + to_string(TestBoundIt.Pos()[1]) + ")   |  ";
      writevalue = writevalue +  plusminus_to_string(comp.GetV()->Cell(TestBoundIt)) + "  ";
    TestBoundIt.Next();
    }
    yellow("        " +writelocation + " \n ");
    blue("     " + writevalue + "\n");
  
    string writelocation1 = "";
    string writevalue1 = "";
    int counter = 0;
    int counter1 = 0;
    int counterfirst = 0;
    
    
    TestBoundIt.SetBoundary(1);
    TestBoundIt2.SetBoundary(3);
    TestBoundIt3.SetBoundary(0);
    TestBoundIt.First();
    TestBoundIt2.First();
    TestBoundIt3.First();
    TestBoundIt22.SetBoundary(3);
    TestBoundIt22.First();
    TestFullIt.First();
    while (TestBoundIt.Valid()){
        counter ++;
        counter1 ++;
        counterfirst ++;
        //red(to_string(TestBoundIt3.Pos()[1]));
        TestBoundIt.Next();
    }

    TestBoundIt.First();
    while (TestBoundIt2.Valid()){

      TestBoundIt22.First();
      //blue(to_string(counterfirst));
      while (TestBoundIt22.Valid()){
        //blue(to_string(counterfirst)); 
        //green(to_string(TestBoundIt22.Pos()[1]));
        if (counterfirst == TestBoundIt22.Pos()[1])
        {
          //blue(to_string(counterfirst));
          writelocation1 = "(" + to_string(TestBoundIt22.Pos()[0])  ;
          writelocation1 = writelocation1 +","+ to_string(TestBoundIt22.Pos()[1]) +")";
          writevalue1 = plusminus_to_string(comp.GetV()->Cell(TestBoundIt22));
          yellow(writelocation1 + " ");
          blue(writevalue1 +" ");
          counterfirst --;    
        }
        TestBoundIt22.Next();
      }

      TestBoundIt.First();
      TestInterIt.First();
      TestFullIt.First();
      TestBoundIt3.First();
        cout << " ";
      while(TestBoundIt3.Valid()){  
            TestBoundIt3.Next();
            TestFullIt.Next();
          }
          TestBoundIt3.First();
      while (TestBoundIt.Valid()){
        if (counter == TestBoundIt.Pos()[1])
        {
          TestBoundIt3.Next();
          TestBoundIt3.Next();
          TestFullIt.Next();
          while(TestBoundIt3.Valid()){
          red(plusminus_to_string(comp.GetV()->Cell(TestFullIt))+"  ");
          TestBoundIt3.Next();
          TestFullIt.Next();
          }
        TestFullIt.Next();
        counter1 --;
        TestBoundIt3.First();
        }
        else{
          while(TestBoundIt3.Valid()){  
            TestBoundIt3.Next();
            TestFullIt.Next();
          }
        counter1 --;
        TestBoundIt3.First();
        }
        TestBoundIt.Next();
      }
        
      TestBoundIt.First();

      while (TestBoundIt.Valid()){
        if (counter == TestBoundIt.Pos()[1])
        {
          writevalue1 = plusminus_to_string(comp.GetV()->Cell(TestBoundIt));
          writelocation1 = "(" + to_string(TestBoundIt.Pos()[0])  ;
          writelocation1 = writelocation1 +","+ to_string(TestBoundIt.Pos()[1]) +")";
          blue(writevalue1 +" ");
          yellow(writelocation1 + "\n");
          counter --;    
        }
        TestBoundIt.Next();
      }
    
    TestBoundIt2.Next();
    }
  
    writelocation = "      ";
    writevalue = "     ";
    TestBoundIt.SetBoundary(0);
    TestBoundIt.First();
    while (TestBoundIt.Valid()){
      writelocation = writelocation +"(" + to_string(TestBoundIt.Pos()[0]) ;
      writelocation = writelocation +"," + to_string(TestBoundIt.Pos()[1]) + ")   |  " ;
      writevalue =  writevalue + plusminus_to_string(comp.GetV()->Cell(TestBoundIt))+ "  ";
    TestBoundIt.Next();
    }
    blue(writevalue + "\n");
    yellow("  " +writelocation + " \n ");
  }
  TestFullIt.First();
#endif //#ifdef USE_DEBUG_PRINT_V
  
#ifdef USE_DEBUG_PRINT_P
//print Coordinates with values for the value p 
  if(true){
    green("Das ist die Ausgabe für P \n");
    TestInterIt.First();
    string writelocation = "";
    string writevalue = "";
    
    TestBoundIt.SetBoundary(2);
    TestBoundIt.First();
    while (TestBoundIt.Valid()){
      writelocation = writelocation + "(" + to_string(TestBoundIt.Pos()[0])+ "," ;
      writelocation = writelocation + to_string(TestBoundIt.Pos()[1]) + ")   |  ";
      writevalue = writevalue +  plusminus_to_string(comp.GetP()->Cell(TestBoundIt)) + "  ";
    TestBoundIt.Next();
    }
    yellow("        " +writelocation + " \n ");
    blue("     " + writevalue + "\n");
  
    string writelocation1 = "";
    string writevalue1 = "";
    int counter = 0;
    int counter1 = 0;
    int counterfirst = 0;

    
    TestBoundIt.SetBoundary(1);
    TestBoundIt2.SetBoundary(3);
    TestBoundIt3.SetBoundary(0);
    TestBoundIt.First();
    TestBoundIt2.First();
    TestBoundIt3.First();
    TestBoundIt22.SetBoundary(3);
    TestBoundIt22.First();
    TestFullIt.First();
    while (TestBoundIt.Valid()){
        counter ++;
        counter1 ++;
        counterfirst ++;
        //red(to_string(TestBoundIt3.Pos()[1]));
        TestBoundIt.Next();
    }

    TestBoundIt.First();
    while (TestBoundIt2.Valid()){

      TestBoundIt22.First();
      //blue(to_string(counterfirst));
      while (TestBoundIt22.Valid()){
        //blue(to_string(counterfirst)); 
        //green(to_string(TestBoundIt22.Pos()[1]));
        if (counterfirst == TestBoundIt22.Pos()[1])
        {
          //blue(to_string(counterfirst));
          writelocation1 = "(" + to_string(TestBoundIt22.Pos()[0])  ;
          writelocation1 = writelocation1 +","+ to_string(TestBoundIt22.Pos()[1]) +")";
          writevalue1 = plusminus_to_string(comp.GetP()->Cell(TestBoundIt22));
          yellow(writelocation1 + " ");
          blue(writevalue1 +" ");
          counterfirst --;    
        }
        TestBoundIt22.Next();
      }

      TestBoundIt.First();
      TestInterIt.First();
      TestFullIt.First();
      TestBoundIt3.First();
        cout << " ";
      while(TestBoundIt3.Valid()){  
            TestBoundIt3.Next();
            TestFullIt.Next();
          }
          TestBoundIt3.First();
      while (TestBoundIt.Valid()){
        if (counter == TestBoundIt.Pos()[1])
        {
          TestBoundIt3.Next();
          TestBoundIt3.Next();
          TestFullIt.Next();
          while(TestBoundIt3.Valid()){
          red(plusminus_to_string(comp.GetP()->Cell(TestFullIt))+"  ");
          TestBoundIt3.Next();
          TestFullIt.Next();
          }
        TestFullIt.Next();
        counter1 --;
        TestBoundIt3.First();
        }
        else{
          while(TestBoundIt3.Valid()){  
            TestBoundIt3.Next();
            TestFullIt.Next();
          }
        counter1 --;
        TestBoundIt3.First();
        }
        TestBoundIt.Next();
      }
        
      TestBoundIt.First();

      while (TestBoundIt.Valid()){
        if (counter == TestBoundIt.Pos()[1])
        {
          writevalue1 = plusminus_to_string(comp.GetP()->Cell(TestBoundIt));
          writelocation1 = "(" + to_string(TestBoundIt.Pos()[0])  ;
          writelocation1 = writelocation1 +","+ to_string(TestBoundIt.Pos()[1]) +")";
          blue(writevalue1 +" ");
          yellow(writelocation1 + "\n");
          counter --;    
        }
        TestBoundIt.Next();
      }
    
    TestBoundIt2.Next();
    }
  
    writelocation = "      ";
    writevalue = "     ";
    TestBoundIt.SetBoundary(0);
    TestBoundIt.First();
    while (TestBoundIt.Valid()){
      writelocation = writelocation +"(" + to_string(TestBoundIt.Pos()[0]) ;
      writelocation = writelocation +"," + to_string(TestBoundIt.Pos()[1]) + ")   |  " ;
      writevalue =  writevalue + plusminus_to_string(comp.GetP()->Cell(TestBoundIt))+ "  ";
    TestBoundIt.Next();
    }
    blue(writevalue + "\n");
    yellow("  " +writelocation + " \n ");
  }

#endif //#ifdef USE_DEBUG_PRINT_P

#ifdef USE_DEBUG_PRINT_ENUM_BORDER
  TestFullIt.First();
  
  
  while (TestFullIt.Valid())
  {
    writelocation1 = "(" + to_string(TestFullIt.Pos()[0])  ;
    writelocation1 = writelocation1 +","+ to_string(TestFullIt.Pos()[1]) +")";
    blue(writelocation1 +" ");
    red(to_string(geom.Cell(TestFullIt).fluid));
    blue(" ");
    if (county == 9)
    {
      blue("\n ");
      county =-1;
    }
    county ++;

    TestFullIt.Next();
  }
  #endif //#ifdef USE_DEBUG_PRINT_ENUM_BORDER

    

    bool printOnlyOnMaster = !comm.getRank();
    comp.TimeStep(printOnlyOnMaster);
  





  }



  /* if (!comm.getRank()) {
    cout << "Total computational time = "
      << zg.Stop() << " µs\n" << endl;
  }*/
  
  return 0;
}
