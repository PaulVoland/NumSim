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

  // Run the time steps until the end is reached
  while (comp.GetTime() < param.Tend()) {
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
    if(true){
    green("Das ist die Ausgabe für u \n");  
    string writelocation = "";
    string writevalue = "";
    bit0.SetBoundary(2);
    bit0.First();
    while (bit0.Valid()){
      writelocation = writelocation + "(" + to_string(bit0.Pos()[0])+ "," ; 
      writelocation = writelocation + to_string(bit0.Pos()[1]) + ")   |  ";
      writevalue = writevalue + plusminus_to_string(comp.GetU()->Cell(bit0))+ "  ";
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
          writevalue1 = plusminus_to_string(comp.GetU()->Cell(bit1));
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

          writevalue1 = plusminus_to_string(comp.GetU()->Cell(init));
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
          writevalue1 = plusminus_to_string(comp.GetU()->Cell(bit0));
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
      writevalue =  writevalue + plusminus_to_string(comp.GetU()->Cell(bit0))+ "  ";
    bit0.Next();
    }
    blue(" " + writevalue + "\n");
    yellow(" " +writelocation + " \n ");
  }
    #endif // USE_DEBUG_PRINT_U

    // Print coordinates with values for the velocity v TODO
    #ifdef USE_DEBUG_PRINT_V
    if(true){
    green("Das ist die Ausgabe für v \n");  
    string writelocation = "";
    string writevalue = "";
    bit0.SetBoundary(2);
    bit0.First();
    while (bit0.Valid()){
      writelocation = writelocation + "(" + to_string(bit0.Pos()[0])+ "," ; 
      writelocation = writelocation + to_string(bit0.Pos()[1]) + ")   |  ";
      writevalue = writevalue + plusminus_to_string(comp.GetV()->Cell(bit0))+ "  ";
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
          writevalue1 = plusminus_to_string(comp.GetV()->Cell(bit1));
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

          writevalue1 = plusminus_to_string(comp.GetV()->Cell(init));
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
          writevalue1 = plusminus_to_string(comp.GetV()->Cell(bit0));
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
      writevalue =  writevalue + plusminus_to_string(comp.GetV()->Cell(bit0))+ "  ";
    bit0.Next();
    }
    blue(" " + writevalue + "\n");
    yellow(" " +writelocation + " \n ");
  }
    #endif // USE_DEBUG_PRINT_V

    // Print coordinates with values for the velocity p TODO
    #ifdef USE_DEBUG_PRINT_P
   if(true){
    green("Das ist die Ausgabe für p \n");  
    string writelocation = "";
    string writevalue = "";
    bit0.SetBoundary(2);
    bit0.First();
    while (bit0.Valid()){
      writelocation = writelocation + "(" + to_string(bit0.Pos()[0])+ "," ; 
      writelocation = writelocation + to_string(bit0.Pos()[1]) + ")   |  ";
      writevalue = writevalue + plusminus_to_string(comp.GetP()->Cell(bit0))+ "  ";
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
          writevalue1 = plusminus_to_string(comp.GetP()->Cell(bit1));
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

          writevalue1 = plusminus_to_string(comp.GetP()->Cell(init));
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
          writevalue1 = plusminus_to_string(comp.GetP()->Cell(bit0));
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
      writevalue =  writevalue + plusminus_to_string(comp.GetP()->Cell(bit0))+ "  ";
    bit0.Next();
    }
    blue(" " + writevalue + "\n");
    yellow(" " +writelocation + " \n ");
  }
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