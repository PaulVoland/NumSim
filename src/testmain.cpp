#include "typedef.hpp"
#include "communicator.hpp"
//#include "grid.hpp"
#include <iostream>
#include <sys/stat.h>
#include <unistd.h>
#include <mpi.h>
#include <stdio.h>

int main(int argc, char **argv) {

  // Create parameter and geometry instances with default values
  Communicator comm(&argc, &argv);
  int rank = 0;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf(" das ist rank %d  \n",rank );
  real_t x = (real_t)(rank);
  real_t sum = 0.0;
  sum = comm.gatherMax(x);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  rank = -1;
  rank = comm.getRank();

  for (int a =0 ; a<11; a++ ){
  	sleep(1);
  	if (a == rank){
  		bool hello = comm.isLeft();
  		printf(" das ist rank %d und der wert lautet %d  \n",rank, hello );
  	}
  }

  






}