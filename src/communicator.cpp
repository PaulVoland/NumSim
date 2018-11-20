#include "communicator.hpp"
#include "geometry.hpp"
#include "iterator.hpp"
#include <iostream>
#include <cmath>
#include <unistd.h>
#include <mpi.h>
#include <stdio.h>

using namespace std;
//------------------------------------------------------------------------------
/** Communicator constructor; initializes MPI Environment
  *
  * \param [in] argc Number of arguments program was started with
  * \param [in] argv Arguments passed to the program on start
  */
Communicator::Communicator(int* argc, char*** argv) {
  // Inizialize size and rank of each run
  MPI_Init(argc,argv);
  MPI_Comm_size(MPI_COMM_WORLD, &_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &_rank);

  /* Partition first in x-direction (vertical slices if odd)
  // Alternative version
  index_t max_div = (index_t)sqrt(_size);
  for (index_t i = max_div; i > 0; i--) {
    if (_size % i == 0) {
      max_div = i;
      break;
    }
  }
  // Evaluation of the numbers of processes in each direction
  _tdim[0] = _size/max_div;
  _tdim[1] = max_div;
  // Set thread indices
  _tidx[1] = (int)(_rank/_tdim[0]);
  _tidx[0] = (int)(_rank - _tidx[1]*_tdim[0]);
  */

  // Partition first in y-direction (horizontal slices if odd)
  for (index_t i= (index_t)sqrt(_size); i > 0; i--) {
    if ((real_t)(_size)/(real_t)(i) == _size/i) {
      _tdim = multi_index_t(i, _size/i);
      break;
    }
  }
  /* Location of main field like this (including numbering, see further down)
  * _tidx[1]
  *   |-------|-------|
  * 2 |   4  x|y  5   |
  *   |-------|-------|
  * 1 |   2  x|y  3   |  _tdim[] = (2, 3)
  *   |-------|-------|
  * 0 |   0  x|y  1   |
  *   |-------|-------|
  *       0       1    _tidx[0]
  */
  // Fix indices for each rank
  _tidx[0] = _rank % _tdim[0];
  _tidx[1] = (index_t)(_rank/_tdim[0]);
  /*
  // Print information on the command line
  printf("Rank %d has number %d in x-direction.\n", _rank, _tidx[0]);
  printf("Rank %d has number %d in y-direction.\n", _rank, _tidx[1]);
  */
  // True/False value which describes if communicator has an even/odd number
  _evenodd = ((_tidx[0]%2==0) && (_tidx[1]%2==0) ||
    (_tidx[0]%2 != 0) && (_tidx[1]%2 != 0) )? true:false;
}
//------------------------------------------------------------------------------
/** Communicator destructor; finalizes MPI Environment
 */
Communicator::~Communicator(){
  MPI_Finalize();
}
//------------------------------------------------------------------------------
/** Returns the position of the current process with respect to the
 *  fields lower left corner
 */
const multi_index_t& Communicator::ThreadIdx() const {return _tidx;}
//------------------------------------------------------------------------------
/** Returns the way the domain is partitioned among all processes
 */
const multi_index_t& Communicator::ThreadDim() const {return _tdim;}
//------------------------------------------------------------------------------
/** Returns whether this process is a red or a black field
 */
const bool& Communicator::EvenOdd() const {return _evenodd;}
//------------------------------------------------------------------------------
/** Gets the sum of all values and distributes the result among all
 *  processes
 *
 * \param [in] val The data over which the sum is to be calculated
 */
real_t Communicator::gatherSum(const real_t& val) const {
 real_t sum = 0;
 MPI_Allreduce(&val,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
 return sum;
}
//------------------------------------------------------------------------------
/** Finds the minimum of the values and distributes the result among
 *  all processes
 *
 * \param [in] val The data over which to find the minimum
 */
real_t Communicator::gatherMin(const real_t& val) const {
 real_t min = 0;
 MPI_Allreduce(&val,&min,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
 return min;
}
//------------------------------------------------------------------------------
/** Finds the maximum of the values and distributes the result among
 *  all processes
 *
 * \param [in] val The data over which to find the maximum
 */
real_t Communicator::gatherMax(const real_t& val) const {
 real_t max = 0;
 MPI_Allreduce(&val,&max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
 return max;
}
//------------------------------------------------------------------------------
/** Synchronizes ghost layer
 *
 * \param [in] grid  The values to sync
 */
void Communicator::copyBoundary(Grid* grid) const {
  // Take care about the numbering of the communicators (even/odd)
  // even grids communicate to the left first, while odd grids to the right first
  if (_evenodd) {
    (!isLeft()   && copyLeftBoundary(grid)); //possibly !isLeft() not necessary, because already in Update_.
    (!isRight()  && copyRightBoundary(grid));
    (!isTop()    && copyTopBoundary(grid));
    (!isBottom() && copyBottomBoundary(grid));
  }
  else {
    (!isRight()  && copyRightBoundary(grid));
    (!isLeft()   && copyLeftBoundary(grid));
    (!isBottom() && copyBottomBoundary(grid));
    (!isTop()    && copyTopBoundary(grid));
  }
}
//------------------------------------------------------------------------------
/** Decide whether our left boundary is a domain boundary
 */
const bool Communicator::isLeft() const {return (_tidx[0] == 0)? true:false;}
//------------------------------------------------------------------------------
/** Decide whether our right boundary is a domain boundary
 */
const bool Communicator::isRight() const {return (_tidx[0] == _tdim[0] - 1)? true:false;}
//------------------------------------------------------------------------------
/** Decide whether our top boundary is a domain boundary
 */
const bool Communicator::isTop() const {return (_tidx[1] == _tdim[1] - 1)? true:false;}
//------------------------------------------------------------------------------
/** Decide whether our bottom boundary is a domain boundary
 */
const bool Communicator::isBottom() const {return (_tidx[1] == 0)? true:false;}
//------------------------------------------------------------------------------
/** Get MPI rank of current process
 */
const int& Communicator::getRank() const {return _rank;}
//------------------------------------------------------------------------------
/** Get number of MPI processes
 */
const int& Communicator::getSize() const {return _size;}
//------------------------------------------------------------------------------
/** Function to sync ghost layer on left boundary:
 *  send values of own left boundary to left neighbor and
 *  and receive values from his right boundary
 *
 *   ------------ ------------
 *  |           x|y           |
 *  |           x|y           |
 *  |           x|y           |
 *  |           x|y           |
 *  |           x|y           |
 *   ------------ ------------
 *
 *   y: values that are sent
 *   x: values that are received
 *
 * \param [in] grid  values whose boundary shall be synced
 */
bool Communicator::copyLeftBoundary(Grid* grid) const {
  MPI_Status stat;
  const index_t messagelength =
    grid->getGeometry()->Size()[1] - 2; // Left boundary has less elements
  real_t message[messagelength];
  BoundaryIterator BoundIt(grid->getGeometry());
  BoundIt.SetBoundary(3);
  int count = 0;
  while (BoundIt.Valid()) {
    message[count] = grid->Cell(BoundIt);
    BoundIt.Next();
    count++;
  }
  int result = MPI_Sendrecv_replace(message, messagelength, MPI_DOUBLE,
    _rank - 1, 1, _rank - 1, 2, MPI_COMM_WORLD, &stat); //_tidx[0] - 1?!

  count = 0;
  BoundIt.First();
  while (BoundIt.Valid()) {
    grid->Cell(BoundIt) = message[count];
    BoundIt.Next();
    count++;
  }
  return (result == MPI_SUCCESS)? true:false;
}
//------------------------------------------------------------------------------
/** Function to sync ghost layer on right boundary
 *  Details analog to left boundary
 *
 * \param [in] grid  values whose boundary shall be synced
 */
bool Communicator::copyRightBoundary(Grid* grid) const {
  MPI_Status stat;
  const index_t messagelength =
    grid->getGeometry()->Size()[1] - 2; // Right boundary has less elements
  real_t message[messagelength];
  BoundaryIterator BoundIt(grid->getGeometry());
  BoundIt.SetBoundary(1);
  int count = 0;
  while (BoundIt.Valid()) {
    message[count] = grid->Cell(BoundIt);
    BoundIt.Next();
    count ++;
  }
  int result = MPI_Sendrecv_replace(message, messagelength, MPI_DOUBLE,
    _rank + 1, 2, _rank + 1, 1, MPI_COMM_WORLD, &stat);

  count = 0;
  BoundIt.First();
  while (BoundIt.Valid()) {
    grid->Cell(BoundIt) = message[count];
    BoundIt.Next();
    count++;
  }
  return (result == MPI_SUCCESS)? true:false;
}
//------------------------------------------------------------------------------
/** Function to sync ghost layer on top boundary
 *  Details analog to left boundary
 *
 * \param [in] grid  values whose boundary shall be synced
 */
bool Communicator::copyTopBoundary(Grid* grid) const {
  MPI_Status stat;
  int toprank = _rank + _tdim[0];
  const index_t messagelength = grid->getGeometry()->Size()[0];
  real_t message[messagelength];
  BoundaryIterator BoundIt(grid->getGeometry());
  BoundIt.SetBoundary(2);
  int count = 0;
  while (BoundIt.Valid()) {
    message[count] = grid->Cell(BoundIt);
    BoundIt.Next();
    count++;
  }
  int result = MPI_Sendrecv_replace(message, messagelength, MPI_DOUBLE,
    toprank, 4, toprank, 3, MPI_COMM_WORLD, &stat);

  count = 0;
  BoundIt.First();
  while (BoundIt.Valid()) {
    grid->Cell(BoundIt) = message[count];
    BoundIt.Next();
    count++;
  }
  return (result == MPI_SUCCESS)? true:false;
}
//------------------------------------------------------------------------------
/** Function to sync ghost layer on bottom boundary
 *  Details analog to left boundary
 *
 * \param [in] grid  values whose boundary shall be synced
 */
bool Communicator::copyBottomBoundary(Grid* grid) const {
  MPI_Status stat;
  int bottomrank = _rank - _tdim[0];
  const index_t messagelength = grid->getGeometry()->Size()[0];
  real_t message[messagelength];
  BoundaryIterator BoundIt(grid->getGeometry());
  BoundIt.SetBoundary(0);
  int count = 0;
  while (BoundIt.Valid()) {
    message[count] = grid->Cell(BoundIt);
    BoundIt.Next();
    count++;
  }
  int result = MPI_Sendrecv_replace(message, messagelength, MPI_DOUBLE,
    bottomrank, 3, bottomrank, 4, MPI_COMM_WORLD, &stat);

  count = 0;
  BoundIt.First();
  while (BoundIt.Valid()) {
    grid->Cell(BoundIt) = message[count];
    BoundIt.Next();
    count++;
  }
  return (result == MPI_SUCCESS)? true:false;
}
//------------------------------------------------------------------------------
