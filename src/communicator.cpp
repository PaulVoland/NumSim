#include "compute.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
#include "grid.hpp"
#include "solver.hpp"
#include "iterator.hpp"

#include <iostream>
#include <cmath>

using namespace std;

//------------------------------------------------------------------------------
class Communicator {
public:
  /** Communicator constructor; initializes MPI Environment
   *
   * \param [in] argc Number of arguments program was started with
   * \param [in] argv Arguments passed to the program on start
   */
  Communicator::Communicator(int *argc, char ***argv){
      MPI_Init(NULL, NULL);
      MPI_Comm_size(MPI_COMM_WORLD, _size);
      MPI_Comm_rank(MPI_COMM_WORLD, _rank);
  }
  /** Communicator destructor; finalizes MPI Environment
   */
  ~Communicator(){
    delete[] _rank;
    delete[] _size;
    delete[] _tidx;
    delete[] _evenodd;
    delete[] _tdim; 
  }

  /** Returns the position of the current process with respect to the
   *  fields lower left corner
   */
  const multi_index_t &ThreadIdx() const{ return _tidx; }

  /** Returns the way the domain is partitioned among all processes
   */
  const multi_index_t &ThreadDim() const{ return _tdim; }

  /** Returns whether this process is a red or a black field
   */
  const bool &EvenOdd() const { return _evenodd; }

  /** Gets the sum of all values and distributes the result among all
   *  processes
   *
   * \param [in] val The data over which the sum is to be calculated
   */
  real_t gatherSum(const real_t &val) const {
    if (_rank == 0){
      
    }
  }
  return max;
  }

  /** Finds the minimum of the values and distributes the result among
   *  all processes
   *
   * \param [in] val The data over which to find the minimum
   */
  real_t gatherMin(const real_t &val) const;

  /** Finds the maximum of the values and distributes the result among
   *  all processes
   *
   * \param [in] val The data over which to find the maximum
   */
  real_t gatherMax(const real_t &val) const;

  /** Synchronizes ghost layer
   *
   * \param [in] grid  The values to sync
   */
  void copyBoundary(Grid *grid) const;

  /** Decide whether our left boundary is a domain boundary
   */
  const bool isLeft() const;

  /** Decide whether our right boundary is a domain boundary
   */
  const bool isRight() const;

  /** Decide whether our top boundary is a domain boundary
   */
  const bool isTop() const;

  /** Decide whether our bottom boundary is a domain boundary
   */
  const bool isBottom() const;

  /** Get MPI rank of current process
   */
  const int &getRank() const;

  /** Get number of MPI processes
   */
  const int &getSize() const;

private:
  multi_index_t _tidx;
  multi_index_t _tdim;
  int _rank;
  int _size;
  bool _evenodd;

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
  bool copyLeftBoundary(Grid *grid) const;

  /** Function to sync ghost layer on right boundary
   *  Details analog to left boundary
   *
   * \param [in] grid  values whose boundary shall be synced
   */
  bool copyRightBoundary(Grid *grid) const;

  /** Function to sync ghost layer on top boundary
   *  Details analog to left boundary
   *
   * \param [in] grid  values whose boundary shall be synced
   */
  bool copyTopBoundary(Grid *grid) const;

  /** Function to sync ghost layer on bottom boundary
   *  Details analog to left boundary
   *
   * \param [in] grid  values whose boundary shall be synced
   */
  bool copyBottomBoundary(Grid *grid) const;
};
//------------------------------------------------------------------------------
