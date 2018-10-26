#include "iterator.hpp"
//------------------------------------------------------------------------------
/// Constructs a new Iterator depending on a geometry
// @param geom input geometry to iterate over
Iterator::Iterator(const Geometry* geom) {
  _geom = geom;
  First();
}

/// Constructs a new Iterator on a geometry with a defined starting value
// @param geom  input geometry to iterate over
// @param value starting index for iteration
Iterator::Iterator(const Geometry* geom, const index_t& value) {
  _geom = geom;
  _value = value;
  _valid = true;
}

/// Returns the current position value
const index_t& Iterator::Value() const {return _value;}

/// Cast operator to convert Iterators to integers
Iterator::operator const index_t&() const {return _value;}

/// Returns the position coordinates
multi_index_t Iterator::Pos() const {
  multi_index_t pos;
  pos[0] = _value % _geom->Size()[0];
  pos[1] = (_value - pos[0]) / _geom->Size()[1];
  return pos;
}

/// Sets the iterator to the first element
void Iterator::First() {
  _value = 0;
  _valid = true;
}

/// Goes to the next element of the iterator, disables it if position is end
void Iterator::Next() {
  multi_index_t size = _geom->Size();
  if (_value + 1 > size[0]*size[1] - 1)
    _valid = false;
  else
    _value++;
}

/// Checks if the iterator still has a valid value
bool Iterator::Valid() const {
  return _valid;
}

/// Returns an Iterator that is located left from this one.
// If we are at the left boundary, the cell sees itself.
Iterator Iterator::Left() const {
  if (_value % _geom->Size()[0] == 0)
    return Iterator(_geom, _value);
  else
    return Iterator(_geom, _value--);
}

/// Returns an Iterator that is located right from this one.
// If we are at the right boundary, the cell sees itself.
Iterator Iterator::Right() const {
  if ((_value + 1) % _geom->Size()[0] == 0)
    return Iterator(_geom, _value);
  else
    return Iterator(_geom, _value++);
}

/// Returns an Iterator that is located above this one.
// If we are at the upper domain boundary, the cell sees itself.
Iterator Iterator::Top() const {
  index_t _value_top = _value + _geom->Size()[1];
  if (_value_top > (_geom->Size()[0]*_geom->Size()[1] - 1))
    return Iterator(_geom, _value);
  else
    return Iterator(_geom, _value_top);
}

/// Returns an Iterator that is located below this one.
// If we are at the lower domain boundary, the cell sees itself.
Iterator Iterator::Down() const {
  int32_t _value_down = _value - _geom->Size()[1];
  if (_value_down < 0)
    return Iterator(_geom, _value);
  else
    return Iterator(_geom, _value_down);
}

//------------------------------------------------------------------------------
/* Iterator for interior cells
*/
/// Constructs a new interior Iterator depending on a geometry
// @param geom input geometry to iterate over
InteriorIterator::InteriorIterator(const Geometry* geom) : Iterator(geom) {
  First();
}

/// Sets the iterator to the first element
void InteriorIterator::First() {
  _value = _geom->Size()[0] + 1;
  _valid = true;
}

/// Goes to the next element of the iterator, disables it if position is end
void InteriorIterator::Next() {
  if ((_value + 1) >= (_geom->Size()[0]*(_geom->Size()[1] - 1) - 1))
    _valid = false;
  else if ((_value + 2) % _geom->Size()[0] == 0)
    _value += 3;
  else
    _value++;
}

//------------------------------------------------------------------------------
/* Iterator for domain boundary cells
*/
/// Constructs a new interior Iterator depending on a geometry
// @param geom input geometry to iterate over
BoundaryIterator::BoundaryIterator(const Geometry* geom) : Iterator(geom) {
  _boundary = 0;
  First();
}

/// Sets the boundary to iterate
//  @param boundary sets which boundary to iterate
//  0 - Bottom boundary
//  1 - Right  boundary
//  2 - Top    boundary
//  3 - Left   boundary
void BoundaryIterator::SetBoundary(const index_t& boundary) {
  _boundary = boundary;
  First();
}

/// Sets the iterator to the first element
void BoundaryIterator::First() {
  _valid = true;

  // Bottom boundary
  if (_boundary == 0)
    _value = 0;
  // Right boundary
  else if (_boundary == 1)
    _value = 2*_geom->Size()[0] - 1;
  // Top boundary
  else if (_boundary == 2)
    _value = _geom->Size()[0]*_geom->Size()[1] - 1;
  // Left boundary
  else if (_boundary == 3)
    _value = _geom->Size()[0]*(_geom->Size()[1] - 2);
}
/// Goes to the next element of the iterator, disables it if position is end
void BoundaryIterator::Next() {
  // Bottom boundary
  if (_boundary == 0) {
    _value++;
    if (_value >= _geom->Size()[0]) {
      _valid = false
    }
  }
  // Right boundary
  if (_boundary == 1) {
    _value += _geom->Size()[1];
    if (_value > _geom->Size()[0]*(_geom->Size()[1] - 1)) {
      _valid = false
    }
  }
  // Top boundary
  if (_boundary == 2) {
    _value--;
    if (_value < _geom->Size()[0]*(_geom->Size()[1] - 1)) {
      _valid = false;
    }
  }
  // Left boundary
  if (_boundary == 3) {
    _value -= _geom->Size()[1];
    if (_value < _geom->Size()[0]) {
      _valid = false;
    }
  }
}
//------------------------------------------------------------------------------