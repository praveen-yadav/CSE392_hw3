#include "point.h"

point& point::operator=(const point &rhs) {

	x = rhs.x;
	y = rhs.y;
	
    return *this;  // Return a reference to myself.
}

point& point::operator/=(const double &rhs) {

	x /= rhs;
	y /= rhs;
	
	return *this;
}
