// 2d point class
#ifndef POINT_H
#define POINT_H

class point
{
public:
	double x;
	double y;

	point(double r, double s){x=r;y=s;}
	point(){;}
	
	point & operator=(const point &rhs);
	point & operator/=(const double &rhs);
};

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


#endif
