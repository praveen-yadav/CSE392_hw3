#include "particle.h"
#include <cmath>

void particle::gen_coords( const double xmin, const double xrange, 
						   const double ymin, const double yrange
						   )
{
	x=((double)rand()/(double)RAND_MAX) * xrange + xmin;
	y=((double)rand()/(double)RAND_MAX) * yrange + ymin;

	return;
}

void particle::gen_coords_cluster( const double xmin, const double xrange, 
								   const double ymin, const double yrange,
								   const double xcenter, const double ycenter,
								   const double rrange
								   )
{
	const double r = ((double)rand()/(double)RAND_MAX) * rrange;
	const double theta = ((double)rand()/(double)RAND_MAX) * 2*pi;

	x = xcenter + r*cos(theta);
	y = ycenter + r*sin(theta);

	if(x>(xmin+xrange))
		x = x*(x/xrange-int(x/xrange))+xmin;
	if(y>(ymin+yrange))
		y = y*(y/yrange-int(y/yrange))+ymin;
}

void particle::get_morton_id( const double xmin,
							  const double ymin,
							  const double x_grid_size,
							  const double y_grid_size,
							  const int max_level )
{
	unsigned int xc_mid = 0;
	unsigned int yc_mid = 0;

	const int xc = (x-xmin)/x_grid_size;
	const int yc = (y-ymin)/y_grid_size;

	for(int i=0; i<max_level; i++){
		xc_mid = xc_mid | (xc<<i) & bits_x[i];
		yc_mid = yc_mid | (yc<<i+1) & bits_y[i];
	}
	
	mt_id = xc_mid | yc_mid;
	
	return;
}
