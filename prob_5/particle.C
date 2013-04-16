#include "particle.h"
#include <cmath>

void particle::gen_coords( const double xmin, const double xrange, 
						   const double ymin, const double yrange,
						   const double mmin, const double mrange )
{
	x = ((double)rand()/(double)RAND_MAX) * xrange + xmin;
	y = ((double)rand()/(double)RAND_MAX) * yrange + ymin;
	m = 1.0;//((double)rand()/(double)RAND_MAX) * mrange + mmin;
	u = 0.0;
	
	return;
}

void particle::gen_coords_cluster( const double xmin, const double xrange, 
								   const double ymin, const double yrange,
								   const double mmin, const double mrange,
								   const double xcenter, const double ycenter,
								   const double rrange
								   )
{
	const double r = ((double)rand()/(double)RAND_MAX) * rrange;
	const double theta = ((double)rand()/(double)RAND_MAX) * 2*pi;

	x = xcenter + r*cos(theta);
	y = ycenter + r*sin(theta);

	// if the point coordinates is out of the domain
	if(x>(xmin+xrange))
		x = x - xrange* int((x-xmin)/xrange);
	if(y>(ymin+yrange))
		y = y - yrange* int((y-ymin)/yrange);
	if(x<xmin)
		x = x + xrange* int((xmin-x)/xrange+1);
	if(y<ymin)
		y = y + yrange* int((ymin-y)/yrange+1);

	m = 1.0;//((double)rand()/(double)RAND_MAX) * mrange + mmin;
	u = 0.0;
	
	return;
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
