// particles class
#ifndef PARTICLES_H
#define PARTICLES_H

#include <iostream>
#include <vector>
#include <cstdlib>
#include "bits_numbers.h"

using namespace std;

class particle
{
public:
	void gen_coords( const double xmin, const double xrange,
					 const double ymin, const double yrange,
					 const double mmin, const double mrange );

	void gen_coords_cluster( const double xmin, const double xrange, 
							 const double ymin, const double yrange,
							 const double mmin, const double mrange,
							 const double xcenter, const double ycenter,
							 const double rrange
							 );
	// coordinates
	double x;
	double y;
	// morton id
	#ifndef LONG_LONG
	unsigned int mt_id;
	#endif
	#ifdef LONG_LONG
	unsigned long long mt_id;
	#endif
	
	// mass
	double m;
	// potential
	double u;
	
	// particle(int);
	void get_morton_id( const double xmin,
						const double ymin,
						const double x_grid_size,
						const double y_grid_size,
						const int max_level );
};

const double pi=3.14159265359;

#endif // PARTICLES_H
