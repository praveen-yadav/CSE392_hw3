// particles class
#ifndef PARTICLES_H
#define PARTICLES_H

#include <iostream>
#include <vector>
#include <cstdlib>

using namespace std;

const int bits_y[13] = {2, 8, 32, 128, 512, 2048, 8192, 32768, 131072,
							524288, 2097152, 8388608, 33554432};
const int bits_x[13] = {1, 4, 16, 64, 256, 1024, 4096, 16384, 65536,
							262144, 1048576, 4194304, 16777216};


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
	long mt_id;
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
