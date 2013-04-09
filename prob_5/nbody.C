#include <omp.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <stack>
#include <cstdlib>
#include <cmath>
#include <bitset>

#include "../msort/msort.h"
#include "particle.h"

using namespace std;

int compar_vector( const particle left, const particle right)
{
	return left.mt_id<right.mt_id;
}

void build_tree( int gids,
				 int** points,
				 int max_pts_per_node,
				 const int max_level)
{


	return;
}
	

int write_points( const vector<particle>& pts,
				  const int np )
{
	ofstream ofile;
	ofile.open ("points.dat");

	for (int i=0; i<np; i++)
		ofile<<pts[i].x<<" "<<pts[i].y<<" "<<pts[i].mt_id<<endl;

	ofile.close();

	return 0;
}

int main()
{
	const double xmin = 0;
	const double xmax = 1;
	const double ymin = 0;
	const double ymax = 1;
	const int np = 1000;
	const int max_level = 4;

	const double xrange = xmax-xmin;
	const double yrange = ymax-ymin;
	const double x_grid_size = (xmax-xmin)/pow(2.,max_level);
	const double y_grid_size = (ymax-ymin)/pow(2.,max_level);
	
	const int nt = 1;
	
	vector<particle> pts;
	pts.resize(np);

	#pragma omp parallel for default(none), shared(pts)
	for(int i=0; i<np; i++)
		pts[i].gen_coords(xmin, xrange, ymin, yrange);

	#pragma omp parallel for default(none), shared(pts)
	for(int i=0; i<np; i++)
		pts[i].get_morton_id(xmin, ymin, x_grid_size, y_grid_size, max_level);
	
	// tmp space for mergesort
	vector<particle> tmp;
	tmp.resize(np);
	
	// parallel merge sort
	mergesort<particle >(&pts[0], nt, np, &tmp[0], compar_vector);
	
	// cout<<"wall clock time = "<<end-start<<endl;
	
	write_points(pts, np);
	
	return 0;
}

