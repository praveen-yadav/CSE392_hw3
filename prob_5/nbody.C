#include <omp.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <stack>
#include <cstdlib>
#include <cmath>
#include <bitset>

using namespace std;

const int bits_y[13] = {2, 8, 32, 128, 512, 2048, 8192, 32768, 131072,
					  54288, 2097152, 8388608, 33554432};
const int bits_x[13] = {1, 4, 16, 64, 256, 1024, 4096, 16384, 65536,
						262144, 1048576, 4194304, 16777216};

int get_morton_id( const unsigned int xc, const unsigned int yc,
				   const int max_level)
{
	unsigned int xc_mid = 0;
	unsigned int yc_mid = 0;

	for(int i=0; i<max_level; i++){
		xc_mid = xc_mid | (xc<<i) & bits_x[i];
		yc_mid = yc_mid | (yc<<i+1) & bits_y[i];
		
		// xc_mid = xc_mid | (xc<<4) & bits_x[4]; // 0b1000000000;
		// xc_mid = xc_mid | (xc<<3) & bits_x[3]; // 0b0010000000;
		// xc_mid = xc_mid | (xc<<2) & bits_x[2]; // 0b0000100000
		// xc_mid = xc_mid | (xc<<1) & bits_x[1]; // 0b0000001000
		// xc_mid = xc_mid | (xc<<0) & bits_x[0]; // 0b0000000010
		
		// yc_mid = yc_mid | (yc<<5) & bits_y[4]; // 0b1000000000;
		// yc_mid = yc_mid | (yc<<4) & bits_y[3]; // 0b0010000000;
		// yc_mid = yc_mid | (yc<<3) & bits_y[2]; // 0b0000100000
		// yc_mid = yc_mid | (yc<<2) & bits_y[1]; // 0b0000001000
		// yc_mid = yc_mid | (yc<<1) & bits_y[0]; // 0b0000000010
	}
	unsigned int mc_id = xc_mid | yc_mid;
		
	return mc_id;
}

void convert_to_morton_id( const vector<vector <double> >& points,
						   const int np,
						   const int max_level,
						   vector<int>& morton_ids,
						   const double xmin,
						   const double xmax,
						   const double ymin,
						   const double ymax )
{
	const double x_grid_size = (xmax-xmin)/pow(2.,max_level);
	const double y_grid_size = (ymax-ymin)/pow(2.,max_level);

#pragma omp parallel for default(none), shared(points, morton_ids)
	for(int i=0; i<np; i++){

		const int xc = (points[i][0]-xmin)/x_grid_size;
		const int yc = (points[i][1]-ymin)/y_grid_size;

		morton_ids[i]=get_morton_id(xc, yc, max_level);


		// cout<<xc<<" "<<yc<<" "<<morton_ids[i]<<endl;
	}
	
	return;
}

void build_tree( int gids,
				 int** points,
				 int max_pts_per_node,
				 const int max_level)
{


	return;
}

void generate_points(const double xmin,
					 const double xmax,
					 const double ymin,
					 const double ymax,
					 const int np,
					 vector<vector <double> >& points)
{
	srand((unsigned)time(NULL));

	const double xrange = xmax-xmin;
	const double yrange = ymax-ymin;

#pragma omp parallel for default(none), shared(points)
	for(int i=0; i<np; i++){
	
		double x=((double)rand()/(double)RAND_MAX) * xrange + xmin;
		double y=((double)rand()/(double)RAND_MAX) * yrange + ymin;

		points[i][0] = x;
		points[i][1] = y;
	}
	
	return;
}
	

int write_points( const vector<vector <double> >& points,
				  const vector<int>& morton_ids,
				  const int np )
{
	ofstream ofile;
	ofile.open ("points.dat");

	for (int i=0; i<np; i++)
		ofile<<points[i][0]<<" "<<points[i][1]<<" "<<morton_ids[i]<<endl;

	ofile.close();

	return 0;
}

int main()
{
	const double xmin = 0;
	const double xmax = 1;
	const double ymin = 0;
	const double ymax = 1;
	const int np = 10000;
	const int max_level = 10;
	
	vector<vector <double> > points;
	points.resize( np, vector<double>( 2 , 0.0 ) );
	
	generate_points(xmin, xmax, ymin, ymax, np, points);

	vector<int> morton_ids(np, 0);

	const double start=omp_get_wtime();
	convert_to_morton_id( points,
						  np,
						  max_level,
						  morton_ids,
						  xmin,
						  xmax,
						  ymin,
						  ymax );
	const double end=omp_get_wtime();

	cout<<"wall clock time = "<<end-start<<endl;
	
	// write_points(points, morton_ids, np);
	
	return 0;
}

