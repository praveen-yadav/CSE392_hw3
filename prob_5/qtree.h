// quad-tree class
#ifndef QTREE_H
#define QTREE_H

#include <iostream>
#include <vector>
#include <cstdlib>
#include <omp.h>
#include "point.h"
#include "particle.h"

using namespace std;

// to check morton ids
const int mt_checker[13] =
	{ 3, 12, 48, 192, 768, 3072, 12288, 49152, 196608, 786432, 3145728,
	  12582912, 50331648};

const int kids_pos[4][2] =
    {
        {0,0},
        {1,0},
		{0,1},
		{1,1}
    };

class qtree
{
public:
	// using nid and level, we can identify the node
	int lid; // local node id [0, 1, 2, 3]
	long gid; // global node id
  	int level;   // level of the node 
	
	// particles data
	vector<particle>* pts;

	// particle global ids that a node has
	vector<int> idx;

	// starting point in the array for parallel
	int st_pt;
	
	int total_np; // number of points

	// pointers to relatives
	qtree* kids; // 4-size array
	qtree* parent;  
	int isleaf;  // true or false
	int kids_exist[4]; // 0 for no 1 for yes
	
	// dimensional data
	point anchor;  // [x;y] coordinates of lower left point of a node
	point width; // width of a box

	// centroid data
	double total_mass;
	point centroid; 
	
	// maximum level and allowed points in a box
	int max_level;
	int max_pts;

	// constructor
	qtree(){;}

	// destructor
	~qtree();
	
	void initialize ( qtree* ini_parent, int ini_level, point ini_anchor,
					  int ini_max_level, int ini_max_pts,
					  point ini_width, int ini_lid,
					  vector<particle>* points, int num_points,
					  long ini_gid, int st_point );
	void initialize_root( qtree* ini_parent, int ini_level, point ini_anchor,
						  int ini_max_level, int ini_max_pts,
						  point ini_width, int ini_lid,
						  vector<particle>* points, int num_points,
						  long ini_gid, int st_point,
						  int np_local);
	
	int insert_points( int n_threads=2 );
	void create_kids();
	void show_kids();
	void points_in_node( vector<int>& idx_parent );
	void show_tree();
	void get_centroid( );

};


#endif // QTREE_H
