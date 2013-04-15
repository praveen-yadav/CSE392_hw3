// quad-tree class
#ifndef QTREE_H
#define QTREE_H

#include <iostream>
#include <vector>
#include <cstdlib>
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

	// maximum level and allowed points in a box
	int max_level;
	int max_pts;

	// constructor
	qtree(){;}

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
};

void qtree::initialize_root( qtree* ini_parent, int ini_level, point ini_anchor,
							 int ini_max_level, int ini_max_pts,
							 point ini_width, int ini_lid,
							 vector<particle>* points, int num_points,
							 long ini_gid, int st_point,
							 int np_proc )
{
	initialize( ini_parent, ini_level, ini_anchor,
				ini_max_level, ini_max_pts,
				ini_width, ini_lid, points, num_points, ini_gid, st_point );
	idx.resize(np_proc,0);
	for(int i=0; i<np_proc; i++){
		idx[i]=st_point+i;
	}
	lid = 0;
	gid = 0;
};

void qtree::initialize ( qtree* ini_parent, int ini_level, point ini_anchor,
						 int ini_max_level, int ini_max_pts,
						 point ini_width, int ini_lid,
						 vector<particle>* points, int num_points,
						 long ini_gid, int st_point )
{
	lid = ini_lid;
	gid = ini_gid;
	level = ini_level;   // level of the node 
	pts = points;
	// idx
	total_np = num_points;
	// kids
	// parent
	isleaf = 1;
	anchor = ini_anchor;  // [x;y] coordinates of lower left point of a node
	width = ini_width; // width of a box
	parent = ini_parent;  // parent of the node
	max_level = ini_max_level;
	max_pts = ini_max_pts;
	st_pt = st_point;
	
	for(int i=0; i<4; i++)
		kids_exist[i] = 1;
}

int n_th;

int qtree::insert_points( int n_threads )
{
	if (isleaf){

		if((level==max_level) || (idx.size()<max_pts)){
			isleaf=1; // becomes leaf
			return 1;
		}
		
		else create_kids();

	}

	if (n_threads > 1){
		// now insert points to kids
		// for(int i=0; i<4; i++){
		#pragma omp parallel sections
		{
			#pragma omp section
			{
				#pragma omp critical
				cout<<omp_get_thread_num()<<" start"<<endl;
				kids[0].points_in_node(idx);
				if(kids[0].insert_points(n_threads/2))
					kids_exist[0] = 0;
				
				kids[1].points_in_node(idx);
				if(kids[1].insert_points(n_threads/2))
					kids_exist[1] = 0;
				#pragma omp critical
				{
				cout<<omp_get_thread_num()<<" "
					<<kids[0].idx.size()+kids[1].idx.size()<<" end"<<endl;
				}
			}

			#pragma omp section
			{
				#pragma omp critical
				cout<<omp_get_thread_num()<<" start"<<endl;
				kids[2].points_in_node(idx);
				if(kids[2].insert_points(n_threads-n_threads/2))
					kids_exist[2] = 0;
				
				kids[3].points_in_node(idx);
				if(kids[3].insert_points(n_threads-n_threads/2))
					kids_exist[3] = 0;

				#pragma omp critical
				{
				cout<<omp_get_thread_num()<<" "
					<<kids[1].idx.size()+kids[2].idx.size()<<" end"<<endl;
				}

			}

		} 
	}

	// serial
	else{
		for(int i=0; i<4; i++){
			kids[i].points_in_node(idx);
			if(kids[i].insert_points(1))
				kids_exist[i] = 0;
		}
			
	}
	// cout<<"RESTL "<<kids[i].lid<<endl;
	// for(int j=0; j<kids[i].idx.size(); j++)
	// cout<<kids[i].idx[j]<<endl;
		
	
	
	//} end for

	
	return 0;

}

void qtree::create_kids()
{
	int kid_level = level+1;
	point kid_width;
	kid_width.x = width.x/2.0;
	kid_width.y = width.y/2.0;
	point kid_anchor;

	// create new kids objects
	kids  = new qtree[4];
	
	for(int i=0; i<4; i++){
		kid_anchor.x = anchor.x + kid_width.x*kids_pos[i][0];
		kid_anchor.y = anchor.y + kid_width.y*kids_pos[i][1];

		// kid_global_id = 4*parent_global_id + kid_local_id
		kids[i].initialize(this, kid_level, kid_anchor,
						   max_level, max_pts, kid_width, i, pts, total_np,
						   i+gid*4, st_pt );
	}

	// not a leaf any more
	isleaf = 0;
			
	// show_kids();
	
	return;
}


void qtree::show_kids()
{
	// for(int i=0; i<4; i++)
	// 	cout<<"kid "<<i<<": "<<kids[i].anchor.x<<" "<<kids[i].anchor.y<<endl;

}

// check which pointes are in nodes
void qtree::points_in_node(vector<int>& idx_parent)
{
	// number of points that a parent has 
	int npts = idx_parent.size();
		
	for(int i=0; i<npts; i++){
		// int mt_id = (*pts)[real_idx].mt_id;
		int mt_id = (*pts)[idx_parent[i]].mt_id;
		int offset = max_level-level;
	
		if(( (mt_id & mt_checker[offset]) >> (offset*2)  ) == lid){
			idx.push_back(idx_parent[i]);
			// cout<<real_idx<<" "<<(pts)[real_idx].x<<" "<<(pts)[real_idx].y
				// <<" "<<(pts)[real_idx].mt_id<<" is in box "<<lid
				// <<" on level "<<level<<endl;
		}
	}
	
	
}

void qtree::show_tree()
{
	// cout<<"level: "<<level<<endl
	// 	// <<"lid:   "<<lid<<endl
	// 	<<"gid:   "<<gid<<endl
		// <<"n_pts: "<<idx.size()<<endl<<endl;
	for(int i=0; i<4; i++){
		if(!kids[i].isleaf){
			// cout<<"not leaf"<<endl;
			kids[i].show_tree();
		}
		else{
			cout<<"level: "<<kids[i].level
				<<" id: "<<kids[i].gid<<endl;
			for(int j=0; j<kids[i].idx.size(); j++)
				cout<<kids[i].idx[j]<<endl;
		}
		// cout<<endl;
	}

	
}


#endif // QTREE_H
