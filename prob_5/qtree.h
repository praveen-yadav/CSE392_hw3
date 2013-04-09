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
	int gid; // global node id
  	int level;   // level of the node 
	
	// particles data
	particle* pts;

	// particle global ids that a node has
	vector<int> idx;

	int np; // number of points

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

	// constructor for the root
	qtree(qtree* pr, int lv, point anch, int mxlv, int mxpts,
		  point wid, int id, particle* pcl, int n  ){
		initialize( pr, lv, anch, mxlv, mxpts,
					wid, id, pcl, n, 0 );
		idx.resize(n,0);
		for(int i=0; i<n; i++)
			idx[i]=i;
		lid = 0;
		gid = 0;
	}
	// constructor for the children
	qtree(){;}

	void initialize(qtree* pr, int lv, point anch, int mxlv, int mxpts,
					point wid, int l_id, particle* pcl, int n, int g_id );
	int insert_points( );
	void create_kids();
	void show_kids();
	void points_in_node( vector<int>& idx_parent );
	void show_tree();
};

void qtree::initialize( qtree* pr, int lv, point anch, int mxlv, int mxpts,
						point wid, int l_id, particle* pcl, int n,
						int g_id)
{
	lid = l_id;
	gid = g_id;
	level = lv;   // level of the node 
	pts = pcl;
	// idx
	np = n;
	// kids
	// parent
	isleaf = 1;
	anchor = anch;  // [x;y] coordinates of lower left point of a node
	width = wid; // width of a box
	parent = pr;  // parent of the node
	max_level = mxlv;
	max_pts = mxpts;
	for(int i=0; i<4; i++)
		kids_exist[i] = 1;
}


int qtree::insert_points(  )
{

	if (isleaf){

		if((level==max_level) || (idx.size()<max_pts)){
			isleaf=1; // becomes leaf
			return 1;
		}
		
		else create_kids();

	}

	return 0;

}

void qtree::create_kids()
{
	int kid_level = level+1;
	point kid_width;
	kid_width.x = width.x/2.0;
	kid_width.y = width.y/2.0;

	point kid_anchor;

	kids  = new qtree[4];
	
	for(int i=0; i<4; i++){
		kid_anchor.x = anchor.x + kid_width.x*kids_pos[i][0];
		kid_anchor.y = anchor.y + kid_width.y*kids_pos[i][1];

		kids[i].initialize(this, kid_level, kid_anchor,
						   max_level, max_pts, kid_width, i, pts, np,
						   i+gid*4 );
	}

	// not leaf any more
	isleaf = 0;
	
	// now insert points to kids
	for(int i=0; i<4; i++){
		kids[i].points_in_node(idx);

		if(kids[i].insert_points())
			kids_exist[i] = 0;
	}
		
	// show_kids();
	
	return;
}


void qtree::show_kids()
{
	// for(int i=0; i<4; i++)
	// 	cout<<"kid "<<i<<": "<<kids[i].anchor.x<<" "<<kids[i].anchor.y<<endl;

}

void qtree::points_in_node(vector<int>& idx_parent)
{
	int npts = idx_parent.size();

	for(int i=0; i<npts; i++){
		int mt_id = (pts)[idx_parent[i]].mt_id;
		int offset = max_level-level;
		if(( (mt_id & mt_checker[offset]) >> (offset*2)  ) == lid){
			idx.push_back(idx_parent[i]);
			// cout<<(*pts)[i].x<<" "<<(*pts)[i].y
			// 	<<" "<<(*pts)[i].mt_id<<" is in box "<<lid
			// 	<<" on level "<<level<<endl;
		}
	}


}

void qtree::show_tree()
{
	cout<<"level: "<<level<<endl
		// <<"lid:   "<<lid<<endl
		<<"gid:   "<<gid<<endl
		<<"n_pts: "<<idx.size()<<endl<<endl;
	for(int i=0; i<4; i++){
		if(! kids[i].isleaf)
			kids[i].show_tree();
	}

	
}


#endif // QTREE_H
