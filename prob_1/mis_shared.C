#include <omp.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <cstdlib>
// #include <cmath>

using namespace std;

// Algorithm 36.1
// 1. Create a set S of candidates for I as follows. For each vertex v in
// parallel, include v \in S with probability 1/2d(v)
// 2. For each edge in E , if both its endpoints are in S , discard the one of
// lower degree; ties are resolved arbitrary (say by vertex number).
// The resulting set is I.

// comparison function for sorting pairs
bool comp_pairs( const pair<int, int>& i, const pair<int, int>& j ) {
    if( i.first < j.first ) return true;
	else if( i.first == j.first) return i.second<j.second;
	else return false;
}

// create sparse symmetric matrix with 0 diagonal
int symm_matrix(const int n, const int ne, 	vector<pair <int, int> >& idx,
				vector<int>& col_ind, vector<int>& row_ptr)
{
	srand(time(0));
	
	for(int i=0; i<ne; i++){

		int rn = rand()%n;
		int cn = rand()%n;

		// no diagonal element
		while(rn==cn){
			rn=rand()%n;
			cn=rand()%n;
		}

		// no duplications
		int flag=1;
		for(int j=0; j<idx.size(); j++){
			if(idx[j].first==rn && idx[j].second==cn){
				flag=0;
				break;
			}
		}

		if(flag){
			pair<int, int> el_1(rn, cn);
			pair<int, int> el_2(cn, rn);
			idx.push_back(el_1);
			idx.push_back(el_2);
		}
	}

	sort(idx.begin(), idx.end(), comp_pairs);

	for(int i=0; i<idx.size(); i++)
		cout<<idx[i].first<<" "<<idx[i].second<<endl;

	col_ind.resize(idx.size(), 0);
	row_ptr.resize(idx[0].first+1,0);

	for(int i=0; i<idx.size()-1; i++){
		col_ind[i] = idx[i].second;
		if(idx[i].first<idx[i+1].first){
			int n_gap = idx[i+1].first-idx[i].first;
			for(int j=0; j<n_gap; j++)
				row_ptr.push_back(i+1);
		}
	}
	col_ind[idx.size()-1] = idx[idx.size()-1].second;

	// for the case that there are empty rows at the end of matrix
	for(int i=row_ptr.size(); i<n+1; i++)
		row_ptr.push_back(idx.size());
	
	// cout<<"col_ind"<<endl;
	// for(int i=0; i<col_ind.size(); i++)
	// 	cout<<col_ind[i]<<" ";
	// cout<<endl;

	// cout<<"row_ptr: "<<row_ptr.size()<<endl;
	// for(int i=0; i<row_ptr.size(); i++)
	// 	cout<<row_ptr[i]<<" ";
	// cout<<endl;
	
	return 0;
}

// E: sparse matrix to describe graph connectivity
// V: vertices
// I: independent set
// n: number of vertices
int mis_shared( const vector<int>& col_ind,
				const vector<int>& row_ptr,
				const int* V, int* I, const int n,
				const int nt )
{
	int* C = new int[n];
	int* r = new int[n];

	srand(time(0));

	// copy vertices
	for (int i=0;i<n;i++){
		C[i]=V[i];
		r[i] = rand()%(2*n);//C[i];
	}
	
	// for(int i=0; i<n; i++){
	// 	cout<<r[i]<<endl;
	// }

	// int r[6]={6,11,9,7,8,8};
	
	// number of elements searched
	int n_done=0;

	// number of independent node
	int n_is=0;

	
	while(n_done<n){		
		// node that are sweeped after each loop
		vector<int> n_sweep[nt];

#pragma omp parallel for shared(I,C,r,n_done,row_ptr,col_ind,n_is, n_sweep) num_threads(nt)
		for(int i=0; i<n; i++){
			
			if(C[i]==-1)
				continue;
	
			// if( omp_get_thread_num()==0)
			// cout<<"working on node "<<i<<" n_done "<<n_done<<endl;		

			const int myrank = omp_get_thread_num();
			
			// get neibhors
			int n_nb = row_ptr[i+1]-row_ptr[i];
			int* nb = new int[n_nb];
			int k=0;
			for(int j=0; j<n_nb; j++){
				int i_nb =  col_ind[j+row_ptr[i]];
				if(C[i_nb]!=-1){
					nb[k] = i_nb;
					k++;
				}
			}
			// number of neighbors
			n_nb=k;
			
			int fl=1;
			for(int j=0; j<n_nb; j++){
				if(r[i] > r[nb[j]]){
					fl = 0;
					break;
				}
				else if(r[i] == r[nb[j]] && i>nb[j]){
					fl = 0;
					break;
				}
			}

			// if the value is smaller than its neibhors
			// or the values are the same but node number is smaller, 
			if(fl){
				I[n_is]=i;
				n_is++;
			
				// C[i] = -1;
				for( int j=0; j<n_nb; j++){
					n_sweep[myrank].push_back(nb[j]);
				}
					// C[nb[j]]=-1;
				n_sweep[myrank].push_back(i);

				//	#pragma omp atomic
				n_done += 1+n_nb;
			}
			
			// if(omp_get_thread_num()==0){
			// for(int j=0; j<n; j++)
			// 	cout<<C[j]<<" ";
			// cout<<endl;
			// }

			delete[] nb;
			
		} //end for loop

		// sweep out all the nodes checked
		for(int j=0; j<nt; j++){
			for(int i=0; i<n_sweep[j].size(); i++){
				C[n_sweep[j][i]]=-1;
			}
		}
	} // while loop

	return n_is;
	
}


// main
int main(int argc, char **argv)
{
	// number of vertices
	const int n = 10;
	// number of edges
	const int ne = 20;
	
	vector<pair <int, int> > idx;
	vector<int> col_ind;
	vector<int> row_ptr;
	symm_matrix(n, ne, idx, col_ind, row_ptr);
	
	// int col_ind_[32] = {3, 4, 9, 4, 5, 9, 5, 7, 0, 4, 5, 0, 1, 3, 8, 1, 2, 3, 6, 8, 5, 7, 8, 2, 6, 8, 4, 5, 6, 7, 0, 1};
	// int row_ptr_[11] = {0, 3, 6, 8, 11, 15, 20, 23, 26, 30, 32};

	// col_ind.assign(col_ind_, col_ind_+32);
	// row_ptr.assign(row_ptr_, row_ptr_+11);

	// initialize vertices
	int V[n];	
	for(int i=0; i<n; i++)
		V[i] = i;

	// initialize independent set 
	int I[n];

	// number of threads
	const int nt=4;
	
	int n_is = mis_shared( col_ind, row_ptr,
						   V, I, n, nt);

	cout<<endl<<"result"<<endl;
	for(int i=0; i<n_is; i++)
		cout<<I[i]<<endl;
	
	return 0;
  
}
