
// p: proessor id
void MxM_MPI(E, V, C, nc, M)
{
	for(all_u_in_V) M(u) = -1; // tc*n

	for (k=0; k<nc-1; k++){ 

		W.empty();
		
		// D(p) is the set of the node in processor p
		// nc = max(#neighbor)+1
		// O(n/p/nc) loops * tc*O(#neighbor) = tc*O(n/p)
		for ( u_st_C(u)==k && M(u)==-1 && u_in_D(p)){ // n/p/nc loops
			NE = E(u,:); // get neighbor
			NE = {v_in_NE && M(v)==-1} // tc*(#neighbor)

			M(u) = v;
			M(v) = u;

			W.push(u);
		}
		
		// (ts+tw*m)*(p-1) for all2all
		// m = n/p
		MPI_Alltoall(M, n/p, M_rec, n/p); 

		// O(n/p)
		for(u_in_D(p)){
			M(u) = i;
			pi = P(i); // get the processor id corresponding to i
			if( M_rec(i)!=-1 )
				if( M(i)==-1 || p<pi ) // higher processor id is stronger
					M(i) = M_rec(i);
		}

		// O(nc/p)
		for( u_in_W ){
			if( M(M(u)) != -1 ) M(u) = -1;
		}
	}
}


// T = nc * (tc*n/p + (ts+tw*n/p)*(p-1))
// nc can be approximated by the average number of neighboring nodes (how dense the graph is)


	
