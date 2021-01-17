//Stub for custom BFS implementations

#include "common.h"
#include "aml.h"
#include "csr_reference.h"
#include "bitmap_reference.h"
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <stdint.h>

//VISITED bitmap parameters
unsigned long *visited;
int64_t visited_size;

int64_t *pred_glob,*column;
int *rowstarts;
oned_csr_graph g;
int* q1, *q2;									// local indices
int64_t* send_buf;
int64_t* recv_buf;
int nglobalverts_fixed;

//user should provide this function which would be called once to do kernel 1: graph convert
void make_graph_data_structure(const tuple_graph* const tg) {
	//graph conversion, can be changed by user by replacing oned_csr.{c,h} with new graph format 
	convert_graph_to_oned_csr(tg, &g);

	column=g.column;
	visited_size = (g.nlocalverts + ulong_bits - 1) / ulong_bits;
	visited = xmalloc(visited_size*sizeof(unsigned long));
	//user code to allocate other buffers for bfs
	rowstarts=g.rowstarts;
  // queue stuff
	q1 = xmalloc(g.nlocalverts*sizeof(int));
	q2 = xmalloc(g.nlocalverts*sizeof(int));
	nglobalverts_fixed = g.nglobalverts - g.nglobalverts%2;
	send_buf = xmalloc(nglobalverts_fixed*sizeof(int64_t));
	recv_buf = xmalloc(nglobalverts_fixed*sizeof(int64_t));

}

//user should provide this function which would be called several times to do kernel 2: breadth first search
//pred[] should be root for root, -1 for unrechable vertices
//prior to calling run_bfs pred is set to -1 by calling clean_pred
void run_bfs(int64_t root, int64_t* pred) {

	//user code to do bfs
	
  int64_t nvisited = 0;
  int q1c = 0, q2c = 0;
  int sum_newly_visited = 0;
  int verts_per_proc = 0;
  int64_t p = 0;
  int prev, next;
 
  // alloc iterators
  unsigned int i, j;
 
	pred_glob=pred;

	CLEAN_VISITED();
	for(i=0;i<g.nlocalverts;i++) q1[i]=0,q2[i]=0;
  for(i=0;i<nglobalverts_fixed;i++) send_buf[i]=-1;

  // MPI SETUP
  int my_rank = 0;
  int num_procs = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Status status;

  if (VERTEX_OWNER(root) == my_rank) {
  	pred[VERTEX_LOCAL(root)] = root;
		SET_VISITED(root);
		q1[0] = VERTEX_LOCAL(root);
		q1c = 1;
		nvisited = 1;
  }
  
  sum_newly_visited = 1;
  verts_per_proc = nglobalverts_fixed / num_procs;
  
	// While there are vertices in current level
	
	int num_round = 0;
	
	while(sum_newly_visited != 0) {

		num_round++;
		
		q2c=0;

		for( i = 0; i < q1c; i++ ) {		  
			for( j = rowstarts[q1[i]]; j < rowstarts[q1[i]+1]; j++ ) {
				send_buf[VERTEX_OWNER(COLUMN(j))*verts_per_proc + VERTEX_LOCAL(COLUMN(j))] = q1[i];				
			}
		}


    //MPI_Alltoall(send_buf, verts_per_proc, MPI_LONG, recv_buf, verts_per_proc, MPI_LONG, MPI_COMM_WORLD);
		
		/*		
		if (my_rank == 0) {
			printf("Rank: %d - send_buf:", my_rank);
			for (i = 0; i < g.nglobalverts-1; i++) {
				printf(" %ld", send_buf[i]);
			}
			printf("\n");
		}
		*/
		
		for (i = 1; i < num_procs+1; i++) {
			prev = (my_rank-i+num_procs) % size;
			next = (my_rank+i) % size;
			MPI_Sendrecv(&send_buf[next*verts_per_proc], verts_per_proc, MPI_LONG, next, 0, &recv_buf[prev*verts_per_proc], verts_per_proc, MPI_LONG, prev, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		}
		
		/*
		if (my_rank == 1) {
			printf("Rank: %d - recv_buf:", my_rank);
			for (i = 0; i < g.nglobalverts-1; i++) {
				printf(" %ld", recv_buf[i]);
			}
			printf("\n");
		}
		*/
		
		for (i = 0; i < num_procs; i++) {
			for(j = 0; j < verts_per_proc; j++) {
		
				p = recv_buf[i*verts_per_proc + j];
		
				if (p != -1) {
					if (!TEST_VISITEDLOC(j)) {
						SET_VISITEDLOC(j);
						pred[j] = VERTEX_TO_GLOBAL(i,p);
						q2[q2c++] = j;
					}
				}
				
				send_buf[i*verts_per_proc + j] = -1;
			}
		}

    MPI_Allreduce(&q2c, &sum_newly_visited, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // swap queues
		q1c = q2c; int *tmp=q1; q1=q2; q2=tmp;
		nvisited += q1c;
		
		if (nvisited == verts_per_proc && q1c == 0) {
			printf("Round %d: %d is dead! x.x\n", num_round, my_rank);
		}
	}
	
	//printf("Rank %d: num_visited: %d\n", my_rank, nvisited);
}

//we need edge count to calculate teps. Validation will check if this count is correct
//user should change this function if another format (not standart CRS) used
void get_edge_count_for_teps(int64_t* edge_visit_count) {
	long i,j;
	long edge_count=0;
	for(i=0;i<g.nlocalverts;i++)
		if(pred_glob[i]!=-1) {
			for(j=g.rowstarts[i];j<g.rowstarts[i+1];j++)
				if(COLUMN(j)<=VERTEX_TO_GLOBAL(my_pe(),i))
					edge_count++;
		}
	aml_long_allsum(&edge_count);
	*edge_visit_count=edge_count;
}

//user provided function to initialize predecessor array to whatevere value user needs
void clean_pred(int64_t* pred) {
	int i;
	for(i=0;i<g.nlocalverts;i++) pred[i]=-1;
}

//user provided function to be called once graph is no longer needed
void free_graph_data_structure(void) {
	free_oned_csr_graph(&g);
	free(visited);
}

//user should change is function if distribution(and counts) of vertices is changed
size_t get_nlocalverts_for_pred(void) {
	return g.nlocalverts;
}
