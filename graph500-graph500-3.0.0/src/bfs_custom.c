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

//user should provide this function which would be called once to do kernel 1: graph convert
void make_graph_data_structure(const tuple_graph* const tg) {
	//graph conversion, can be changed by user by replacing oned_csr.{c,h} with new graph format 
	convert_graph_to_oned_csr(tg, &g);

	column=g.column;
	visited_size = (g.nlocalverts + ulong_bits - 1) / ulong_bits;
	visited = xmalloc(visited_size*sizeof(unsigned long));
	//user code to allocate other buffers for bfs
	rowstarts=g.rowstarts;
}

//user should provide this function which would be called several times to do kernel 2: breadth first search
//pred[] should be root for root, -1 for unrechable vertices
//prior to calling run_bfs pred is set to -1 by calling clean_pred
void run_bfs(int64_t root, int64_t* pred) {

	//user code to do bfs
	
  int64_t nvisited = 0;
  
  // queue stuff
	int* q1 = xmalloc(g.nlocalverts*sizeof(int));
	int* q2 = xmalloc(g.nlocalverts*sizeof(int));
  int q1c, q2c;
 
  // alloc iterators
  unsigned int i, j;
 
	pred_glob=pred;

	CLEAN_VISITED();
	for(i=0;i<g.nlocalverts;i++) q1[i]=0,q2[i]=0;

	nvisited=1;
  pred[VERTEX_LOCAL(root)]=root;
	SET_VISITED(root);
	q1[0]=VERTEX_LOCAL(root);
	q1c = 1; 

	// While there are vertices in current level
	while(q1c != 0) {

		q2c=0;

		//for all vertices in current level send visit AMs to all neighbours
		for( i = 0; i < q1c; i++ ) {
		  //printf("Checking neighbours for %d!\n", q1[i]);
			for( j = rowstarts[q1[i]]; j < rowstarts[q1[i]+1]; j++ ) {
			  //printf("Found neighbour %d \n", COLUMN(j));
				if (!TEST_VISITED(COLUMN(j))) {
				  //printf("Not visitid\n");
					SET_VISITED(COLUMN(j));
					pred[COLUMN(j)] = VERTEX_LOCAL(q1[i]);
					q2[q2c++] = COLUMN(j);
				}
			}
		}

    // swap queues
		q1c = q2c; int *tmp=q1; q1=q2; q2=tmp;
		nvisited += q1c;

	}

  free(q1); free(q2);
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
