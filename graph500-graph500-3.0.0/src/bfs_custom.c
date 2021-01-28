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
#include <stdbool.h>

#include <time.h>

//VISITED bitmap parameters
unsigned long *visited;
int64_t visited_size;

int64_t *pred_glob,*column;
int *rowstarts;
oned_csr_graph g;
int* q1, *q2;									// local indices
int* send_buf_visitors, * send_buf_visited;
int* recv_buf_visitors, * recv_buf_visited;
int* send_buf, * recv_buf;
int* num_visited_this_round;
int* num_visitors_this_round;
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
	
//	send_buf_visitors = xmalloc(nglobalverts_fixed*sizeof(int));
//	send_buf_visited = xmalloc(nglobalverts_fixed*sizeof(int));
//	recv_buf_visitors = xmalloc(nglobalverts_fixed*sizeof(int));
//	recv_buf_visited = xmalloc(nglobalverts_fixed*sizeof(int));

    send_buf = xmalloc(2*nglobalverts_fixed*sizeof(int));
    recv_buf = xmalloc(2*nglobalverts_fixed*sizeof(int));
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

   // printf("1\n");

    // alloc iterators
  unsigned int i, j;
 
	pred_glob=pred;

	CLEAN_VISITED();
	for(i=0;i<g.nlocalverts;i++) q1[i]=0,q2[i]=0;


    // MPI SETUP
  int my_rank = 0;
  int num_procs = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Status status;
  MPI_Status* statuss = xmalloc(2*num_procs*sizeof(MPI_Status));
  MPI_Request* requests = xmalloc(2*num_procs*sizeof(MPI_Request));
	int request_index = 0;

  if (VERTEX_OWNER(root) == my_rank) {
  	pred[VERTEX_LOCAL(root)] = root;
		SET_VISITED(root);
		q1[0] = VERTEX_LOCAL(root);
		q1c = 1;
		nvisited = 1;
  }


    sum_newly_visited = 1;
  verts_per_proc = nglobalverts_fixed / num_procs;
	num_visited_this_round = xmalloc(num_procs*sizeof(int));
	num_visitors_this_round = xmalloc(num_procs*sizeof(int));


	for (i=0;i<num_procs;i++) num_visited_this_round[i]=0, num_visitors_this_round[i]=0;
//  for (i=0;i<nglobalverts_fixed;i++) send_buf_visited[i]=-1, send_buf_visitors[i]=-1;
//  for (i=0;i<nglobalverts_fixed;i++) recv_buf_visited[i]=-1 , recv_buf_visitors[i]=-1;

    for (i=0;i<2*nglobalverts_fixed;i++) send_buf[i]=-1;
    for (i=0;i<2*nglobalverts_fixed;i++) recv_buf[i]=-1;


	// While there are vertices in current level
	
	int num_round = 0;
	
	while(sum_newly_visited != 0) {

		num_round++;
		
		q2c=0;


		for( i = 0; i < q1c; i++ ) {		  
			for( j = rowstarts[q1[i]]; j < rowstarts[q1[i]+1]; j++ ) {
				//send_buf[2 * VERTEX_OWNER( COLUMN(j))*verts_per_proc + VERTEX_LOCAL(COLUMN(j))] = q1[i];
				//printf("asdf %d\n", VERTEX_OWNER(COLUMN(j))*verts_per_proc + num_visited_this_round[VERTEX_OWNER(COLUMN(j))]);
				
				//printf("visitor: %ld, visited: %ld\n", q1[i], VERTEX_LOCAL(COLUMN(j)));
			
				bool flag = true;
			
				for (int k = 0; k < 2* num_visited_this_round[VERTEX_OWNER(COLUMN(j))]; k+=2) {
					//if (send_buf_visited[VERTEX_OWNER(COLUMN(j))*verts_per_proc + k] == VERTEX_LOCAL(COLUMN(j))) {
                    if (send_buf[VERTEX_OWNER(COLUMN(j))*2*verts_per_proc + k] == VERTEX_LOCAL(COLUMN(j))) {
						flag = false;
						break;
					}
				}
				
				if (flag) {
                    send_buf[ VERTEX_OWNER(COLUMN(j))*2*verts_per_proc + 2*num_visited_this_round[VERTEX_OWNER(COLUMN(j))]] = VERTEX_LOCAL(COLUMN(j));
					//send_buf_visitors[VERTEX_OWNER(COLUMN(j))*verts_per_proc + num_visited_this_round[VERTEX_OWNER(COLUMN(j))]] = q1[i];

					//store the visitors
                    send_buf[(VERTEX_OWNER(COLUMN(j)))*2*verts_per_proc + 2*num_visited_this_round[VERTEX_OWNER(COLUMN(j))] +1] = q1[i];

                    num_visited_this_round[VERTEX_OWNER(COLUMN(j))]++;
				}
			}
		}

//        printf("Rank: %d - send_buf:", my_rank);
//        for (i = 0; i < 2* nglobalverts_fixed; i++) {
//            printf(" %d", send_buf[i]);
//        }
////        printf("\nrecv_buf:", my_rank);
////        for (i = 0; i < 2* nglobalverts_fixed; i++) {
////            printf(" %d", recv_buf[i]);
////        }
//        printf("\n");

    //clock_t start =  clock();

    //MPI_Alltoall(num_visited_this_round, 1, MPI_INT, num_visitors_this_round, 1, MPI_INT, MPI_COMM_WORLD);

    for (int i = 1; i< num_procs +1; i++){
        prev = (my_rank-i+num_procs) % num_procs;
        next = (my_rank+i) % num_procs;
        MPI_Sendrecv(&num_visited_this_round[next], 1, MPI_INT, next, 0, &num_visitors_this_round[prev], 1, MPI_INT, prev, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    //double time = (clock() - start ) * 1000 / CLOCKS_PER_SEC;
    //if (my_rank < 5) printf("Alltoall time : %f \n", time);

        //printf("Rank %d: visited: %d\n", my_rank, num_visited_this_round[0]);
		//printf("Rank %d: visitors: %d\n", my_rank, num_visitors_this_round[0]);

		//start = clock();

        for (i = 1; i < num_procs+1; i++) {
			prev = (my_rank-i+num_procs) % num_procs;
			next = (my_rank+i) % num_procs;

            MPI_Isend(&send_buf[next*2*verts_per_proc], 2*num_visited_this_round[next], MPI_INT, next, 0, MPI_COMM_WORLD, &requests[request_index++]);
			MPI_Irecv(&recv_buf[prev*2*verts_per_proc], 2*num_visitors_this_round[prev], MPI_INT, prev, MPI_ANY_TAG, MPI_COMM_WORLD, &requests[request_index++]);

//			MPI_Isend(&send_buf_visited[next*verts_per_proc], num_visited_this_round[next], MPI_INT, next, 0, MPI_COMM_WORLD, &requests[request_index++]);
//			MPI_Irecv(&recv_buf_visited[prev*verts_per_proc], num_visitors_this_round[prev], MPI_INT, prev, MPI_ANY_TAG, MPI_COMM_WORLD, &requests[request_index++]);
			//MPI_Isend(&send_buf_visitors[next*verts_per_proc], num_visited_this_round[next], MPI_INT, next, 0, MPI_COMM_WORLD, &requests[request_index++]);
			//MPI_Irecv(&recv_buf_visitors[prev*verts_per_proc], num_visitors_this_round[prev], MPI_INT, prev, MPI_ANY_TAG, MPI_COMM_WORLD, &requests[request_index++]);
		}

		MPI_Waitall(2*num_procs, requests, statuss);
        //double time2 = (clock() - start) * 1000 / CLOCKS_PER_SEC;
        //if (my_rank < 5) printf("Send and recv time : %f \n", time2);

//			printf("Rank: %d - send_buf:", my_rank);
//			for (i = 0; i < 2* nglobalverts_fixed; i++) {
//				printf(" %d", send_buf[i]);
//			}
//			printf("\nRank: %d recv_buff", my_rank);
//			for (i = 0; i < 2* nglobalverts_fixed; i++) {
//				printf(" %d", recv_buf[i]);
//			}
//			printf("\n");


        for (i = 0; i < num_procs; i++) {
			for(j = 0; j < num_visitors_this_round[i]; j++) {
		
				//p = recv_buf_visitors[i*verts_per_proc + j];
                p = recv_buf[i*2*verts_per_proc + 2*j + 1 ];
		
				if (p != -1) {
//					if (!TEST_VISITEDLOC(recv_buf_visited[i*verts_per_proc + j])) {
//						SET_VISITEDLOC(recv_buf_visited[i*verts_per_proc + j]);
//						pred[recv_buf_visited[i*verts_per_proc + j]] = VERTEX_TO_GLOBAL(i,p);
//						q2[q2c++] = recv_buf_visited[i*verts_per_proc + j];
//						//printf("Rank %d: added %d to the queue!\n", my_rank, q2[q2c-1]);
//						nvisited++;
//					}
                    if (!TEST_VISITEDLOC(recv_buf[i*2*verts_per_proc + 2*j])) {
                        SET_VISITEDLOC(recv_buf[i*2*verts_per_proc + 2*j]);
                        pred[recv_buf[i*2*verts_per_proc + 2*j]] = VERTEX_TO_GLOBAL(i,p);
                        q2[q2c++] = recv_buf[i*2*verts_per_proc + 2*j];
                        //printf("Rank %d: added %d to the queue!\n", my_rank, q2[q2c-1]);
                        nvisited++;
                    }
				}
				
			}
			
			num_visited_this_round[i] = 0;
			num_visitors_this_round[i] = 0;
		}

    MPI_Allreduce(&q2c, &sum_newly_visited, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // swap queues
		q1c = q2c; int *tmp=q1; q1=q2; q2=tmp;
		request_index = 0;
		//nvisited += q1c;
		//if (num_round == 2) exit(-1);
	}

	//printf("Rank %d: num_visited: %d\n", my_rank, nvisited);
	free(requests);
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
	free(q1);
	free(q2);
	free(send_buf_visitors);
	free(send_buf_visited);
	free(recv_buf_visitors);
	free(recv_buf_visited);
	free(num_visited_this_round);
	free(num_visitors_this_round);
}

//user should change is function if distribution(and counts) of vertices is changed
size_t get_nlocalverts_for_pred(void) {
	return g.nlocalverts;
}
