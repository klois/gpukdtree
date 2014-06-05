#include "device.h"
#include "structs.h"

//#include "scan.cl"

/**
 * Re-arrange the particle arrays, according to the computed segmented-scanned indeces. 
 * @param nodelist an array with all nodes
 * @oaram particles an array with all particles
 * @param activelist an array with the indices in nodelist of all active nodes
 * @param nextlist an array with the indices in nodelist of all nodes created in the current iteration
 * @param nActive the number of active nodes = size of activelist
 */
__kernel void sortP_postscan(__global struct Node* nodelist, __global struct Particle* particles, __global NodeId* activelist,
		__global NodeId* nextlist, UINT nActive,

		__global struct Particle* buffered_particles,
		__global int *scan_data		
	)
{
		//UINT activeId = get_global_id(0);
	UINT id = get_local_id(0);
	UINT activeId = get_group_id(0);
	if(activeId >= nActive)		{
//		printf("Error in sortP_postscan\n");
		return;
	}

	struct Node parent = nodelist[activelist[activeId]];
	NodeId leftId = parent.left_child;
	NodeId rightId = parent.right_child;

	UINT pLow = parent.particlesLow;
	UINT pHigh = parent.particlesHigh;
	UINT splitDim = parent.split_dim;

	const FLOAT split = nodelist[rightId].bounding_box.boxArr[0][splitDim];

	// calculating the pivot position (exclusive scan does not include last elemet)
//	const float ex_value = buffered_particles[pHigh-1].posArr[splitDim];
	const int ex_flag = 1;//(ex_value < split);
	const int rel_pivot = scan_data[pHigh-1] + ex_flag;
	const UINT localSize = get_local_size(0);

	// re-arranging particle in [pLow,pHigh) according to the pre-scanned segmented indeces
	for(UINT offset=pLow, block = 0; offset<pHigh; offset+=localSize, block++)
	{		
		const int rel_index = id + block * localSize;
		const int abs_index = pLow + rel_index;	
		if(abs_index >= pHigh) break;				 

		// 1. read current particle from buffer
		struct Particle temp = buffered_particles[abs_index];

		// 2. calculate new index
		UINT new_index;
		if(temp.posArr[splitDim] >= split) /// XXX Note(Biagio): maybe we can improve that by using the data vector? (it can be accessed coalesced)
			new_index = pLow + rel_pivot + ( rel_index - scan_data[pLow + rel_index]);//rel_pivot + scan_data0[pLow + rel_index];
		else			
			new_index = pLow + scan_data[pLow + rel_index];

//		int flag = (temp.posArr[splitDim] >= split);
//		new_index = pLow + (rel_pivot + ( rel_index - scan_data[pLow + rel_index])) * flag + (scan_data[pLow + rel_index]) *  (1- flag);
/*				
// debugging code
		if(new_index >= pHigh)
			printf("Error new index: %d >= %d \n", new_index, pHigh);
		if(new_index < pLow)
			printf("Error new index: %d <  %d \n", new_index, pLow);
*/

		// 3. write in a new position
		particles[new_index] = temp;
	}	

/*	
// debugging code
	barrier(CLK_LOCAL_MEM_FENCE); /// for debugging purpose, this code will be executed serially
	if(id==0){
		printf("\n pLow %d - pHigh %d - pivot %d", pLow, pHigh, pLow+rel_pivot);
		
		int bound = 10; //pHigh;
		printf("\nindex: ");
		for(int i=pLow; i<pLow+bound; i++) printf("%4d ", i);
		printf("\ndata0  ");
		for(int i=pLow; i<pLow+bound; i++) printf("%4d ", scan_data0[i]);
		printf("\ndata1  ");
		for(int i=pLow; i<pLow+bound; i++) printf("%4d ", scan_data1[i]);

		printf("\nlast elements\n");				
		printf("\nindex: ");
		for(int i=pHigh-bound; i<pHigh; i++) printf("%4d ", i);
		printf("\ndata0  ");
		for(int i=pHigh-bound; i<pHigh; i++) printf("%4d ", scan_data0[i]);
		printf("\ndata1  ");
		for(int i=pHigh-bound; i<pHigh; i++) printf("%4d ", scan_data1[i]);
		printf("\n");				

		printf("check low\n");				
		for(int i=pLow; i<pLow + rel_pivot; i++)
			if(particles[i].posArr[splitDim] > split)
				printf("Error in particles[%d] => %f with split %f \n", i, particles[i].posArr[splitDim], split);
		printf("check high\n");				
		for(int i=pLow + rel_pivot; i<pHigh; i++)
			if(particles[i].posArr[splitDim] < split)
				printf("Error in particles[%d] => %f with split %f \n", i, particles[i].posArr[splitDim], split);
		
		printf("double check\n");
		for(int i=pLow; i<pHigh; i++) {// for each input particle, is it present in the output?
			bool checkFlag = 0;
			for(int j=pLow; j<pHigh; j++){
				if(particles[i].posArr[splitDim] == particles[j].posArr[splitDim]){
					checkFlag = 1;
					break;
				}
			}
			if(!checkFlag)
				printf("Error: input particle %d is not in output\n", i);
		}
		printf("------\n");				
	}
*/

//	barrier(CLK_GLOBAL_MEM_FENCE);
	// update pLow
	if(id==0) {		
		//printf("\n pLow %d - pHigh %d - pivot %d", pLow, pHigh, pLow+rel_pivot);

		int newLow = pLow + rel_pivot;
/*
		for(UINT i = parent.particlesLow; i < newLow; ++i)
			if(particles[i].posArr[splitDim] > split)
				printf("%d: %f > %f\n", i, particles[i].posArr[splitDim], split);
		for(UINT i = newLow; i < parent.particlesHigh; ++i)
			if(particles[i].posArr[splitDim] < split)
				printf("%d: %f < %f\n", i, particles[i].posArr[splitDim], split);
*/
		//printf(" %d, ", newLow);
//		printf("(%d %d %d), ", pLow, newLow, pHigh);

		nodelist[leftId].particlesLow = parent.particlesLow;
		nodelist[leftId].particlesHigh = newLow;
		nodelist[rightId].particlesLow = newLow;
		nodelist[rightId].particlesHigh = parent.particlesHigh;	
	}
}


/**
 * Re-arrange the particle arrays, according to the computed segmented-scanned indices.
 * @param nodelist an array with all nodes
 * @param chunks array holding all chunks
 * @oaram particles an array with all particles
 * @param activelist an array with the indices in nodelist of all active nodes
 * @param nextlist an array with the indices in nodelist of all nodes created in the current iteration
 * @param nActive the number of active nodes = size of activelist
 */
__kernel void sortP_postscan_chunked(__global struct Node* nodelist, __global struct Chunk* chunks, __global struct Particle* particles, __global NodeId* activelist,
		UINT nActive,

		__global struct Particle* buffered_particles,
		__global int *scan_data
	)
{
	//UINT activeId = get_global_id(0);
	UINT id = get_local_id(0);
	UINT chunkId = get_group_id(0);

	struct Chunk chunk = chunks[chunkId];
	UINT activeId = chunk.parentNode;
	if(activeId >= nActive || chunk.end == 0) {
	//		printf("Error in sortP_postscan\n");
		return;
	}

	struct Node parent = nodelist[activelist[activeId]];
	NodeId leftId = parent.left_child;
	NodeId rightId = parent.right_child;

	UINT pLow = parent.particlesLow;
	UINT pHigh = parent.particlesHigh;
	UINT splitDim = parent.split_dim;

	const FLOAT split = nodelist[rightId].bounding_box.boxArr[0][splitDim];

	// calculating the pivot position (exclusive scan does not include last element)
//	const float ex_value = buffered_particles[pHigh-1].posArr[splitDim];
//	const int ex_flag = (ex_value < split);
	const int rel_pivot = scan_data[pHigh-1] + 1;

	// re-arranging particle in [pLow,pHigh) according to the pre-scanned segmented indices
//	for(UINT offset=pLow; offset<pHigh; offset+=get_local_size(0))
	UINT offset = chunk.start + id;
	if(offset < pHigh)
	{
		const int abs_index = offset;
		const int rel_index = abs_index - pLow;
//		if(abs_index >= pHigh) break;

		// 1. read current particle from buffer
		struct Particle temp = buffered_particles[abs_index];

		// 2. calculate new index
		UINT new_index;
		if(temp.posArr[splitDim] >= split) /// XXX Note(Biagio): maybe we can improve that by using the data vector? (it can be accessed coalesced)
			new_index = pLow + rel_pivot + ( rel_index - scan_data[pLow + rel_index]) -1;//rel_pivot + scan_data0[pLow + rel_index];
		else
			new_index = pLow + scan_data[pLow + rel_index];

//if(get_global_id(0) == 0) printf("idx: %d\n", new_index);
		// 3. write in a new position
		particles[new_index] = temp;
	}

//		barrier(CLK_GLOBAL_MEM_FENCE);
	// update pLow
	if(id==0 && parent.particlesLow == chunk.start) {
		int newLow = pLow + rel_pivot;

//if(parent.particlesLow == chunk.start)
//	printf("%d - %d - %d %d %f\n", pLow, newLow, pHigh,splitDim, split);

		nodelist[leftId].particlesLow = parent.particlesLow;
		nodelist[leftId].particlesHigh = newLow;
		nodelist[rightId].particlesLow = newLow;
		nodelist[rightId].particlesHigh = parent.particlesHigh;
	}
}
