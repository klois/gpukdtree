#include "device.h"
#include "structs.h"

#if 0
/**
 * assignes the particles of the nodes in activelist to their child nodes
 * @param nodelist an array with all nodes
 * @oaram particles an array with all particles
 * @param activelist an array with the indices in nodelist of all active nodes
 * @param nextlist an array with the indices in nodelist of all nodes created in the current iteration
 * @param nActive the number of active nodes = size of activelist
 */
__kernel void sortP_prescan(__global struct Node* nodelist, __global struct Particle* particles, __global NodeId* activelist,		
		UINT nActive,		
		__global int *scan_data,
		__global int *scan_flag				
	)
{
	//UINT activeId = get_global_id(0);
	UINT id = get_local_id(0);
	UINT activeId = get_group_id(0);
	if(activeId >= nActive)		{
//		printf("Error in sortP_prescan\n");
		return;
	}

	struct Node parent = nodelist[activelist[activeId]];
	NodeId leftId = parent.left_child;
	NodeId rightId = parent.right_child;

	if(id == 0)
	{
		printf("%d: %d - %d\n", activelist[activeId], parent.particlesLow, parent.particlesHigh);
	}

//	UINT pLow = parent.particlesLow;
//	UINT pHigh = parent.particlesHigh;
	UINT splitDim = parent.split_dim;

	const FLOAT split = nodelist[rightId].bounding_box.boxArr[0][splitDim];

	// for each particle in [pLow,pHigh)
	for(UINT index=parent.particlesLow + id; index<parent.particlesHigh; index+=get_local_size(0))
	{
		const int value = particles[index].posArr[splitDim] < split; // value will be 0 if on left, 1 if on right
//		const int notValue = 1 - value;
		scan_data[index] = value;

//		const int flag = (index==pLow); // 1 for the first element, 0 0therwise
		scan_flag[index] = (index==parent.particlesLow);
	}	

	// this avoid that the prefix scan includes also other nodes
//	if(id==0 && pHigh < nParticles) scan_flag[pHigh] = 1;	
	
}
#endif
/**
 * assignes the particles of the nodes in activelist to their child nodes
 * @param nodelist an array with all nodes
 * @param chunks array holding all chunks
 * @oaram particles an array with all particles
 * @param activelist an array with the indices in nodelist of all active nodes
 * @param nextlist an array with the indices in nodelist of all nodes created in the current iteration
 * @param nActive the number of active nodes = size of activelist
 */
__kernel void sortP_prescan_chunked(__global struct Node* nodelist, __global struct Chunk* chunks, __global struct Particle* particles,
		__global NodeId* activelist, UINT nActive,
		__global int *scan_data,
		__global int *scan_flag
	)
{
	//UINT activeId = get_global_id(0);
	UINT id = get_local_id(0);
	UINT chunkId = get_group_id(0);

	struct Chunk chunk = chunks[chunkId];
	UINT activeId = chunk.parentNode;
	if(activeId >= nActive || chunk.end == 0)		{
//		printf("Error in sortP_prescan %d \t", chunk.end);
		return;
	}

	struct Node parent = nodelist[activelist[activeId]];
//	NodeId leftId = parent.left_child;
	NodeId rightId = parent.right_child;
//printf("%d\t", activeId);
//	UINT pLow = parent.particlesLow;
//	UINT pHigh = parent.particlesHigh;
	UINT splitDim = parent.split_dim;

	const FLOAT split = nodelist[rightId].bounding_box.boxArr[0][splitDim];

	// for each particle in [pLow,pHigh)
//	for(UINT index=parent.particlesLow + id; index<parent.particlesHigh; index+=get_local_size(0))
	UINT index = chunk.start + id;
	if(index >= chunk.end) return;

	const int value = particles[index].posArr[splitDim] < split; // value will be 0 if on left, 1 if on right
	scan_data[index] = value;

	scan_flag[index] = parent.particlesLow;

	// this avoid that the prefix scan includes also other nodes
//	if(id==0 && pHigh < nParticles) scan_flag[pHigh] = 1;

}
