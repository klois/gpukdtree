#include "device.h"
#define LOCAL_SIZE 256

#include "structs.h"

//#include "scan.cl"

/**
 * assignes the particles of the nodes in activelist to their child nodes
 * @param nodelist an array with all nodes
 * @oaram particles an array with all particles
 * @param activelist an array with the indices in nodelist of all active nodes
 * @param nextlist an array with the indices in nodelist of all nodes created in the current iteration
 * @param nActive the number of active nodes = size of activelist
 */
__kernel void sortParticlesToChilds(__global struct Node* nodelist, __global struct Particle* particles, __global NodeId* activelist,
		__global NodeId* nextlist, UINT nActive)
{
	UINT activeId = get_global_id(0);
	//UINT activeId = get_group_id(0);
	if(activeId >= nActive)		{
//		printf("Error in sortP\n");
		return;
	}

	struct Node parent = nodelist[activelist[activeId]];
	NodeId leftId = parent.left_child;
	NodeId rightId = parent.right_child;

	int pLow = parent.particlesLow;
	int pHigh = parent.particlesHigh - 1; // need to be int for size check
	UINT splitDim = parent.split_dim;

	FLOAT split = nodelist[rightId].bounding_box.boxArr[0][splitDim];



#if 1
	while(pLow < pHigh)
	{
		while((particles[pLow].posArr[splitDim] < split) && (pLow < pHigh))
		{
			++pLow;
//			printf("low %d\n", pLow);
		}
		while((particles[pHigh].posArr[splitDim] >= split) && (pHigh > pLow))
		{
			--pHigh;
//			printf("high %d\n", pHigh);
		}

//printf("*/& %d %f, %d %f\n", pLow, pHigh, particles[pLow].pos.arr[splitDim], particles[pHigh].pos.arr[splitDim]);
		if(pLow < pHigh)
		{
			struct Particle tmp = particles[pLow];
			particles[pLow] = particles[pHigh];
			particles[pHigh] = tmp;
		}
	}

//printf(" %d, ", pLow);

//printf("%d - %d - %d %d %f\n", _pLow, pLow, _pHigh, splitDim, split);

	nodelist[leftId].particlesLow = parent.particlesLow;
	nodelist[leftId].particlesHigh = pLow;
//	nodelist[leftId].small = pLow - parent.particlesLow <= T;

	nodelist[rightId].particlesLow = pLow;
	nodelist[rightId].particlesHigh = parent.particlesHigh;
//	nodelist[rightId].small = pLow - parent.particlesLow <= T;

/// printf("( CHECK  split: %f - first %f - last %f )", split, particles[_pLow].posArr[splitDim], particles[_pHigh-1].posArr[splitDim]);

#endif


#if 0	
	/// Note(Biagio)
	// this version uses additional local-memory size aux arrays and a reduction to quickly
	// sort particle two stes: left and right, acocrding to the current pivot 

	int id = get_local_id(0); 
	int wg = get_local_size(0); // workgroup size = block size, power of 2 (256, right?)
	
	// local memory 
	__local short leftFlag[LOCAL_SIZE];	
	__local int partialSum[LOCAL_SIZE];	
	__local struct Particle left[LOCAL_SIZE];
	__local struct Particle right[LOCAL_SIZE];

	__local int partialSum[LOCAL_SIZE];	

	
	int size = pHigh - pLow;
	int b = xxx;
	int pass = xxx;

	for(b..
		kernel__scan_block_anylength(partialSum, particles, b, size, pass);


	//for each LOCAL_SIZE elements
	for(int offset = leftId; offset < rightId; offset += LOCAL_SIZE)
	{
		int index = offset + id;

		// step 1, filling left & right
		if(particles[index].posArr[splitDim] < split) {
			left[id] = particles[index];
			partialSum[id] = leftFlag[id] = 1;
		}
		else {
			right[id] = particles[index];
			partialSum[id] = leftFlag[id] = 0;
		}		

		// step 2, reduction to know the size of left
	    for (unsigned int stride = LOCAL_SIZE; stride >= 1; stride >>= 1) {
			barrier(CLK_LOCAL_MEM_FENCE);
			if (id < stride)
				partialSum[id] += partialSum[id+stride];
		}

		// step 3, writing back to nodelist
		barrier(CLK_LOCAL_MEM_FENCE);
		int size = leftFlag[0];
		if(leftFlag[id])
			particles[index] = left[id];
		else
			particles[index] = right[id];

		pLow = offset + size;

printf(" %d, ", pLow);
		nodelist[leftId].particlesLow = parent.particlesLow;
		nodelist[leftId].particlesHigh = pLow;
		nodelist[rightId].particlesLow = pLow;
		nodelist[rightId].particlesHigh = parent.particlesHigh;
		
	}

#endif

}
