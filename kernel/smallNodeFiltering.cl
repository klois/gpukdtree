#include "device.h"
#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics: enable

#include "structs.h"

/**
 * identifies all small nodes in nextlist, adds them to smalllist sets nextlist at their position to 0
 * @param nodelist an array with all nodes
 * @param nextlist an array with the indices in nodelist of all nodes created in the current iteration
 * @param smalllist an array with the indices in nodelist of all small nodes
 * @param sizes array with the number of active nodes at posititon 0, the total number of nodes at postition 1 and the number of smallnodes at position 2.
 * 		position 1 and 2 will be updated inside this kernel
 */
__kernel void smallNodeFiltering(__global struct Node* nodelist, __global NodeId* nextlist, 
								 __global NodeId* smalllist,     __global UINT* sizes
								 //, __global UINT *numParticles
								 )
{
	UINT nNext = sizes[0] * 2;
	UINT nextId = get_global_id(0);
	if(nextId >= nNext)
		return;

	//macro to debug memory accesses
#if 0
	int maxNodes = numParticles * 2 - 1;
	if(nextId >= numParticles || nextId <0)
		printf("ERROR(1) in snf for thread %d - %d\n", nextId, numParticles);
	
	if(nextlist[nextId] >= maxNodes || nextlist[nextId]<0)
		printf("ERROR(2) in snf for thread %d\n", nextlist[nextId] );	

	if(smalllist[nextId] >= maxNodes || smalllist[nextId]<0)
		printf("ERROR(3) in snf for thread %d\n", smalllist[nextId] );

#endif

	if((nodelist[nextlist[nextId]].particlesHigh - nodelist[nextlist[nextId]].particlesLow) <= T)
	{
		UINT idx = atomic_inc(&sizes[2]);
		smalllist[idx] = nextlist[nextId]; // add node to smallnodes
//printf("%d small with %d\n", nextlist[nextId], (nodelist[nextlist[nextId]].particlesHigh - nodelist[nextlist[nextId]].particlesLow));
		nextlist[nextId] = 0; // mark position in nextlist as unused
	} 
//	else
//		printf("%d large with %d\n", nextlist[nextId],(nodelist[nextlist[nextId]].particlesHigh - nodelist[nextlist[nextId]].particlesLow));

	if(nextId == 0)
		sizes[1] += nNext; // update size of nodelist
}
