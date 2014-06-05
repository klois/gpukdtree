#include "device.h"
#include "structs.h"

/**
 * writes initial values in the sizes and activelist array
 * @param activelist array holding the indices in nodelist of the active nodes. For first iteration the first element will be set to 0: the root is active
 * @param sizes array continging
 * 		0 activelistsize: will be set to 1 (only the root is active)
 * 		1 total number of nodes: will be set to 1 (the root is the only node)
 * 		2 smalllistsize: will be set to 0 (no smallnodes identified so far)
 */
__kernel void init(__global NodeId* activelist, __global UINT* sizes, __global struct Particle* particles)
{
//printf("BASE %d, particlesLow %f, offset %d\n", particles, &particles[0].totAcc, (&(particles[0].totAcc) - particles));

//printf("x %f %f, y %f - %f, z %f - %f\n", nodelist[0].bounding_box.box[0].x, nodelist[0].bounding_box.box[1].x, 
//   nodelist[0].bounding_box.box[0].y, nodelist[0].bounding_box.box[1].y, nodelist[0].bounding_box.box[0].z, nodelist[0].bounding_box.box[1].z);

	// add root node to activelist
	activelist[0] = 0;

	// set size of activelist to 1
	sizes[0] = 1;

	// set size of nodelist to 1
	sizes[1] = 1;

	// set size of smalllist to 0
	sizes[2] = 0;	

	// set three depth
	sizes[3] = 1;
	sizes[4] = 1;
}


__kernel void memset_int(__global int *mem, __private int val) { 
	mem[get_global_id(0)] = val; 
}

__kernel void memset_int_s(__global int *mem, __private int val, UINT size) {//, __global int *scan) {
	const int id = get_global_id(0);
	if(id < size) {
		mem[id] = val;
//		scan[id] = -10000000;
	}
}

__kernel void memset_chunks(__global struct Chunk* chunks, UINT size) {
	const int id = get_global_id(0);
	if(id < size) {
		chunks[id].end = 0;
	}
}
