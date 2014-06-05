#include "device.h"
#include "structs.h"

/**
 * Calculates a bounding box for each chunk
 * @param nodelist an array with all nodes
 * @param nNodes number of nodes in nodelist
 * @param chunks array holding all chunks
 * @param boxOut output array that will be populated with one bounding box per chunk
 * @param nodeOut output array that will be filled with pointers to nodes so that nodeOut[x] points to the node that is associated with boxOut[x]
 */
__kernel void chunkedBBox(__global struct Node* nodelist, __global NodeId* activelist, __global struct Particle* particles, __global struct Chunk* chunks,
		__global struct BBox* boxOut)
{
	UINT chunkId = get_group_id(0);

	UINT elemId = get_local_id(0);

	__local struct BBox buf[chunk_size];

	// create dummy bounding box as neutral element for bounding box calculation
	struct BBox dummy;

	// CUDA compatibility forces me to do this in that shitty way
	dummy.box[0].x = FMAX;
	dummy.box[0].y = FMAX;
	dummy.box[0].z = FMAX;
	dummy.box[1].x = FMIN;
	dummy.box[1].y = FMIN;
	dummy.box[1].z = FMIN;
	
	struct Chunk chunk = chunks[chunkId];
//printf("%d + %d = %d\n", chunk.start, elemId, chunkId);
//if(elemId == 0) printf("R cid: %d\n", chunkId);
	if(chunk.end == 0)
	{
		if(elemId == 0) // entire chunk is empty -> bbox should be empty
		{
//printf("ending group %d\n", chunkId);
			boxOut[chunkId] = dummy;
		}
		return;
	}

//	__local struct BBox buf[chunk_size];

	if(chunk.start + elemId  < chunk.end) {
//printf("%f %f %f\n", particles[0].pos.x, particles[0].pos.y, particles[0].pos.z);

		FLOAT3 pos = particles[chunk.start + elemId].pos;

		buf[elemId].box[0] = pos;
		buf[elemId].box[1] = pos;
	}
	else
		buf[elemId] = dummy;

	barrier(CLK_LOCAL_MEM_FENCE);

	for(UINT currsize = chunk_size >> 1; currsize > 0; currsize >>= 1)
	{
		// TODO check speed when dropping this check
		if(elemId < currsize)// && chunks[chunkId].start + elemId  < chunks[chunkId].end)
		{
			buf[elemId].box[0] = min(buf[elemId].box[0], buf[elemId+currsize].box[0]);
			buf[elemId].box[1] = max(buf[elemId].box[1], buf[elemId+currsize].box[1]);
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	// TODO add optimizations for GPUs

	if(elemId != 0) return;

//if(chunkId == 0) 
//printf("Chunk %d x %f %f, y %f - %f, z %f - %f\n", chunkId, buf[elemId].box[0].x, buf[elemId].box[1].x, 
//	   buf[elemId].box[0].y, buf[elemId].box[1].y, buf[elemId].box[0].z, buf[elemId].box[1].z);

	boxOut[chunkId] = buf[0];
//	nodeOut[chunkId] = chunks[0].parentNode;
}
