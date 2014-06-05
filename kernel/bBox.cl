#include "device.h"

#include "structs.h"

/**
 * Calculates a single bounding box for each node using it's chunks bounding boxes
 * @param nodelist an array with all nodes
 * @param in an array of bounding boxes, one for each chunk
 * @param nNodes number of nodes in nodelist
 */
__kernel void bBox(__global struct Node* nodelist, __global NodeId* activelist, __global struct BBox* in, UINT nNodes)
{
	UINT nodeId = get_group_id(0);

	if(nodeId >= nNodes)
		return;

	__local struct BBox buf[chunk_size];

	UINT elemId = get_local_id(0);
	UINT offset = (nodelist[activelist[nodeId]].particlesLow / min(T, chunk_size)) + nodeId; // starting index of current node

	// create dummy bounding box as neutral element for bounding box calculation
	struct BBox dummy;

	dummy.box[0].x = FMAX;
	dummy.box[0].y = FMAX;
	dummy.box[0].z = FMAX;
	dummy.box[1].x = FMIN;
	dummy.box[1].y = FMIN;
	dummy.box[1].z = FMIN;

//	__local struct BBox buf[chunk_size];

	buf[elemId] = dummy;

//UINT upperBound = (nodeId+1) < nNodes ?
//			((nodelist[activelist[nodeId+1]].particlesLow / min(T, chunk_size)) + nodeId + 1) : // starting index of the next node's chunks
//			((nodelist[activelist[nodeId]].particlesHigh / min(T, chunk_size)) + nodeId + 1);
UINT upperBound = (nodelist[activelist[nodeId]].particlesHigh / min(T, chunk_size)) + nodeId + 1;
//if(elemId == 0)
//printf("%d goes to %d or %d \n", nodeId, upperBound, ll); 

	for(UINT idx = offset + elemId; idx < upperBound; idx += chunk_size) {
//		buf[elemId] = in[idx];
//if(nodelist[activelist[nodeId]].bounding_box.box[1].x < in[idx].box[1].x)
//	printf("fucked up at %d with %f against %f\n", idx, nodelist[activelist[nodeId]].bounding_box.box[1].x, in[idx].box[1].x);
		buf[elemId].box[0] = min(buf[elemId].box[0], in[idx].box[0]);
		buf[elemId].box[1] = max(buf[elemId].box[1], in[idx].box[1]);
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	for(UINT currsize = get_local_size(0)/2; currsize > 0; currsize /= 2)
	{
		if(elemId < currsize) {
			buf[elemId].box[0] = min(buf[elemId].box[0], buf[elemId+currsize].box[0]);
			buf[elemId].box[1] = max(buf[elemId].box[1], buf[elemId+currsize].box[1]);
		}

		barrier(CLK_LOCAL_MEM_FENCE);
	}

// FIXME this newline is important for functionality on GTX480 with 310.19
	// TODO add optimizations for GPUs
	if(elemId == 0)
	{
//printf("%d: [%f %f], [%f %f], [%f %f]\n", nodeId, buf[0].boxArr[0][0], buf[0].boxArr[1][0],
//		buf[0].boxArr[0][1], buf[0].boxArr[1][1], buf[0].boxArr[0][2], buf[0].boxArr[1][2]);
		nodelist[activelist[nodeId]].bounding_box = buf[0];
	}
}

