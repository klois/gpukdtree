#include "device.h"
#include "structs.h"

// swap Node** and its size
void swap(struct Node*** a, UINT* aN, struct Node*** b, UINT* bN)
{
	struct Node** tmp = *a;
	UINT tmpN = *aN;
	*a = *b;
	*aN = *bN;
	*b = tmp;
	*bN = tmpN;
}

/**
 * divides the particles in a node into fixed sized chunks
 * @param nodelist an array with all nodes
 * @param activelist an array with the indices in nodelist of all nodes to be processed
 * @param chunks array which will be filled with the generated chunks
 * @param nActive the number of active nodes = size of activelist
 */
__kernel void groupParticlesIntoChunks(__global struct Node* nodelist, __global NodeId* activelist, __global struct Chunk* chunks, UINT activeN) {
//for(int i = 0; i < 55; ++i)
//	printf("%d: %d\n", i, check[i]);
	UINT nodeId = get_group_id(0);

	struct Node node = nodelist[activelist[nodeId]];

	UINT chunkIdT = get_local_id(0);
	UINT nParticles = node.particlesHigh - node.particlesLow;
	UINT offset =  (node.particlesLow / min(T, chunk_size)) + nodeId; // starting index of current node's chunks
	UINT nextOffset = (nodeId+1) < activeN ?
			((nodelist[activelist[nodeId+1]].particlesLow / min(T, chunk_size)) + nodeId + 1) : // starting index of the next node's chunks
			((node.particlesHigh / min(T, chunk_size)) + nodeId + 1);
	UINT nChunks = (nParticles + chunk_size - 1) / chunk_size;

	for(UINT chunkId = chunkIdT; offset + chunkId < nextOffset; chunkId += get_local_size(0))
	{
		if(chunkId >= nChunks)
			return;
	//printf("W cid: %d\n", offset + chunkId);

		UINT lowerBound = node.particlesLow;
	//	printf("node %d\n", offset + chunkId);
		// create chunk
		struct Chunk newChunk;
		newChunk.parentNode = nodeId;

//printf("%d\t", nodeId);
		newChunk.start = lowerBound + chunkId * chunk_size;
	//	newChunk.end = (lowerBound + (chunkId+1) * chunk_size > node.particlesHigh) ?
	//			node.particlesHigh : lowerBound + (chunkId+1) * chunk_size;

		UINT flag = ((lowerBound + (chunkId+1) * chunk_size > node.particlesHigh));
		newChunk.end = (node.particlesHigh * flag) + ((lowerBound + (chunkId+1) * chunk_size) * (1 - flag));
		chunks[offset + chunkId] = newChunk;
	}
//if(activelist[chunks[offset + chunkId].parentNode] == 52)
//printf("NodeId: %d, %d : %d\n", activelist[nodeId], chunks[offset + chunkId].start, chunks[offset + chunkId].end);
}
