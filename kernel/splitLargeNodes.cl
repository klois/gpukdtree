#include "device.h"
#include "structs.h"

/**
 * Creates two child nodes of each node in activelist. Activelist does not contain any leaves since all nodes are big nodes at this point
 * @param nodelist an array with all nodes
 * @param activelist an array with the indices in nodelist of all active nodes
 * @param nextlist an array with the indices in nodelist of all nodes created in the current iteration
 * @param sizes array with the number of active nodes at posititon 0, the total number of nodes at postition 1 and the number of smallnodes at position 2.
 * @param nActive the number of active nodes = size of activelist
 */
__kernel void splitLargeNodes(__global struct Node* nodelist, __global NodeId* activelist, __global NodeId* nextlist, __global UINT* sizes, UINT nActive)
{
	UINT activeId = get_global_id(0);
	if(activeId >= nActive)
		return;
	else
	{
		nextlist[activeId * 2] = 0;
		nextlist[activeId * 2 + 1] = 0;
	}

	NodeId parentId = activelist[activeId];
	struct Node node = nodelist[parentId];

	// get longest dimension to split node
	FLOAT sides[3] = {node.bounding_box.boxArr[1][0]-node.bounding_box.boxArr[0][0],
					  node.bounding_box.boxArr[1][1]-node.bounding_box.boxArr[0][1],
					  node.bounding_box.boxArr[1][2]-node.bounding_box.boxArr[0][2]};

	UINT dim = sides[0] > sides[1] ? (sides[0] > sides[2] ? 0 : (sides[1] > sides[2] ? 1 : 2)) : (sides[1] > sides[2] ? 1 : 2); //dim where to cut

	//create two new child nodes
	struct Node left_child, right_child;

	left_child.bounding_box = node.bounding_box;
	left_child.bounding_box.boxArr[1][dim] = node.bounding_box.boxArr[0][dim] + (node.bounding_box.boxArr[1][dim]-node.bounding_box.boxArr[0][dim])/2.0;
	left_child.level = node.level + 1;

	right_child.bounding_box = node.bounding_box;
	right_child.bounding_box.boxArr[0][dim] = left_child.bounding_box.boxArr[1][dim];
	right_child.level = node.level + 1;

	UINT leftIdx = sizes[1] + activeId * 2;
	UINT rightIdx = leftIdx + 1;

	// store newly created nodes in nodelist
	nodelist[leftIdx] = left_child;
	nodelist[rightIdx] = right_child;

	//add to the parent's next fields
	nodelist[parentId].left_child = leftIdx;
	nodelist[parentId].right_child = rightIdx;

	// add nodes to nextlist
	nextlist[activeId * 2] = leftIdx;
	nextlist[activeId * 2 + 1] = rightIdx;
//printf("%d: %f %f, left %f - %f, right %f - %f\n", dim, node.bounding_box.boxArr[0][dim], node.bounding_box.boxArr[1][dim],
//	   left_child.bounding_box.boxArr[0][dim], left_child.bounding_box.boxArr[1][dim], right_child.bounding_box.boxArr[0][dim], right_child.bounding_box.boxArr[1][dim]);
	// save splitting dimension in parent node
	nodelist[parentId].split_dim = dim;
//printf("nnodes %d, parent %d, left %d, right %d\n", sizes[1], parentId, leftIdx, rightIdx);
}
