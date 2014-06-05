#include "device.h"
#include "structs.h"

/**
 * copies all nodes on level currLevel on the coorect position in tree and calculates the position for all nodes on level currLevel+1
 * @param nodelist an array with all nodes
 * @param tree array that will be filled with nodes from nodelist in depth-first order
 * @param currLevel the level in the three that will be processed
 * @param nNodes number of nodes in nodelist
 */
__kernel void downPass(__global struct Node* nodelist, __global struct KdNode* tree, UINT currLevel, UINT nNodes)
{
	UINT idx = get_global_id(0);

	if(idx >= nNodes)
		return;

	struct Node node = nodelist[idx];
	struct KdNode kdNode;

	kdNode.center_of_mass = node.center_of_mass;
	kdNode.mass = node.mass;
	kdNode.center_geometric = node.center_geometric;
	kdNode.l = node.l;
	kdNode.size = node.size;
	kdNode.leaf = (node.particlesHigh - node.particlesLow) == 1;

	if(node.level != currLevel)
		return;

	if(node.particlesHigh - node.particlesLow > 1) { // leaf node check
		nodelist[node.left_child].address = node.address + 1;
		nodelist[node.right_child].address = node.address + 1 + nodelist[node.left_child].size;
		// store node it at FINAL format to node.address
		kdNode.left_child = node.address + 1;
		kdNode.right_child = node.address + 1 + nodelist[node.left_child].size;
	} else {
		kdNode.left_child = 0;
	}

	// store node it at FINAL format to it.address
	tree[node.address] = kdNode;

}
