#include "device.h"
#include "structs.h"

/**
 * calculates the size (in nodes) for each nodes including it's subtree
 * @param nodelist an array with all nodes
 * @param tree an array with the nodes of a previously build tree
 * @param particles an array with all particles
 * @param nNodes number of nodes in nodelist
 */
__kernel void updateTree(__global struct Node* nodelist, __global struct KdNode* tree, __global struct Particle* particles, int currLevel, UINT nNodes)
{
	UINT idx = get_global_id(0);

	if(idx >= nNodes)
		return;

	struct Node node = nodelist[idx];
	if(node.level != currLevel)
		return;

	if(node.particlesHigh - node.particlesLow <= 1)
	{
		// size, mass and l of leafs does not change	
	
		struct Particle particle = particles[node.particlesLow];
	
		node.center_of_mass = particle.pos;
		node.center_geometric = particle.pos;
		
		node.bounding_box.box[0] = particle.pos;
		node.bounding_box.box[1] = particle.pos;
	}
	else
	{
		// size and mass of nodes does not change

		struct Node leftChild = nodelist[node.left_child];
		struct Node rightChild = nodelist[node.right_child];

// this expression does not work in CUDA
//		node.center_of_mass = (leftChild.center_of_mass * leftChild.mass +
//				rightChild.center_of_mass * rightChild.mass ) / node.mass;

		// even CUDA can compile and execute this block correctly
		node.center_of_mass.x = (leftChild.center_of_mass.x * leftChild.mass +
				rightChild.center_of_mass.x * rightChild.mass ) / node.mass;
		node.center_of_mass.y = (leftChild.center_of_mass.y * leftChild.mass +
				rightChild.center_of_mass.y * rightChild.mass ) / node.mass;
		node.center_of_mass.z = (leftChild.center_of_mass.z * leftChild.mass +
				rightChild.center_of_mass.z * rightChild.mass ) / node.mass;
				
		// update bounding box
		node.bounding_box.box[0] = min(leftChild.bounding_box.box[0], rightChild.bounding_box.box[0]);
		node.bounding_box.box[1] = max(leftChild.bounding_box.box[1], rightChild.bounding_box.box[1]);

		FLOAT3 l = node.bounding_box.box[1] - node.bounding_box.box[0];
		node.center_geometric = node.bounding_box.box[0] + l / 2;
		node.l = max(l.x, max(l.y, l.z));
/*
		FLOAT lx = (node.bounding_box.box[1].x - node.bounding_box.box[0].x);
		FLOAT ly = (node.bounding_box.box[1].y - node.bounding_box.box[0].y);
		FLOAT lz = (node.bounding_box.box[1].z - node.bounding_box.box[0].z);

		node.center_geometric.x = node.bounding_box.box[0].x + lx / 2;
		node.center_geometric.y = node.bounding_box.box[0].y + ly / 2;
		node.center_geometric.z = node.bounding_box.box[0].z + lz / 2;

		//node.l = (lx + ly + lz) / 3;
		node.l = max(lx, max(ly, lz));
*/
	}

	// store updated nodes at old position
	nodelist[idx] = node;
	
	// construct kdNode
	struct KdNode kdNode = tree[node.address];

	kdNode.center_of_mass = node.center_of_mass;
	kdNode.mass = node.mass;
	kdNode.center_geometric = node.center_geometric;
	kdNode.l = node.l;
//	kdNode.size = node.size;
//	kdNode.leaf = (node.particlesHigh - node.particlesLow) == 1;

	// write node to tree, position did not change and is therefore known
	tree[node.address] = kdNode;
}
