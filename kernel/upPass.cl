#include "device.h"
#include "structs.h"

/**
 * calculates the size (in nodes) for each nodes including it's subtree
 * @param nodelist an array with all nodes
 * @param particles an array with all particles
 * @param nNodes number of nodes in nodelist
 */
__kernel void upPass(__global struct Node* nodelist, __global struct Particle* particles, int currLevel, UINT nNodes)
{
	UINT idx = get_global_id(0);

	if(idx >= nNodes)
		return;

	struct Node node = nodelist[idx];
	if(node.level != currLevel)
		return;

	if(node.particlesHigh - node.particlesLow <= 1)
	{
		node.size = 1;
		node.mass = particles[node.particlesLow].mass;
//printf("%d leaf\n", node.level);
		node.center_of_mass = particles[node.particlesLow].pos;
		node.center_geometric = particles[node.particlesLow].pos;
/*
		node.center_of_mass.x = particles[node.particlesLow].pos.x;
		node.center_of_mass.y = particles[node.particlesLow].pos.y;
		node.center_of_mass.z = particles[node.particlesLow].pos.z;

		node.center_geometric.x = particles[node.particlesLow].pos.x;
		node.center_geometric.y = particles[node.particlesLow].pos.y;
		node.center_geometric.z = particles[node.particlesLow].pos.z;
*/
		node.l = 0;
	}
	else
	{
		node.size = nodelist[node.left_child].size + nodelist[node.right_child].size + 1;
//printf("%d %d\n", node.level, node.size);
		node.mass = nodelist[node.left_child].mass + nodelist[node.right_child].mass;

// this expression does not work in CUDA
//		node.center_of_mass = (nodelist[node.left_child].center_of_mass * nodelist[node.left_child].mass +
//				nodelist[node.right_child].center_of_mass * nodelist[node.right_child].mass ) / node.mass;

		// even CUDA can compile and execute this block correctly
		node.center_of_mass.x = (nodelist[node.left_child].center_of_mass.x * nodelist[node.left_child].mass +
				nodelist[node.right_child].center_of_mass.x * nodelist[node.right_child].mass ) / node.mass;
		node.center_of_mass.y = (nodelist[node.left_child].center_of_mass.y * nodelist[node.left_child].mass +
				nodelist[node.right_child].center_of_mass.y * nodelist[node.right_child].mass ) / node.mass;
		node.center_of_mass.z = (nodelist[node.left_child].center_of_mass.z * nodelist[node.left_child].mass +
				nodelist[node.right_child].center_of_mass.z * nodelist[node.right_child].mass ) / node.mass;

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
}
