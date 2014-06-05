#include "device.h"
#include "structs.h"

/**
 * Creates two child nodes of each non-leave node in activelist and assignes the node's particles to the newly created child nodes
 * @param nodelist an array with all nodes
 * @param activelist an array with the indices in nodelist of all active nodes
 * @param nextlist an array with the indices in nodelist of all nodes created in the current iteration
 * @oaram particles an array with all particles
 * @param sizes array with the number of active nodes at posititon 0, the total number of nodes at postition 1 and the number of newly created nodes in this
 * 		iteration at position 2, the current maximum level of the tree at position 3 and 4. Position 4 will be increased by 1 if the maximum level grows in
 * 		this iteration.
 */
__kernel void splitSmallNodes(__global struct Node* nodelist, __global NodeId* activelist, __global NodeId* nextlist, __global struct Particle* particles,
		__global UINT* sizes)
{
	UINT activeId = get_global_id(0);
	if(activeId >= sizes[0]) 
		return;

	NodeId parentId = activelist[activeId];

	struct Node node = nodelist[parentId];

	// get longest dimension to split node
	FLOAT sides[3] = {node.bounding_box.boxArr[1][0]-node.bounding_box.boxArr[0][0],
					  node.bounding_box.boxArr[1][1]-node.bounding_box.boxArr[0][1],
					  node.bounding_box.boxArr[1][2]-node.bounding_box.boxArr[0][2]};
	UINT d = sides[0] > sides[1] ? (sides[0] > sides[2] ? 0 : (sides[1] > sides[2] ? 1 : 2)) : (sides[1] > sides[2] ? 1 : 2); //dim where to cut
	FLOAT V0 = sides[d];//sides[0] * sides[1] * sides[2];

	// compute SAH and determine the split plane
	FLOAT SAH = FMAX; // always split until only one particle per node left
	FLOAT ppos = 17.0f;
//	UINT leftmass; // number of particles in child to calculate mass
//	UINT rightmass;
	// check for leafs
	if(node.particlesHigh - node.particlesLow <= 1u) {
		return;
	}

	//create two new child nodes
	struct Node left_child, right_child;

	// sort particles to child nodes
	int pLow = node.particlesLow;
	int pHigh = node.particlesHigh - 1;

	// special handling for nodes with only two particles inside
	if(pHigh - pLow == 1) {
		if(particles[pLow].posArr[d] > particles[pHigh].posArr[d])
		{
			struct Particle tmp = particles[pLow];
			particles[pLow] = particles[pHigh];
			particles[pHigh] = tmp;
		}
		++pLow;
	} else {
		// loops over splitting points
		for(UINT s = node.particlesLow; s < node.particlesHigh; ++s)
		{
	//		for(UINT d = 0; d < 3; ++d) // every particle yields to a splitting point
			{
				// number of nodes in potential left an right child
				UINT Cl = 0;
				UINT Cr = 0;

				for(UINT p = node.particlesLow; p < node.particlesHigh; ++p)
				{
					if(particles[p].posArr[d] < particles[s].posArr[d])
						++Cl;
					else
						++Cr;
				}

				// all nodes have to contain particles after a split
				if((Cl * Cr) == 0)
					continue;

				// volume of potential left and right child
				FLOAT Vl = particles[s].posArr[d] - node.bounding_box.boxArr[0][d];
				/*V0 * ( // volume of parent node
						(particles[s].posArr[d] - node.bounding_box.boxArr[0][d]) / // side length of child
						sides[d] ); // side length of parent */
				FLOAT Vr = V0 - Vl;

				FLOAT SAHj = (Cl*Vl + Cr*Vr);

				if(SAH > SAHj) {
					SAH = SAHj;
					ppos = particles[s].posArr[d];
	//				leftmass = Cl;
	//				rightmass = Cr;
				}
			}
		}

	//	if(pHigh -pLow == 2)
	//		printf("A %f: %f %f %f \n", ppos, particles[pLow].posArr[pdim], particles[pLow+1].posArr[pdim], particles[pLow+2].posArr[pdim]);
		while(pLow < pHigh)
		{
			while((particles[pLow].posArr[d] < ppos) && (pLow < pHigh))
			{
				++pLow;
			}
			while((particles[pHigh].posArr[d] >= ppos) && (pHigh > pLow))
			{
				--pHigh;
			}

			if(pLow < pHigh)
			{
				struct Particle tmp = particles[pLow];
				particles[pLow] = particles[pHigh];
				particles[pHigh] = tmp;
			}
		}
	}

	left_child.particlesLow = node.particlesLow;
	left_child.particlesHigh = pLow;
//	nodelist[leftId].small = pLow - parent.particlesLow <= T;

	right_child.particlesLow = pLow;
	right_child.particlesHigh = node.particlesHigh;
//	nodelist[rightId].small = pLow - parent.particlesLow <= T;


	left_child.level = node.level + 1;
	right_child.level = node.level + 1;

//#define THIGHT
#ifndef THIGHT
	left_child.bounding_box = node.bounding_box;
	left_child.bounding_box.boxArr[1][d] = ppos;
	right_child.bounding_box = node.bounding_box;
	right_child.bounding_box.boxArr[0][d] = ppos;

#else
	struct BBox right_box;
	struct BBox left_box;

	right_box.box[0] = particles[right_child.particlesLow].pos;
	right_box.box[1] = particles[right_child.particlesLow].pos;

	for(UINT s = right_child.particlesLow + 1; s < right_child.particlesHigh; ++s)
	{
		right_box.box[0] = min(right_box.box[0], particles[s].pos);
		right_box.box[1] = max(right_box.box[1], particles[s].pos);
	}

	left_box.box[0] = particles[left_child.particlesLow].pos;
	left_box.box[1] = particles[left_child.particlesLow].pos;

	for(UINT s = left_child.particlesLow + 1; s < left_child.particlesHigh; ++s)
	{
		left_box.box[0] = min(left_box.box[0], particles[s].pos);
		left_box.box[1] = max(left_box.box[1], particles[s].pos);
	}

	left_child.bounding_box = node.bounding_box;
	left_child.bounding_box.boxArr[1][d] = ppos;
	right_child.bounding_box = node.bounding_box;
	right_child.bounding_box.boxArr[0][d] = ppos;

	FLOAT3 left_diffs = fabs((left_child.bounding_box.box[1] - left_child.bounding_box.box[0]) - (left_box.box[1] - left_box.box[0]));
	FLOAT3 right_diffs = fabs((right_child.bounding_box.box[1] - right_child.bounding_box.box[0]) - (right_box.box[1] - right_box.box[0]));

	FLOAT left_diff = max(max(left_diffs.x, left_diffs.y), left_diffs.z);
	FLOAT right_diff = max(max(right_diffs.x, right_diffs.y), right_diffs.z);

	FLOAT threshold = 100.0f;

	if(left_diff > threshold)
		left_child.bounding_box = left_box;

	if(right_diff > threshold)
		right_child.bounding_box = right_box;

#endif

	// save splitting dimension in parent node
	nodelist[parentId].split_dim = d;

	// add nodes to nextlist
	UINT childId = atomic_add(&sizes[2], 2u);
	UINT leftIdx = sizes[1] + childId;
	UINT rightIdx = leftIdx + 1;
//if(childId >= 250000) printf("child is gone %d\n", childId);
	nextlist[childId] = leftIdx;
	nextlist[childId + 1] = rightIdx;
//if(rightIdx >= 500000) printf("right is gone %d\n", rightIdx);

	//add to the parent's next fields
	nodelist[parentId].left_child = leftIdx;
	nodelist[parentId].right_child = rightIdx;

	// save newly created nodes in nodelsit
	nodelist[leftIdx] = left_child;
	nodelist[rightIdx]  = right_child;

	//increase max level if necessary
	if(node.level == sizes[3]) 
		sizes[4] = node.level + 1;

//	nodelist[parentId].center_of_mass[3] = SAH;

//printf("nnodes %d, parent %d, left %d, right %d\n", sizes[1], parentId, leftIdx, rightIdx);
}

