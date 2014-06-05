#include "device.h"

// I'm sick of this
#ifndef __device__
#define __device__
#endif

#include "structs.h"

#define SOFTENING_SPLINEKERNEL
//#define VECTOR

__device__ void dvec(const FLOAT3 a, const FLOAT3 b, FLOAT res[3])
{
	res[0] = b.x - a.x;
	res[1] = b.y - a.y;
	res[2] = b.z - a.z;
}


__device__ FLOAT3 calcAcc(const FLOAT3 ppos, const FLOAT3 npos, const FLOAT node_mass, FLOAT eps)
{
	FLOAT3 acc;
	acc.x = 0.0f;
	acc.y = 0.0f;
	acc.z = 0.0f;

//	FLOAT d[3];
	FLOAT3 d;
	FLOAT fac;
//	dvec(ppos, npos, d);
	d = npos - ppos;

#ifdef SOFTENING_PLUMMER
//	fac = node_mass / pow(d.x*d.x + d.y*d.y + d.z*d.z + eps*eps, (FLOAT)(3.0f/2.0f)); //thats the Plummer sphere softening
	FLOAT dist = d.x*d.x + d.y*d.y + d.z*d.z;

/*	FLOAT disteps = dist + eps*eps;
	if(dist > eps*eps)
		fac = node_mass / dist * rsqrt(dist);
	else
		fac = node_mass / disteps * rsqrt(disteps);
*/

        dist = sqrt(dist + eps*eps);
	fac = node_mass / (dist*dist*dist);
//	fac = node_mass / dist * rsqrt(dist);
#elif defined(SOFTENING_SPLINEKERNEL)
	//NOTE: gadget/arepo: 	1) far away: use just newtonian point mass potential
	//					2) inside softening: spline kernel potential
	//					3) in tree-node: max softening of all particles in the node
	//					4)!!!! open cell if spline kernels intersect!!!
//	FLOAT dist = (d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
	FLOAT dist = d.x*d.x + d.y*d.y + d.z*d.z;

	if(dist == 0.0f) return acc;

//	dist = sqrt(dist);
	
	if(dist > eps)
	{
		fac = node_mass / dist * rsqrt(dist);
	}
//		fac = node_mass / (dist*dist*dist);
	else
	{
		FLOAT eps_inv = 1.0f/eps;
		FLOAT u = dist*eps_inv;
		FLOAT eps3_inv = eps_inv*eps_inv*eps_inv;

		FLOAT u2 = u * u;

		if(u < 0.5f)
			fac = node_mass * eps3_inv * (10.666666666667f + u2 * (32.0f * u - 38.4f));
		else //u >= 0.5
			fac = node_mass * eps3_inv * (21.333333333333f - 48.0f * u + 38.4f * u2 - 10.666666666667f * u2 * u - 0.066666666667f / (u2 * u));
	}
/*		FLOAT u = dist/eps;
		FLOAT eps_inv = 1.0f/eps;
		FLOAT eps3_inv = eps_inv*eps_inv*eps_inv;

		FLOAT u2 = u*u;
		FLOAT u3 = u2 * u;
		FLOAT nme3i = node_mass * eps3_inv;

		FLOAT fac2 = nme3i * (10.666666666667f + u2 * (32.0f * u - 38.4f));
		FLOAT fac3 = nme3i * (21.333333333333f - 48.0f * u + 38.4f * u2 - 10.666666666667f * u3 - 0.066666666667f / (u3));

		FLOAT flag2 = (float)(u < 0.5f);

		if(u < 0.5f) 
			fac = fac2;
		else
			fac = fac3;
	}*/
#else
	printf("choose a softening type!\n");
	return acc;
#endif

//	acc.x = d[0] * fac;
//	acc.y = d[1] * fac;
//	acc.z = d[2] * fac;
	acc = d * fac;

//	if(acc.x != 0 && acc.y != 0 && acc.z != 0)
//		printf("acc1 %f %f %f\n", acc.x, acc.y, acc.z);

	return acc;
}

__device__ bool openCell(const FLOAT3 ppos, const struct KdNode n, const FLOAT totAcc, FLOAT eps)
{
	const FLOAT alpha = 0.001f; //this is ErrTolForceAcc
	//FLOAT d = dist(n.center_geometric, ppos);
//	FLOAT d[3];
	FLOAT3 d;
	FLOAT r2;		//distance to center of mass of the node!!!!!
//	dvec(ppos, n.center_of_mass, d);
	d = n.center_of_mass - ppos;

//	r2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
	r2 = d.x*d.x + d.y*d.y + d.z*d.z;

	//cell can be used if these three conditions are fulfilled

	if(n.mass * n.l * n.l > r2 * r2 * alpha*totAcc)
		return true;
#ifdef VECTOR
        int3 flags = (fabs(ppos.x - n.center_geometric.x) < 0.6*n.l);
        if(all(flags))
#else
	else if((fabs(ppos.x - n.center_geometric.x) < 0.6*n.l) &&
		(fabs(ppos.y - n.center_geometric.y) < 0.6*n.l) &&
		(fabs(ppos.z - n.center_geometric.z) < 0.6*n.l))
#endif
		return true; 
	else if(sqrt(r2) < 2*eps)
		return true;
	return false;
/*
	bool c1 = (n.mass * n.l * n.l > r2 * r2 * alpha*totAcc);
//if(c1) return true;
	int3 flags = (fabs(ppos.x - n.center_geometric.x) < 0.6*n.l);
	bool c2 = //all(flags);
			((fabs(ppos.x - n.center_geometric.x) < 0.6f*n.l) &&
		(fabs(ppos.y - n.center_geometric.y) < 0.6f*n.l) &&
		(fabs(ppos.z - n.center_geometric.z) < 0.6f*n.l));
//if(c2) return true;
	bool c3 = r2 < (4*eps*eps);
//return c3;

	return (c1 || c2 || c3);*/
}

//#define CACHING

/**
 * performs the tree walk to calculate the acceleration for each particle
 * @param kdtree an array of nodes holding all nodes of the kdtree in depth first ordering
 * @param particles an array with all particles
 * @param nParticles the total number of particles
 * @param eps the softening length
 */
__kernel void treeWalk(__global struct KdNode* kdtree, __global struct Particle* particles, UINT nParticles, FLOAT eps) {
	UINT p = get_global_id(0);

	if(p >= nParticles)
		return;

	struct Particle particle = particles[p];
#ifdef CACHING
	__local struct KdNode nodeCache[256];
	uint id = get_local_id(0);
	nodeCache[id] = kdtree[id];
//	nodeCache[id + get_local_size(0)] = kdtree[id + get_local_size(0)];
	barrier(CLK_LOCAL_MEM_FENCE);
#endif
//uint cnt = 0;
//uint offset = 0;
	UINT current = 0;
	const UINT end = kdtree->size;
	FLOAT3 acc;
	acc.x = 0.0f;
	acc.y = 0.0f;
	acc.z = 0.0f;

	while(current < end)
	{
#ifdef CACHING
		struct KdNode node = (current < 256) ? nodeCache[current] : kdtree[current];
#else
		struct KdNode node = kdtree[current];
#endif

//		bool stop = node.leaf != 0 || !openCell(particle.pos, node, particle.totAcc, eps);
//		acc += stop ? calcAcc(particle.pos, node.center_of_mass, node.mass, eps) : 0;
//		current = stop ? node.size + current : node.left_child;
#if 1
		//stop decending and calculate force on leaves or nodes which do not fulfill the cell opening criterion
		if(node.leaf != 0 || !openCell(particle.pos, node, particle.totAcc, eps))
		{
//printf("leaf    acc: %f %f %f\n", node.cmArr[0], node.cmArr[1], node.cmArr[2]);
			acc += calcAcc(particle.pos, node.center_of_mass, node.mass, eps);
			current += node.size;
		}
		//open the cell and walk further
		else 
		{
//printf("opended acc: %f %f %f\n", node.cmArr[0], node.cmArr[1], node.cmArr[2]);
//printf("%d left child: %d\n", current, node.left_child);
			current = node.left_child;
		}
//++cnt;
//if(cnt != current) {
//printf("c %d %d\n", cnt, current);
//offset = current - cnt;
//}
#endif
	}

	particles[p].totAcc = sqrt(pow(acc.x,2) + pow(acc.y,2) + pow(acc.z,2));
	particles[p].acc.x = acc.x * G;
	particles[p].acc.y = acc.y * G;
	particles[p].acc.z = acc.z * G;
}
