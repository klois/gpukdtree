#include "device.h"

// I'm sick of this
#ifndef __device__
#define __device__
#endif

#include "structs.h"

#define SOFTENING_SPLINEKERNEL

//#define ErrTolIntAccuracy 0.025 //put this in parameterfile or whatever


/**
 * calculates the timestep for each particle
 * @param particles an array with all particles
 * @param nParticles the total number of particles
 * @param eps the softening length
 * @param ErrTolIntAccuracy accuracy switch for time integration
 */


__kernel void drift_particles(__global struct Particle* particles, UINT nParticles, FLOAT dt)
{
	UINT p = get_global_id(0);

	if(p >= nParticles)
		return;

	struct Particle particle = particles[p];

/*	FLOAT3 tmp = particle.pos;
	int3 tmpI;
	tmpI.x = (int)particle.pos.x;
        tmpI.y = (int)particle.pos.y;
        tmpI.z = (int)particle.pos.z;

	tmp.x = tmp.x < 0 ?  tmp.x - tmpI.x : tmp.x + tmpI.x;
        tmp.y = tmp.y < 0 ?  tmp.y - tmpI.y : tmp.y + tmpI.y;
        tmp.z = tmp.z < 0 ?  tmp.z - tmpI.z : tmp.x + tmpI.z;

	tmp += dt *  particle.vel;
	particle.pos.x = tmpI.x < 0 ? tmp.x + tmpI.x : tmp.x - tmpI.x;
        particle.pos.y = tmpI.y < 0 ? tmp.y + tmpI.y : tmp.y - tmpI.y;
        particle.pos.z = tmpI.z < 0 ? tmp.z + tmpI.z : tmp.z - tmpI.z;
*/
	
	particle.pos.x += dt * particle.vel.x;
	particle.pos.y += dt * particle.vel.y;
	particle.pos.z += dt * particle.vel.z;

	particles[p] = particle;

}


__kernel void kick_particles(__global struct Particle* particles, UINT nParticles, FLOAT dt)
{
	UINT p = get_global_id(0);

	if(p >= nParticles)
		return;

	struct Particle particle = particles[p];
	
	particle.vel.x += dt * particle.acc.x;
	particle.vel.y += dt * particle.acc.y;
	particle.vel.z += dt * particle.acc.z;

	particles[p] = particle;

//if(p == 0)
//        printf("Kick Vel: %f\n",  particles[0].vel.x);
}

