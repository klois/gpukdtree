#pragma once


#include "gpukdtree.h"
#include "gpuKdtreeWalk.h"

struct forcetest_data
{
	UINT type, id;
	FLOAT time, time_tree;
	FLOAT pos[3];
	FLOAT accDirect[3];
	FLOAT acc[3];
	
	/*
	fprintf(FdForceTest, "%d %d %g %g %g %g %g %g %g %g %g %g %g\n",
	P[i].Type, P[i].ID, All.Time, All.Time - TimeOfLastTreeConstruction,
	P[i].Pos[0], P[i].Pos[1], P[i].Pos[2],
	P[i].GravAccelDirect[0], P[i].GravAccelDirect[1], P[i].GravAccelDirect[2],
	P[i].GravAccel[0], P[i].GravAccel[1], P[i].GravAccel[2]);
	*/
};


void extern check_force(const char* fn_in, const char* fn_out, struct Particle* particles, UINT* particleIds, UINT nParticles);
void check_force_internal(const char* fn_out, struct Particle* refParticles, struct Particle* particles, UINT* particleIds, UINT nParticles);
