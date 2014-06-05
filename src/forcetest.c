#include <math.h>
#include <stdio.h>

#include "forcetest_c.h"

void dvec(const FLOAT a[3], const FLOAT b[3], FLOAT res[3])
{
	res[0] = b[0] - a[0];
	res[1] = b[1] - a[1];
	res[2] = b[2] - a[2];
}

FLOAT3 dvecVec(const FLOAT3 a, const FLOAT b[3])
{
	FLOAT3 res;
	res.x = b[0] - a.x;
	res.y = b[1] - a.y;
	res.z = b[2] - a.z;
	return res;
}

FLOAT distance(const FLOAT a[3], const FLOAT b[3])
{
	FLOAT d[3]; //= {b[0]-a[0], b[1]-a[1], b[2]-a[2]};
	dvec(a, b, d);
	return sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
}

FLOAT distanceVec(const FLOAT3 a, const FLOAT b[3])
{
	FLOAT3 d; //= {b[0]-a[0], b[1]-a[1], b[2]-a[2]};
	d = dvecVec(a, b);
	return sqrt(d.x*d.x + d.y*d.y + d.z*d.z);
}

FLOAT distanceVecVec(const FLOAT3 a, const FLOAT3 b)
{
	FLOAT3 d; //= {b[0]-a[0], b[1]-a[1], b[2]-a[2]};
	d.x = b.x - a.x;
	d.y = b.y - a.y;
	d.z = b.z - a.z;
	return sqrt(d.x*d.x + d.y*d.y + d.z*d.z);
}

void check_force(const char* fn_in, const char* fn_out, struct Particle* particles, UINT* particleIds, UINT nParticles)
{
	//error estimator
	FLOAT var_gadget_tree = 0.0;
	FLOAT var_kd_tree = 0.0;
	
	//read gadget force file
	FILE* in = fopen(fn_in, "r");
	FILE* out = fopen(fn_out, "w");
	
//	fprintf(out, "# (Gadget tree <-> Gadget direct) (kd tree <-> Gadget tree) (kd tree <-> Gadget direct)\n");
FLOAT max_ = 0.0f;
	for(unsigned j = 0; j < nParticles; ++j)
	{
		struct forcetest_data ftd;
		
		fscanf(in, "%u", &ftd.type);
		fscanf(in, "%u", &ftd.id);
		fscanf(in, "%f", &ftd.time);
		fscanf(in, "%f", &ftd.time_tree);
		fscanf(in, "%f", &ftd.pos[0]);
		fscanf(in, "%f", &ftd.pos[1]);
		fscanf(in, "%f", &ftd.pos[2]);
		fscanf(in, "%f", &ftd.accDirect[0]);
		fscanf(in, "%f", &ftd.accDirect[1]);
		fscanf(in, "%f", &ftd.accDirect[2]);
		fscanf(in, "%f", &ftd.acc[0]);
		fscanf(in, "%f", &ftd.acc[1]);
		fscanf(in, "%f", &ftd.acc[2]);
			
		//norm(ftd.acc);
		//norm(ftd.accDirect);

		//find corresponding particle
		for(UINT i = 0; i < nParticles; ++i)
		{
			if(ftd.id == particleIds[i])
			{
				FLOAT ref =sqrt(ftd.accDirect[0]*ftd.accDirect[0] + ftd.accDirect[1]*ftd.accDirect[1] + ftd.accDirect[2]*ftd.accDirect[2]);
				
				FLOAT tg_d = distance(ftd.accDirect, ftd.acc) / ref;
				FLOAT tkd_tg = distanceVec(particles[i].acc, ftd.acc) / ref;
				FLOAT tkd_d = distanceVec(particles[i].acc, ftd.accDirect) / ref;

//				fprintf(out, "%u %f %f %f :::: %f %f %f\n", ftd.id, tg_d, tkd_tg, tkd_d, particles[i].acc.x, particles[i].acc.y, particles[i].acc.z);
				fprintf(out, "%u %f %f %f \n", ftd.id,  particles[i].acc.x, particles[i].acc.y, particles[i].acc.z);


				var_gadget_tree += tg_d*tg_d;
				var_kd_tree += tkd_d*tkd_d;

				max_ = max_ < (tkd_d) ? tkd_d : max_;
				
				break;
			}
		}
	}

	fflush(out);	
	printf("mean sq. error: GADGET: %f KDTREE: %f %f\n", var_gadget_tree/nParticles, var_kd_tree/nParticles, max_);
	
}

void check_force_internal(const char* fn_out, struct Particle* refParticles, struct Particle* particles, UINT* particleIds, UINT nParticles)
{
	//error estimator
	FLOAT var_kd_tree = 0.0;

	//read gadget force file
	FILE* out = fopen(fn_out, "w");

	FLOAT max_ = 0.0f;
	for(unsigned i = 0; i < nParticles; ++i)
	{
		FLOAT ref =sqrt(refParticles[i].acc.x*refParticles[i].acc.x +
				refParticles[i].acc.y*refParticles[i].acc.y + refParticles[i].acc.z*refParticles[i].acc.z);

		FLOAT tkd_d = distanceVecVec(particles[i].acc, refParticles[i].acc) / ref;

		fprintf(out, "%u %f %f %f \n", particleIds[i],  particles[i].acc.x, particles[i].acc.y, particles[i].acc.z);

		var_kd_tree += tkd_d*tkd_d;

		max_ = max_ < (tkd_d) ? tkd_d : max_;
	}

	fflush(out);
	printf("KDTREE: mean sq. error: %f max: %f\n", var_kd_tree/nParticles, max_);

}

