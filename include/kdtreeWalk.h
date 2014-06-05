//kdtreeWalk.h

#ifndef KDTREEWALK_H
#define KDTREEWALK_H

#include "kdtree.h"

class TreeWalker
{
private:
	FLOAT eps; 		//softening length
	
	bool openCell(const FLOAT *const ppos, const Node *const n, const FLOAT totAcc);
	void walkParticle(const Particle *const p, const Node *const n, FLOAT acc[3]);
	void walkParticleIt(const Particle *const p, const Node *const kdTree, const unsigned totParticles, FLOAT acc[3]);
	void calcAcc(const FLOAT ppos[3], const FLOAT npos[3], const FLOAT node_mass, FLOAT acc[3]);
	
	
public:
	TreeWalker(const FLOAT eps);
	void walk(Tree *tree, const Node *const kdTree, const unsigned totParticles);
	
	static FLOAT distance(const FLOAT a[3], const FLOAT b[3]);
	static void dvec(const FLOAT a[3], const FLOAT b[3], FLOAT res[3]);
};

#endif
