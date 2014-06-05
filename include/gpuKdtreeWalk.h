//gpuKdtreeWalk.h

#ifndef GPUKDTREEWALK_H
#define GPUKDTREEWALK_H

#include "gpukdtree.h"

void walk(icl_buffer* kdTree, icl_buffer* particles, const UINT nParticles, FLOAT eps, icl_device* dev);

#endif
