//gpuKdtreeWalk.h

//#ifndef GPUKDTREE_H
//#define GPUKDTREE_H
#pragma once

#include <float.h>
#include <stdio.h>

#ifdef CUDA
#include "lib_icl_cuda.h"
#else
#include "lib_icl.h"
#endif
#include "lib_icl_ext.h"
#include "timer.h"

#ifdef CUDA
#define FLOAT 		float
#define FLOAT3 		float4
#define INT 		int
#define UINT 		unsigned int
#else
#define FLOAT 		cl_float
#define FLOAT3 		cl_float3
#define INT 		cl_int
#define UINT 		cl_uint
#define x s[0]
#define y s[1]
#define z s[2]
#endif
#define FMAX		FLT_MAX

#ifndef max 
#define max(a, b) (a > b ? a : b)
#endif
#ifndef min
#define min(a, b) (a < b ? a : b)
#endif

#include "structs.h"

FLOAT getBoxVolume(struct BBox bbox);

void getBoxCenter(FLOAT3 center, struct BBox bbox);

void initBBox(struct BBox* bbox);

void printBox(struct BBox box);
//#endif

#define timing 1
