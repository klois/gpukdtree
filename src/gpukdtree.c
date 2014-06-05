#include <math.h>
#include <assert.h>

#include "gpukdtree.h"

#include "gpuKdtreeWalk.h"
#include "forcetest_c.h"

#include "scan.h" // efficient scan & prefix scan implementation

#include "snapshot_io.h"
#include "data_ops.h"

//#define KERNEL_BUILD_MACRO "-cl-mad-enable -Iinclude -cl-unsafe-math-optimizations -cl-fast-relaxed-math"
#define KERNEL_BUILD_MACRO "-Iinclude "
#define DEVICE ICL_GPU

#ifdef _WIN32
#define fmin min
#endif

void run(const FLOAT t_max, const FLOAT eps, const FLOAT ErrTolIntAccuracy, struct Particle *particles_host, icl_buffer* particles_device, const UINT nParticles, icl_device* dev, struct Tree *tree, struct Particle *ref);


FLOAT getBoxVolume(struct BBox bbox)
{
	return (bbox.box[1].x - bbox.box[0].x) * (bbox.box[1].y - bbox.box[0].y) * (bbox.box[1].z - bbox.box[0].z);
}

void getBoxCenter(FLOAT3 center, struct BBox bbox)
{
	center.x = (bbox.box[1].x + bbox.box[0].x)/2.0;
	center.y = (bbox.box[1].y + bbox.box[0].y)/2.0;
	center.z = (bbox.box[1].z + bbox.box[0].z)/2.0;
}

void initBBox(struct BBox* bbox)
{
	bbox->box[0].x = bbox->box[0].y = bbox->box[0].z = FMAX;
	bbox->box[1].x = bbox->box[1].y = bbox->box[1].z = FMAX * -1.0f;
};

void printBox(struct BBox box) 
{
	printf("bbox: %0.2f %0.2f : %0.2f %0.2f : %0.2f %0.2f\n",
			box.box[0].x, box.box[1].x, box.box[0].y, box.box[1].y, box.box[0].z, box.box[1].z);
}
/// Round up to next higher power of 2 (return x if it's already a power of 2).
// http://stackoverflow.com/questions/364985/algorithm-for-finding-the-smallest-power-of-two-thats-greater-or-equal-to-a-giv
UINT pow2roundup (UINT a)
{
//    if (x == 0)
//        return 2;
	--a;
	a |= a >> 1;
	a |= a >> 2;
	a |= a >> 4;
	a |= a >> 8;
	a |= a >> 16;
	return a+1;
}

// swap buffers 
void swap(icl_buffer** a, icl_buffer** b)
{
	icl_buffer* tmp = *a;
	*a = *b;
	*b = tmp;
}

UINT buildTree(icl_buffer *nodelist, icl_buffer *particlesD, icl_buffer *treeD, UINT nParticles, icl_device* dev)
{
	UINT level = 1;
	UINT nNodes = nParticles * 2 - 1;

	icl_timer* timer = icl_init_timer(ICL_MILLI);
//	void icl_start_timer(icl_timer* timer);
	double time = 0;

	// overapproximate size of temporal lists
/*	struct Node** activelist = (struct Node**)malloc(nParticles * sizeof(struct Node*)); UINT activeN = 0;
	struct Node** smalllist = (struct Node**)malloc(nParticles * sizeof(struct Node*)); UINT smallN = 0;
	struct Node** nextlist = (struct Node**)malloc(nParticles * sizeof(struct Node*)); UINT nextN = 0;*/
	icl_buffer* activelist = icl_create_buffer(dev, CL_MEM_READ_WRITE, sizeof(NodeId) * nParticles);
	icl_buffer* smalllist = icl_create_buffer(dev, CL_MEM_READ_WRITE, sizeof(NodeId) * nParticles);
	icl_buffer* nextlist = icl_create_buffer(dev, CL_MEM_READ_WRITE, sizeof(NodeId) * nParticles);
	icl_buffer* sizes = icl_create_buffer(dev, CL_MEM_READ_WRITE, sizeof(UINT) * 5); // holds the current size of each of 3 buffers:
		// 0 activelist
		// 1 nodelist
		// 2 smalllist
		// 3 old max level
		// 4 new max level

	UINT maxNchunks = ((nParticles / fmin(T, chunk_size)) * 2) -1;
//	assert(maxNchunks <= 256 && "adapt implementation"); // TODO allow more than 256 chunks per node
	icl_buffer* chunks = icl_create_buffer(dev, CL_MEM_READ_WRITE, sizeof(struct Chunk) * maxNchunks);
	icl_buffer* bboxes = icl_create_buffer(dev, CL_MEM_READ_WRITE, sizeof(struct BBox) * maxNchunks);

	size_t localSize = 1;
	size_t globalSize = 1;

/*
struct Particle* particles = (struct Particle*)malloc(3000 * sizeof(struct Particle));
icl_read_buffer(particlesD, CL_TRUE, sizeof(struct Particle) * 3000, &particles[0], NULL, NULL);
printf("%f %f %f\n", particles[0].pos.x, particles[0].pos.y, particles[0].pos.z);
*/
	// compile OpenCL kernels
	icl_kernel* init = icl_create_kernel(dev, "kernel/init.cl", "init", KERNEL_BUILD_MACRO, ICL_SOURCE);
	icl_kernel* resetChunks = icl_create_kernel(dev, "kernel/init.cl", "memset_chunks", KERNEL_BUILD_MACRO, ICL_SOURCE);
	icl_kernel* gp2c = icl_create_kernel(dev, "kernel/groupToChunks.cl", "groupParticlesIntoChunks", KERNEL_BUILD_MACRO, ICL_SOURCE);
	icl_kernel* cBBox = icl_create_kernel(dev, "kernel/chunkedBBox.cl", "chunkedBBox", KERNEL_BUILD_MACRO, ICL_SOURCE);
	icl_kernel* bBox = icl_create_kernel(dev, "kernel/bBox.cl", "bBox", KERNEL_BUILD_MACRO, ICL_SOURCE);
	icl_kernel* sln = icl_create_kernel(dev, "kernel/splitLargeNodes.cl", "splitLargeNodes", KERNEL_BUILD_MACRO, ICL_SOURCE);
	icl_kernel* sortP = icl_create_kernel(dev, "kernel/sortParticlesToChilds.cl", "sortParticlesToChilds", KERNEL_BUILD_MACRO, ICL_SOURCE);
	icl_kernel* snf = icl_create_kernel(dev, "kernel/smallNodeFiltering.cl", "smallNodeFiltering", KERNEL_BUILD_MACRO, ICL_SOURCE);
	icl_kernel* pnl = icl_create_kernel(dev, "kernel/packNextlist.cl", "packNextlist", KERNEL_BUILD_MACRO, ICL_SOURCE);	

	//////////////////////////////////////////////////////////////////////////
	icl_kernel* preScan  = icl_create_kernel(dev, "kernel/sortP_prescan.cl", "sortP_prescan_chunked", KERNEL_BUILD_MACRO, ICL_SOURCE);
	icl_kernel* postScan = icl_create_kernel(dev, "kernel/sortP_postscan.cl", "sortP_postscan_chunked", KERNEL_BUILD_MACRO, ICL_SOURCE);
	segmented_scan_init(nParticles, dev, KERNEL_BUILD_MACRO, ICL_SOURCE);
	icl_kernel* memset_int_s  = icl_create_kernel(dev, "kernel/init.cl", "memset_int_s", KERNEL_BUILD_MACRO, ICL_SOURCE);

	// approach with segmented scan
	icl_buffer *scan_data = icl_create_buffer(dev, CL_MEM_READ_WRITE, sizeof(int) * nParticles);
	icl_buffer *scan_flag = icl_create_buffer(dev, CL_MEM_READ_WRITE, sizeof(int) * nParticles);
	icl_buffer* buffered_particles = icl_create_buffer(dev, CL_MEM_READ_WRITE, sizeof(struct Particle) * nParticles);


#if timing == 1
	icl_timer* timer_gp2c =  icl_init_timer(ICL_MILLI);
	icl_timer* timer_cBBox =  icl_init_timer(ICL_MILLI);
	icl_timer* timer_bBox =  icl_init_timer(ICL_MILLI);
	icl_timer* timer_sln =  icl_init_timer(ICL_MILLI);
	icl_timer* timer_sortP =  icl_init_timer(ICL_MILLI);
	icl_timer* timer_snf =  icl_init_timer(ICL_MILLI);
	icl_timer* timer_pnl =  icl_init_timer(ICL_MILLI);
	icl_timer* timer_ran =  icl_init_timer(ICL_MILLI);


	icl_timer* timer_prescan =  icl_init_timer(ICL_MILLI);
	icl_timer* timer_scan =  icl_init_timer(ICL_MILLI);
	icl_timer* timer_postscan =  icl_init_timer(ICL_MILLI);

#endif
	//add root node to the activelist and initialize size lists
	icl_run_kernel(init, 1, &globalSize, &localSize, NULL, NULL, 3,
			(size_t)0, (void *)activelist,
			(size_t)0, (void *)sizes,
			(size_t)0, (void *)particlesD);

	UINT activeN = 1;

	icl_finish(dev);

	// smallest power of 2 bigger or equal to maxxNchnunks
	UINT pow2maxNchunks = pow2roundup(maxNchunks);
	// processLargeNode
	while(activeN != 0)
	{

		icl_start_timer(timer);
		// group triangles into chunks
		size_t localSize1 = min(pow2maxNchunks, 256);

#if timing == 1
		icl_start_timer(timer_gp2c);
#endif
		size_t globalSize1 = ((maxNchunks + localSize1 -1) / localSize1) * localSize1;
		// reset chunks
		icl_run_kernel(resetChunks, 1, &globalSize1, &localSize1, NULL, NULL, 2,
				(size_t)0, (void *)chunks,
				sizeof(UINT), &maxNchunks);

		globalSize1 = localSize1 * activeN;
		// split every node in chunk of chunk_size
		icl_run_kernel(gp2c, 1, &globalSize1, &localSize1, NULL, NULL, 4,
				(size_t)0, (void *)nodelist,
				(size_t)0, (void *)activelist,
				(size_t)0, (void *)chunks,
				sizeof(UINT), &activeN);
#if timing == 1
		icl_finish(dev);
		icl_stop_timer(timer_gp2c);
#endif

		// compute per chunk bounding box
		size_t localSize2 = chunk_size;
		size_t globalSize2 = maxNchunks * chunk_size;
#if timing == 1
		icl_start_timer(timer_cBBox);
#endif

		icl_run_kernel(cBBox, 1, &globalSize2, &localSize2, NULL, NULL, 5,
				(size_t)0, (void *)nodelist,
				(size_t)0, (void *)activelist,
				(size_t)0, (void *)particlesD,
				(size_t)0, (void *)chunks,
				(size_t)0, (void *)bboxes);
#if timing == 1
		icl_finish(dev);
		icl_stop_timer(timer_cBBox);
#endif

		// compute per node bounding box
		size_t localSize3 = min(pow2maxNchunks, 256);
		size_t globalSize3 = localSize3 * activeN;
#if timing == 1
		icl_start_timer(timer_bBox);
#endif
		icl_run_kernel(bBox, 1, &globalSize3, &localSize3, NULL, NULL, 4,
						(size_t)0, (void *)nodelist,
						(size_t)0, (void *)activelist,
						(size_t)0, (void *)bboxes,
						sizeof(UINT), &activeN);
#if timing == 1
		icl_finish(dev);
		icl_stop_timer(timer_bBox);
#endif

		// split large nodes
		size_t localSize4 = 256;
		size_t globalSize4 = ((activeN + 255) / 256) * 256;
#if timing == 1
		icl_start_timer(timer_sln);
#endif
		icl_run_kernel(sln, 1, &globalSize4, &localSize4, NULL, NULL, 5,
						(size_t)0, (void *)nodelist,
						(size_t)0, (void *)activelist,
						(size_t)0, (void *)nextlist,
						(size_t)0, (void *)sizes,
						sizeof(UINT), &activeN);
#if timing == 1
		icl_finish(dev);
		icl_stop_timer(timer_sln);
#endif

		///////////////////////////////////////////////////////////////////////////////
		// XXx replaced with segmented scan

		//globalSize = (activeN+1) * 256;
#if timing == 1
		icl_start_timer(timer_sortP);
#endif


#if DEVICE == ICL_CPU 
		// sort particles to child nodes
		size_t localSize5 = 256;
		size_t globalSize5 = ((activeN + 255) / 256) * 256;
		icl_run_kernel(sortP, 1, &globalSize5, &localSize5, NULL, NULL, 5,
			(size_t)0, (void *)nodelist,
			(size_t)0, (void *)particlesD,
			(size_t)0, (void *)activelist,
			(size_t)0, (void *)nextlist,
			sizeof(UINT), &activeN
		);
#else		
		// init scan_flag to 1
		cl_int initFlag = 1;
		size_t np = (size_t)((nParticles + localSize4 -1 ) / localSize4) * localSize4;
		icl_run_kernel(memset_int_s, 1, &np, &localSize4, NULL, NULL, 3,
			(size_t)0, (void *)scan_flag,
			sizeof(cl_int), &initFlag,
			sizeof(UINT), &nParticles
		);

#if timing == 1
		icl_start_timer(timer_prescan);
#endif
		// pre-scan fills data0 and data1 with 1 and 0 whenever value < pivot
		localSize = chunk_size;
//		globalSize = activeN * 256;
		globalSize = maxNchunks * chunk_size;
		icl_run_kernel(preScan, 1, &globalSize, &localSize, NULL, NULL, 7,
			(size_t)0, (void *)nodelist,
			(size_t)0, (void *)chunks,
			(size_t)0, (void *)particlesD,
			(size_t)0, (void *)activelist,
			sizeof(UINT), &activeN,
			(size_t)0, (void *)scan_data,
			(size_t)0, (void *)scan_flag						
		);
#if timing == 1
		icl_finish(dev);
		icl_stop_timer(timer_prescan);
#endif

#if timing == 1
		icl_start_timer(timer_scan);
#endif
		// scan for
		scan(scan_data, scan_flag, nParticles);
#if timing == 1
		icl_finish(dev);
		icl_stop_timer(timer_scan);
#endif
		// copy partially sorted data to the final
		icl_copy_buffer(particlesD, buffered_particles, sizeof(struct Particle) * nParticles, NULL, NULL);
//		swap(particlesD, buffered_particles);
#if timing == 1
		icl_start_timer(timer_postscan);
#endif

		localSize = chunk_size;
		globalSize = maxNchunks * chunk_size;
		icl_run_kernel(postScan, 1, &globalSize, &localSize, NULL, NULL, 7,
			(size_t)0, (void *)nodelist,
			(size_t)0, (void *)chunks,
			(size_t)0, (void *)particlesD,
			(size_t)0, (void *)activelist,
			sizeof(UINT), &activeN,

			(size_t)0, (void *)buffered_particles,
			(size_t)0, (void *)scan_data						
		);
#if timing == 1
		icl_finish(dev);
		icl_stop_timer(timer_postscan);
#endif

icl_finish(dev);
#endif
		///////////////////////////////////////////////////////////////////////////////

#if timing == 1
		icl_finish(dev);
		icl_stop_timer(timer_sortP);
#endif
		
		// small node filtering
		size_t localSize6 = 256;
		size_t globalSize6 = ((activeN*2 + 255) / 256) * 256;
#if timing == 1
		icl_start_timer(timer_snf);
#endif
		icl_run_kernel(snf, 1, &globalSize6, &localSize6, NULL, NULL, 4,
			(size_t)0, (void *)nodelist,
			(size_t)0, (void *)nextlist,
			(size_t)0, (void *)smalllist,
			(size_t)0, (void *)sizes
			//, sizeof(UINT), &nParticles
		);
#if timing == 1
		icl_finish(dev);
		icl_stop_timer(timer_snf);
#endif

		// packing of nextlist
		size_t localSize7 = 1;
		size_t globalSize7 = 1;
#if timing == 1
		icl_start_timer(timer_pnl);
#endif
		icl_run_kernel(pnl, 1, &globalSize7, &localSize7, NULL, NULL, 2,
			(size_t)0, (void *)nextlist,
			(size_t)0, (void *)sizes
		);
#if timing == 1
		icl_finish(dev);
		icl_stop_timer(timer_pnl);
#endif

		// swap nextlist and activelist
		swap(&nextlist, &activelist);

#if timing == 1
		icl_start_timer(timer_ran);
#endif
		// read size of next activelist set in kernel
		icl_read_buffer(sizes, CL_TRUE, sizeof(UINT), &activeN, NULL, NULL);
		icl_finish(dev);
#if timing == 1
		icl_stop_timer(timer_ran);
#endif

		++level;
//printf("%d: ActiveN %d\n", level, activeN);
		time = icl_stop_timer(timer);

	}

	icl_release_kernel(init);
	icl_release_kernel(gp2c);
	icl_release_kernel(cBBox);
	icl_release_kernel(bBox);
	icl_release_kernel(sln);
	icl_release_kernel(sortP);
	icl_release_kernel(snf);
	icl_release_kernel(pnl);
//////////////////////////////////////////////////////////////////////////
	icl_release_kernel(preScan);
	icl_release_kernel(postScan);
	segmented_scan_release();

	icl_release_buffers(3, scan_data, scan_flag, buffered_particles);

#if timing == 1
	printf("gp2c %f\ncBBox %f\nbBox %f\nsln  %f\nsortP %f\nsnf %f\npnl %f\nran %f\n\n",
			timer_gp2c->current_time,
			timer_cBBox->current_time,
			timer_bBox->current_time,
			timer_sln->current_time,
			timer_sortP->current_time,
			timer_snf->current_time,
			timer_pnl->current_time,
			timer_ran->current_time);
	icl_release_timer(timer_gp2c);
	icl_release_timer(timer_cBBox);
	icl_release_timer(timer_bBox);
	icl_release_timer(timer_sln);
	icl_release_timer(timer_sortP);
	icl_release_timer(timer_snf);
	icl_release_timer(timer_pnl);

	printf("prescan %f\nscan %f\npostscan %f\n\n", timer_prescan->current_time, timer_scan->current_time, timer_postscan->current_time);

	icl_release_timer(timer_prescan);
	icl_release_timer(timer_scan);
	icl_release_timer(timer_postscan);
#endif

/*
icl_read_buffer(nodelist, CL_TRUE, sizeof(struct Node) * 6000, tree->nodelist, NULL, NULL);
printf("node: %d, left %d, right %d", tree->nodelist[0].particlesHigh - tree->nodelist[0].particlesLow,
		tree->nodelist[1].particlesHigh - tree->nodelist[1].particlesLow, tree->nodelist[2].particlesHigh - tree->nodelist[2].particlesLow);

printBox(tree->nodelist[49].bounding_box);
printBox(tree->nodelist[53].bounding_box);
printBox(tree->nodelist[54].bounding_box);


for(int i = 0; i < 6000; ++i)
	if(tree->nodelist[i].bounding_box.box[0].x != 0.0)
		printBox(tree->nodelist[i].bounding_box);
*/
	//small nodes stage
//	preprocessSmallNodes(smalllist);
	icl_release_buffers(3, activelist, chunks, bboxes);
	icl_kernel* sasl = icl_create_kernel(dev, "kernel/swapActiveAndSmalllist.cl", "swapActiveAndSmalllist", KERNEL_BUILD_MACRO, ICL_SOURCE);
	icl_kernel* ssn = icl_create_kernel(dev, "kernel/splitSmallNodes.cl", "splitSmallNodes", KERNEL_BUILD_MACRO, ICL_SOURCE);

#if timing == 1
	icl_timer* timer_ssn = icl_init_timer(ICL_MILLI);
	icl_timer* timer_sasl = icl_init_timer(ICL_MILLI);
	icl_timer* timer_rsn = icl_init_timer(ICL_MILLI);
#endif

	size_t localSize8 = 1;
	size_t globalSize8 = 1;
	UINT setMaxLevel = 0;
	icl_run_kernel(sasl, 1, &globalSize8, &localSize8, NULL, NULL, 2,
					(size_t)0, (void *)sizes,
					sizeof(UINT), &level);
	// get number of small nodes
	icl_read_buffer(sizes, CL_TRUE, sizeof(UINT), &activeN, NULL, NULL);

	while(activeN != 0)
	{
		icl_start_timer(timer);
		// compute SVH and determine the split plane
		size_t localSize9 = 256;
		size_t globalSize9 = ((activeN + 255) / 256) * 256;
#if timing == 1
		icl_start_timer(timer_ssn);
#endif
		icl_run_kernel(ssn, 1, &globalSize9, &localSize9, NULL, NULL, 5,
						(size_t)0, (void *)nodelist,
						(size_t)0, (void *)smalllist,
						(size_t)0, (void *)nextlist,
						(size_t)0, (void *)particlesD,
						(size_t)0, (void *)sizes);
#if timing == 1
		icl_finish(dev);
		icl_stop_timer(timer_ssn);
#endif

		size_t localSizeA = 1;
		size_t globalSizeA = 1;
#if timing == 1
		icl_start_timer(timer_sasl);
#endif
		icl_run_kernel(sasl, 1, &globalSizeA, &localSizeA, NULL, NULL, 2,
						(size_t)0, (void *)sizes,
						sizeof(UINT), &setMaxLevel);
#if timing == 1
		icl_finish(dev);
		icl_stop_timer(timer_sasl);
#endif

		swap(&nextlist, &smalllist);

		// read size of next activelist set in kernel
#if timing == 1
		icl_start_timer(timer_rsn);
#endif
		icl_read_buffer(sizes, CL_TRUE, sizeof(UINT), &activeN, NULL, NULL);
//printf("small size %d\n", activeN);
		icl_finish(dev);
#if timing == 1
		icl_stop_timer(timer_rsn);
#endif
		time = icl_stop_timer(timer);
	}

	icl_release_buffer(smalllist);
	icl_release_buffer(nextlist);

	icl_release_kernel(sasl);
	icl_release_kernel(ssn);

#if timing == 1
	printf("ssn %f\nsasl %f\nrsn %f\n\n", timer_ssn->current_time, timer_sasl->current_time, timer_rsn->current_time);
	icl_release_timer(timer_ssn);
	icl_release_timer(timer_sasl);
	icl_release_timer(timer_rsn);
#endif

	UINT s[5];
	icl_read_buffer(sizes, CL_TRUE, sizeof(UINT) * 5, &s, NULL, NULL);
	icl_release_buffer(sizes);

	icl_kernel* upPass = icl_create_kernel(dev, "kernel/upPass.cl", "upPass", KERNEL_BUILD_MACRO, ICL_SOURCE);
	icl_kernel* downPass = icl_create_kernel(dev, "kernel/kdDownPass.cl", "downPass", KERNEL_BUILD_MACRO, ICL_SOURCE);

#if timing == 1
	icl_timer* timer_upPass = icl_init_timer(ICL_MILLI);
	icl_timer* timer_downPass = icl_init_timer(ICL_MILLI);
	icl_timer* timer_rt = icl_init_timer(ICL_MILLI);
#endif

	UINT treeHeight = s[4];
	printf("Tree height: %d\n", treeHeight);

	size_t localSizeB = 256;
	size_t globalSizeB = ((nNodes + 255) / 256) * 256;
	icl_start_timer(timer);
#if timing == 1
		icl_start_timer(timer_upPass);
#endif

	for(int l = (int)treeHeight; l >= 0; --l)
	{
		icl_run_kernel(upPass, 1, &globalSizeB, &localSizeB, NULL, NULL, 4,
						(size_t)0, (void *)nodelist,
						(size_t)0, (void *)particlesD,
						sizeof(int), &l,
						sizeof(UINT), &nNodes);
	}
#if timing == 1
		icl_finish(dev);
		icl_stop_timer(timer_upPass);
		icl_start_timer(timer_downPass);
#endif
	for(UINT l = 0; l <= treeHeight; ++l)
	{
		icl_run_kernel(downPass, 1, &globalSizeB, &localSizeB, NULL, NULL, 4,
						(size_t)0, (void *)nodelist,
						(size_t)0, (void *)treeD,
						sizeof(UINT), &l,
						sizeof(UINT), &nNodes);
	}
	icl_finish(dev);
#if timing == 1
	icl_stop_timer(timer_downPass);
#endif
	time = icl_stop_timer(timer);

	icl_release_kernel(upPass);
	icl_release_kernel(downPass);

#if timing == 1
	icl_start_timer(timer_rt);
#endif

#if timing == 1
	icl_finish(dev);
	icl_stop_timer(timer_rt);

	printf("upPass %f\ndownPass %f\nread Tree %f\n\n", timer_upPass->current_time, timer_downPass->current_time, timer_rt->current_time);
	icl_release_timer(timer_upPass);
	icl_release_timer(timer_downPass);
	icl_release_timer(timer_rt);
#endif

	//	struct Node* kdTree = (struct Node*)malloc(sizeof(struct Node) * nNodes);
//	icl_read_buffer(treeD, CL_TRUE, sizeof(struct Node) * nNodes, kdTree, NULL, NULL);
//	printf("%d", tree->nodelist[0].left_child);

	printf("\nTime: %f\n", time);
	icl_release_timer(timer);

	return treeHeight;
}


void updateTree(icl_buffer *nodelist, icl_buffer *particlesD, icl_buffer *treeD, UINT nParticles, UINT treeHeight, icl_device* dev) {
	// compile OpenCL kernels
	icl_kernel* update = icl_create_kernel(dev, "kernel/updateTree.cl", "updateTree", KERNEL_BUILD_MACRO, ICL_SOURCE);
	UINT nNodes = nParticles * 2 - 1;

#if timing == 1
	icl_timer* timer_update = icl_init_timer(ICL_MILLI);
#endif
	
	size_t localSizeB = 256;
	size_t globalSizeB = ((nNodes + 255) / 256) * 256;

#if timing == 1
		icl_start_timer(timer_update);
#endif

	for(int l = (int)treeHeight; l >= 0; --l)
	{
		icl_run_kernel(update, 1, &globalSizeB, &localSizeB, NULL, NULL, 5,
						(size_t)0, (void *)nodelist,
						(size_t)0, (void *)treeD,
						(size_t)0, (void *)particlesD,
						sizeof(int), &l,
						sizeof(UINT), &nNodes);
	}
#if timing == 1
	icl_finish(dev);
	icl_stop_timer(timer_update);

	printf("tree update %f\n\n", timer_update->current_time);
#endif

}


void initBBox2(struct BBox* bbox)
{
	bbox->box[0].x = bbox->box[0].y = bbox->box[0].z = FMAX;
	bbox->box[1].x = bbox->box[1].y = bbox->box[1].z = 0.0;
};

int main (int argc, char **argv)
{
	struct BBox box;
	struct Particle *particles;
	particle_data *P;
	io_header header;
	int tot = snapshotLoader(argv[1], &header, &P);
	int k = 0;

	if(tot <= 0)
	{
		printf("error while loading snapshot file\n");
		return -1;
	}

	initBBox2(&box);

	particles = (struct Particle*)malloc(header.npartTotal[1] * sizeof(struct Particle));
	UINT* particleIds = (UINT*)malloc(header.npartTotal[1] * sizeof(UINT));
	//for(int j = header.npartTotal[0]; j < header.npartTotal[0]+header.npartTotal[1]; ++j)
#define F 1
	for(int j = header.npartTotal[0]; j < header.npartTotal[0]+header.npartTotal[1]; j += F)
	{
		particles[k].pos.x = P[j].Pos[0];
		particles[k].pos.y = P[j].Pos[1];
		particles[k].pos.z = P[j].Pos[2];
		particles[k].vel.x = P[j].Vel[0];
		particles[k].vel.y = P[j].Vel[1];
		particles[k].vel.z = P[j].Vel[2];
		particles[k].mass = P[j].Mass;
		particles[k].id = P[j].Id;
		particleIds[k] = P[j].Id;

		//printf("%f %f %f\n",  particles[k].pos.x, particles[k].pos.y, particles[k].pos.z);

		//get bbox
/*		
		if(particles[k].pos.x < box.box[0].x)
			box.box[0].x = particles[k].pos.x;

		if(particles[k].pos.y < box.box[0].y)
			box.box[0].y = particles[k].pos.y;

		if(particles[k].pos.z < box.box[0].z)
			box.box[0].z = particles[k].pos.z;

		if(particles[k].pos.x >= box.box[1].x)
			box.box[1].x = particles[k].pos.x;

		if(particles[k].pos.y >= box.box[1].y)
			box.box[1].y = particles[k].pos.y;

		if(particles[k].pos.z >= box.box[1].z)
			box.box[1].z = particles[k].pos.z;
*/
		++k;
	}

	free(P);
	
	header.npartTotal[1] /= F;
	struct Tree tree;
	tree.nodelist = (struct Node*)malloc(2*header.npartTotal[1] * sizeof(struct Node));
	struct Particle* ref = (struct Particle*)malloc(sizeof(struct Particle) * header.npartTotal[1]);

	// init ocl
	icl_init_devices(DEVICE);

	if (icl_get_num_devices() != 0)
	{
		icl_device* dev = icl_get_device(0);

		icl_print_device_short_info(dev);

		icl_buffer* particlesD = icl_create_buffer(dev, CL_MEM_READ_WRITE, sizeof(struct Particle) * header.npartTotal[1]);

		// copy particles to ocl device
		icl_write_buffer(particlesD, CL_TRUE, sizeof(struct Particle) * header.npartTotal[1], &particles[0], NULL, NULL);
		
		
		run(1, 0.00001, 0.0025, particles, particlesD, header.npartTotal[1], dev, &tree, ref);
		
		
		
//		icl_buffer* kdTree = icl_create_buffer(dev, CL_MEM_READ_WRITE, sizeof(struct KdNode) * (header.npartTotal[1] * 2 - 1));
//		buildTree(tree.nodelist, particlesD, kdTree, header.npartTotal[1], dev);
 
 		// TODO particles have been resorted during tree construction. Upload them in original sorting for comparison, not needed for correctness REMOVE IT!
//		icl_write_buffer(particlesD, CL_TRUE, sizeof(struct Particle) * header.npartTotal[1], &particles[0], NULL, NULL);

// 		walk(kdTree, particlesD, header.npartTotal[1], 0.00001, dev);
 
// 		// read particles from device, used as reference for correctness check
// 		icl_read_buffer(particlesD, CL_TRUE, sizeof(struct Particle) * header.npartTotal[1], &ref[0], NULL, NULL);
// 
// 		printf("Walk second time with last acceleration of particles\n");
// 
// 		walk(kdTree, particlesD, header.npartTotal[1], 0.00001, dev, particles);
// 
// 		// read particles from device
 		icl_read_buffer(particlesD, CL_TRUE, sizeof(struct Particle) * header.npartTotal[1], &particles[0], NULL, NULL);
//
//		icl_release_buffer(kdTree);
		icl_release_buffer(particlesD);
		icl_release_devices();

		printf("\nSUCCESS\n");
	} else {
		printf("ERROR! Cannot find requested device\n");
		return -1;
	}

//	check_force("forcetest_1e5.txt", "result.txt", particles, particleIds, header.npartTotal[1]);
//	check_force("forcetest.txt", "result.txt", particles, particleIds, header.npartTotal[1]);
#if timing == 1
	check_force_internal("result.txt", ref, particles, particleIds, header.npartTotal[1]);
#endif

// display interactions, stored in each particle at acc.x
/*FLOAT average = 0;
for(UINT i = 0; i < header.npartTotal[1]; ++i) {
	average += particles[i].acc.x;
}
printf("Average number of interactions: %f\n", average/header.npartTotal[1]);
*/
	free(tree.nodelist);
	free(particles);
	free(ref);
	

	return 0;
}


///////////////////////////////////////////////////////////////////////////////////////
//TIME INTEGRATION
///////////////////////////////////////////////////////////////////////////////////////

float calcTimestep(FLOAT eps, const FLOAT ErrTolIntAccuracy, icl_buffer* particles, const UINT nParticles, icl_device* dev)
{
#ifdef SOFTENING_PLUMMER
	eps = eps;
#elif defined(SOFTENING_SPLINEKERNEL)
	eps = 2.8f * eps; //used as a scale length for the spline kernel
#endif
	
	// compile opencl kernel
/*
	icl_kernel* tw = icl_create_kernel(dev, "kernel/timestep.cl", "calcTimestep", "-Iinclude", ICL_SOURCE);
	
	size_t localSize = 256; 
	size_t globalSize = ((nParticles + 255) / 256) * 256;
	icl_run_kernel(tw, 1, &globalSize, &localSize, NULL, NULL, 4,
					(size_t)0, (void *)particles,
					sizeof(UINT), &nParticles,
					sizeof(FLOAT), &eps,
					sizeof(FLOAT), &ErrTolIntAccuracy);

	icl_finish(dev);
	icl_release_kernel(tw);
	
	//TODO: implement individual timestepping (power of 2 timebins)
*/
	return 1e-5;
}
UINT compute_acceleration(UINT mode, icl_buffer* nodelist, icl_buffer* particles, const UINT nParticles, const FLOAT eps, UINT treeHeight, icl_device* dev, struct Particle* ref)
{
	static icl_buffer* kdTree = NULL;
		
	if(kdTree == NULL) //just for the first time to initialize acceleration for opening criterion
	{
		kdTree = icl_create_buffer(dev, CL_MEM_READ_WRITE, sizeof(struct KdNode) * (nParticles * 2 - 1));
		treeHeight = buildTree(nodelist, particles, kdTree, nParticles, dev);
		// first walk through the entire tree since acceleration is 0
		printf("First tree walk\n");
		walk(kdTree, particles, nParticles, eps, dev);
#if timing == 1
		// store acceleration of particles by unfolding the entire tree = direct force summation, used for correctness validation
		printf("Pruned tree walk\n");
                icl_read_buffer(particles, CL_TRUE, sizeof(struct Particle) * nParticles, &ref[0], NULL, NULL);
		// second walk with pruned tree
		walk(kdTree, particles, nParticles, eps, dev);
#endif
	}
	
	else if(mode == 1) //rebuild the tree
	{
		treeHeight = buildTree(nodelist, particles, kdTree, nParticles, dev);
		walk(kdTree, particles, nParticles, eps, dev);
	}
	else if(mode == 2) //TODO: dynamic tree update
	{
		updateTree(nodelist, particles, kdTree, nParticles, treeHeight, dev);
		walk(kdTree, particles, nParticles, eps, dev);
	}
	else //end of sim, release kdTree buffer
	{
		icl_release_buffer(kdTree);
	}

	return treeHeight;
}

void kick(const FLOAT dt, icl_buffer* particles, const UINT nParticles, icl_device* dev)
{
	icl_kernel* tw = icl_create_kernel(dev, "kernel/timestep.cl", "kick_particles", "-Iinclude", ICL_SOURCE);
	
	size_t localSize = 256; 
	size_t globalSize = ((nParticles + 255) / 256) * 256;
	icl_run_kernel(tw, 1, &globalSize, &localSize, NULL, NULL, 3,
					(size_t)0, (void *)particles,
					sizeof(UINT), &nParticles,
					sizeof(FLOAT), &dt);

	icl_finish(dev);
	icl_release_kernel(tw);
}

void drift(const FLOAT dt, icl_buffer* particles, const UINT nParticles, icl_device* dev)
{
	icl_kernel* tw = icl_create_kernel(dev, "kernel/timestep.cl", "drift_particles", "-Iinclude", ICL_SOURCE);
	
	size_t localSize = 256; 
	size_t globalSize = ((nParticles + 255) / 256) * 256;
	icl_run_kernel(tw, 1, &globalSize, &localSize, NULL, NULL, 3,
					(size_t)0, (void *)particles,
					sizeof(UINT), &nParticles,
					sizeof(FLOAT), &dt);

	icl_finish(dev);
	icl_release_kernel(tw);	
}

void energy_statistic(struct Particle *particles_host, const UINT nParticles, const FLOAT time)
{
	int i, j;
	double Ekin = 0.0;
	double Epot = 0.0;

	for(j = 0; j < nParticles; ++j)
	{
		if(j % 1000 == 0) {
			printf("calculating engergy statistic %d / %d\r", j, nParticles);
			fflush(stdout);
		}
		float v2 = particles_host[j].vel.x * particles_host[j].vel.x + particles_host[j].vel.y * particles_host[j].vel.y + particles_host[j].vel.z * particles_host[j].vel.z;
		Ekin += particles_host[j].mass * v2;
	
#pragma omp parallel for reduction(- : Epot)	
		for(i = 0; i < j; ++i)
			Epot -= particles_host[j].mass * particles_host[i].mass / dist(particles_host[j].posArr, particles_host[i].posArr);
	}
	
	Ekin *= 0.5;
	Epot *= G;
	
	FILE *fd = fopen("energy.txt", "a");
	fprintf(fd, "%g %g %g %g\n", time, Ekin, Epot, Ekin+Epot);
	printf("energy: %g %g %g %g                                    \n", time, Ekin, Epot, Ekin+Epot);
	fclose(fd);
}

void out_snapshot(struct Particle *particles_host, icl_buffer* particles_device, const UINT nParticles, icl_device* dev, const FLOAT current_time)
{
	int j;
	particle_data *P = (particle_data*)malloc(sizeof(particle_data) * nParticles);
	io_header header;
	static int cs = 0;
	char fn[200];

	sprintf(fn, "./output/snapshot_%03d", cs);
	
	icl_read_buffer(particles_device, CL_TRUE, sizeof(struct Particle) * nParticles, particles_host, NULL, NULL);
	energy_statistic(particles_host, nParticles, current_time);
///
//	printf("\tid: %d pos: %g %g %g vel: %g %g %g\n", particles_host[0].id, particles_host[0].pos.x, particles_host[0].pos.y, particles_host[0].pos.z, particles_host[0].vel.x, particles_host[0].vel.y, particles_host[0].vel.z);
///
	for(j = 0; j < nParticles; ++j)
	{
		P[j].Pos[0] = particles_host[j].pos.x;
		P[j].Pos[1] = particles_host[j].pos.y;
		P[j].Pos[2] = particles_host[j].pos.z;
		
		P[j].Vel[0] = particles_host[j].vel.x;
		P[j].Vel[1] = particles_host[j].vel.y;
		P[j].Vel[2] = particles_host[j].vel.z;
		
		//P[j].Mass = particles_host[j].mass; //in header
		P[j].Id = particles_host[j].id;
	
		P[j].Accel[0] = particles_host[j].acc.x;
		P[j].Accel[1] = particles_host[j].acc.y;
		P[j].Accel[2] = particles_host[j].acc.z;
	}
	
	memset(&header, 0, sizeof(io_header));
	header.time = current_time;
	header.num_files = 1;
	header.npart[1] = header.npartTotal[1] = nParticles;
	header.mass[1] = particles_host[0].mass; //for now all particles have the same mass
	
	unsigned blocks = 8199; //0b10000000000111;
	write_snapshot_format2(fn, &header, P, blocks);
	
	cs++;
	
	free(P);
}


//main simulation loop
void run(const FLOAT t_max, const FLOAT eps, const FLOAT ErrTolIntAccuracy, struct Particle *particles_host, icl_buffer* particles_device, const UINT nParticles, icl_device* dev, struct Tree *tree, struct Particle *ref)
{
	UINT k = 0;
	UINT treeHeight = 0;
	//FLOAT dt=1.5e-6;
	//FLOAT dt=1.220703125e-5;
	FLOAT dt=3.05176e-06;
	FLOAT current_time = 0.0; //time before current full timestep (drift), kicks are at current_time-+dt/2.0
	FLOAT timeBetSnapshot = 1e-3;
	FLOAT timeLastSnapshot = 0.0;


/*
	tree->nodelist[0].center_of_mass.x = 1;
	tree->nodelist[0].center_of_mass.y = 2;
	tree->nodelist[0].center_of_mass.z = 3;

	tree->nodelist[0].center_geometric.x = 1.3;
	tree->nodelist[0].center_geometric.y = 2.2;
	tree->nodelist[0].center_geometric.z = 3.1;

	tree->nodelist[0].mass = 77.0;
	tree->nodelist[0].l = 42.7;

	tree->nodelist[0].bounding_box.box[0].x = -3.0;
	tree->nodelist[0].bounding_box.box[0].y = -3.2;
	tree->nodelist[0].bounding_box.box[0].z = -3.3;
	tree->nodelist[0].bounding_box.box[1].x = 3.0;
	tree->nodelist[0].bounding_box.box[1].y = 3.2;
	tree->nodelist[0].bounding_box.box[1].z = 3.3;

	tree->nodelist[0].size = 8;
	tree->nodelist[0].level = 7;
	tree->nodelist[0].address = 17;

	tree->nodelist[0].left_child = 1;
	tree->nodelist[0].right_child = 2;
	tree->nodelist[0].split_dim = 3;
*/
	// create root node in nodelist
	tree->nodelist[0].particlesLow = 0;
	tree->nodelist[0].particlesHigh = nParticles;
	tree->nodelist[0].level = 0;
	//tree->nodelist[0].bounding_box = bounding_box;
	tree->nodelist[0].address = 0;

	icl_buffer* nodelist = icl_create_buffer(dev, CL_MEM_READ_WRITE, sizeof(struct Node) * (nParticles * 2 - 1));
	// copy root node to ocl device
	icl_write_buffer(nodelist, CL_TRUE, sizeof(struct Node), &(tree->nodelist[0]), NULL, NULL);

	//snapshot_000 IC with computed code properties (acceleration...)
//	out_snapshot(particles_host, particles_device, nParticles, dev, 0.0); 

	// TODO particles have been resorted during tree construction. Upload them in original sorting for comparison, not needed for correctness REMOVE IT!
//	icl_write_buffer(particles_device, CL_TRUE, sizeof(struct Particle) * nParticles, &particles_host[0], NULL, NULL);

	//just for the first time, kick to half timestep
	treeHeight = compute_acceleration(1, nodelist, particles_device, nParticles, eps, treeHeight, dev, ref);
#if timing == 1
	return;		
#endif
	//dt = calcTimestep(eps, ErrTolIntAccuracy, particles, nParticles);
	kick(dt/2.0, particles_device, nParticles, dev);



//	out_snapshot(particles_host, particles_device, nParticles, dev, 0.0); 



	//TEST
	icl_read_buffer(particles_device, CL_TRUE, sizeof(struct Particle) * nParticles, particles_host, NULL, NULL);
	printf("\tid: %d pos: %g %g %g vel: %g %g %g\n", particles_host[0].id,  particles_host[0].pos.x, particles_host[0].pos.y, particles_host[0].pos.z, particles_host[0].vel.x, particles_host[0].vel.y, particles_host[0].vel.z);
	printf("\tacc: %g %g %g\n", particles_host[0].acc.x, particles_host[0].acc.y, particles_host[0].acc.z);
	//
	
	while(current_time < t_max)
	{
		current_time += dt;
		printf("___step: %d time: %g\n", k++, current_time);
//		printf("\tid: %d pos: %g %g %g vel: %g %g %g\n", )
		
		//drift to next full timestep at current_time
		drift(dt, particles_device, nParticles, dev);
		
		//get new accelerations
		treeHeight = compute_acceleration(2, nodelist, particles_device, nParticles, eps, treeHeight, dev, ref); //TODO: mode 2: implement dynamic tree update
		
		//kick particles to current_time+dt/2.0
		kick(dt, particles_device, nParticles, dev);
		
		//output & energy statistic
		if(current_time-timeLastSnapshot > timeBetSnapshot)
		{
			out_snapshot(particles_host, particles_device, nParticles, dev, current_time);
			timeLastSnapshot = current_time;
		}

		//TEST
		icl_read_buffer(particles_device, CL_TRUE, sizeof(struct Particle) * nParticles, particles_host, NULL, NULL);
		printf("\tid: %d pos: %g %g %g vel: %g %g %g\n", particles_host[0].id,  particles_host[0].pos.x, particles_host[0].pos.y, particles_host[0].pos.z, particles_host[0].vel.x, particles_host[0].vel.y, particles_host[0].vel.z);
		printf("\tacc: %g %g %g\n", particles_host[0].acc.x, particles_host[0].acc.y, particles_host[0].acc.z);
	}
	
	out_snapshot(particles_host, particles_device, nParticles, dev, current_time); //write a snapshot also for the final time
	printf("final time reached: %g\n", current_time);
	
	compute_acceleration(0, nodelist, particles_device, nParticles, eps, treeHeight, dev, ref); //clean up
	
	
	icl_release_buffer(nodelist);

}





