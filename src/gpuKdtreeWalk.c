//#include "gpukdtree.h"
#include "gpuKdtreeWalk.h"

void walk(icl_buffer* kdTree, icl_buffer* particles, const UINT nParticles, FLOAT eps, icl_device* dev)
{
#ifdef SOFTENING_PLUMMER
	eps = eps;
#elif defined(SOFTENING_SPLINEKERNEL)
	eps = 2.8f * eps; //used as a scale length for the spline kernel
#endif

	// compile opencl kernel
	icl_kernel* tw = icl_create_kernel(dev, "kernel/kdTreeWalk.cl", "treeWalk", "-Iinclude", ICL_SOURCE);

	icl_timer* timer = icl_init_timer(ICL_MILLI);
	icl_start_timer(timer);
/*
	cl_int status = CL_SUCCESS;
	cl_image_format format;
	format.image_channel_order = CL_RGBA;
	format.image_channel_data_type = CL_FLOAT;
	cl_image_desc descriptor;
	descriptor.image_type = CL_MEM_OBJECT_IMAGE1D_BUFFER;
	descriptor.image_width = 3 * (nParticles * 2 - 1);
	descriptor.image_height = 1;
	descriptor.image_depth = 1;
	descriptor.image_array_size = 0;
	descriptor.image_row_pitch = 0;
	descriptor.image_slice_pitch = 0;
	descriptor.num_mip_levels = 0;
	descriptor.num_samples = 0;
	descriptor.buffer = getBuffer(kdTree);
	
	icl_buffer* image = icl_create_buffer(dev, CL_MEM_READ_WRITE, 1);
	icl_local_device* ld = getLocalDevice(dev->device_id);

	((icl_local_buffer*)(image->buffer_add))->mem = clCreateImage(ld->context, ICL_MEM_READ_ONLY ,
		&format, &descriptor, NULL, &status);

	ICL_ASSERT(status == CL_SUCCESS, "Error setting image: \"%s\"", icl_error_string(status));
*/

//	UINT cnt = 0;
	size_t localSize = 256; 
	size_t globalSize = ((nParticles + 255) / 256) * 256;

	icl_run_kernel(tw, 1, &globalSize, &localSize, NULL, NULL, 4,
					(size_t)0, (void *)kdTree,
//					(size_t)0, (void *)image,
					(size_t)0, (void *)particles,
					sizeof(UINT), &nParticles,
					sizeof(FLOAT), &eps);

	icl_finish(dev);

	double walk = icl_stop_timer(timer);
	printf("Walk: %f\n", walk);

	icl_release_kernel(tw);
	icl_release_timer(timer);
//	printf("calls %d\n", cnt);
}
