/* Kernel code taken from "Scan Primitives for GPU Computing" (Sengupta et al.) */

#pragma once
#ifdef CUDA
#include "lib_icl_cuda.h"
#else
#include "lib_icl.h"
#endif



void segmented_scan_init(UINT maxN, icl_device *_dev, char* build_options, icl_create_kernel_flag flag);

void segmented_scan_release();

void scan_host(UINT *data, UINT *flag, UINT n);
void scan(icl_buffer *data, icl_buffer *flag, UINT n);

void recursive_scan(icl_buffer *d_data, icl_buffer *d_part, icl_buffer *d_flag, UINT n);


