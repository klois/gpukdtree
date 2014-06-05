gpukdtree
=========

OpenCL:

Export OPENCL_ROOT to the root directory of your OpenCL installation.
Call make from the the root directory of the git repository
run with bin/gpukdtree [input_file] from the root directory of the git repository

CUDA:

Export CUDA_PATH to the root directory of your CUDA installation.
Call make cuda from the the root directory of the git repository
run with bin/gpukdtree [input_file] from the root directory of the git repository

Example input files can be found in the folder IC

To benchmark a single iteration of tree construction/walk set "#define timing 1" in includes/gpukdtree.h. To do a full nbody simulation set "#define timing 0"
