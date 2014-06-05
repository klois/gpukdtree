#include "device.h"
#include "structs.h"

/**
 * Copies smalllistsize (sizes[2]) to activelistsize (sizes[0]). Sets maxlevel (size[2]). Depending on setMaxLevel, smalllistsize and activelistsize are also
 * swapped and new maxlevel sizes[4] is set
 * @param sizes array continging
 * 		0 activelistsize
 * 		1 total number of nodes
 * 		2 smalllistsize
 * 		3 old maxlevel
 * 		4 new maxlevel
 * @param setMaxlevel flag and maximum tree level. If 0 size[2] will be reset to 0 and sizes[4] will be copied in sizes[3]. If bigger 0, sizes[2] will be set
 * 		previous value of sizes[0] and sizes[3] as well as sizes[4] are set to setMaxLevel
 */
__kernel void swapActiveAndSmalllist(__global UINT* sizes, UINT setMaxLevel)
{
	// swap size of activelist and smalllist
	if(!setMaxLevel)
	{
		sizes[0] = sizes[2];
		sizes[1] += sizes[0];
		sizes[2] = 0;
		// copy the new size to the old size
		sizes[3] = sizes[4];
	}
	else
	{
		UINT tmp = sizes[2];
		sizes[2] = sizes[0];
		sizes[0] = tmp;
		sizes[3] = setMaxLevel;
		sizes[4] = setMaxLevel;
	}
}
