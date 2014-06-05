#include "device.h"
#include "structs.h"

/**
 * packs nextlist to remove all 0 values from the array
 * @param nextlist holds the indices in nodeslist that will be active in the next iteration + some zero values
 * @param sizes array with the number of active nodes at posititon 0. This value will be set to the numbe of active nodes in the next iteration
 */
__kernel void packNextlist(__global NodeId* nextlist, __global UINT* sizes) {
	UINT nNext = sizes[0] * 2;
	UINT cnt = 0;
//printf("-> %d %d\n", nextlist[0], nextlist[1]);
	while(nextlist[cnt] != 0 && cnt < nNext) {
		++cnt;
	}

	if(cnt == nNext) {
		sizes[0] = nNext; // save size of nextlist
		return;
	}

	for(UINT i = cnt+1; i < nNext; ++i)
	{
		UINT curr = nextlist[i];
		if(curr != 0)
		{
//printf("%d -> %d\n", i, cnt);
			nextlist[cnt] = curr;
			++cnt;
		}
	}

//	printf("activesize %d, nodesize %d, smallsize %d\n", sizes[0], sizes[1], sizes[2]);

	sizes[0] = cnt; // save size of nextlist

}
