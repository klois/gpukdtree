#include <stdlib.h>
#include "timer.h"

INLINE
icl_timer* icl_init_timer(icl_time_flag time_flag){
	icl_timer* timer = (icl_timer*)malloc(sizeof(icl_timer));
	timer->start = 0;
	timer->clocks = 0;
	timer->current_time = 0;
#ifdef _WIN32
//	timer->freq = 0;
	QueryPerformanceFrequency((LARGE_INTEGER *)&(timer->freq));
#endif
	timer->time_flag = time_flag;
	return timer;
}

INLINE
void icl_start_timer(icl_timer* timer){
#ifdef _WIN32
	QueryPerformanceCounter((LARGE_INTEGER *)&(timer->start));
#else
	struct timespec tmp;
	clock_gettime(CLOCK_REALTIME, &tmp);
	timer->start = (time_int)tmp.tv_sec * 1e9 + (time_int)tmp.tv_nsec;
#endif
}

INLINE
void icl_restart_timer(icl_timer* timer){
	icl_start_timer(timer);
	timer->clocks = 0;
	timer->current_time = 0;
}

INLINE
double icl_stop_timer(icl_timer* timer){
	time_int stop;
	double time = 0.0;

#ifdef _WIN32
	QueryPerformanceCounter((LARGE_INTEGER *)&stop);
#else
	struct timespec tmp;
	clock_gettime(CLOCK_REALTIME, &tmp);
	stop = (time_int)tmp.tv_sec * 1e9 + (time_int)tmp.tv_nsec;
#endif
	timer->clocks += (stop - timer->start);
	timer->start = 0;

	switch(timer->time_flag) {
		case ICL_NANO:
			#ifdef _WIN32
				time = (double)timer->clocks * (double) 1e9 / (double) timer->freq;
			#else
				time = (double)timer->clocks;
			#endif
			break;
		case ICL_MILLI:
			#ifdef _WIN32
				time = (double)timer->clocks * (double) 1e3 / (double) timer->freq;
			#else
				time = (double)timer->clocks / (double) 1e6;
			#endif
			break;
		case ICL_SEC:
			#ifdef _WIN32
				time = (double)timer->clocks / (double) timer->freq;
			#else
				time = (double)timer->clocks / (double) 1e9;
			#endif
			break;
	}
	
	timer->current_time = time;
	return time;
}

INLINE
void icl_release_timer(icl_timer* timer){
	free(timer);
}
