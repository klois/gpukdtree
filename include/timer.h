#pragma once

/*
#ifdef __cplusplus
#define INLINE inline
#else
#define INLINE __inline
#endif
*/
#define INLINE

#ifndef _WIN32
	#include <time.h>
	#include <stdbool.h>
#else
	#define WINDOWS_LEAN_AND_MEAN
	#define NOMINMAX
	#include <windows.h>
#endif

typedef enum {ICL_SEC, ICL_MILLI, ICL_NANO} icl_time_flag;

#ifdef _WIN32
typedef __int64 time_int;
#else
typedef long long time_int;
#endif

typedef struct _icl_timer {
	time_int start;
	time_int clocks;
#ifdef _WIN32
	time_int freq;
#endif
	double current_time;
	icl_time_flag time_flag;
} icl_timer;


//#ifdef __cplusplus 
//#ifndef _WIN32 
//extern "C" {
//#endif
//#endif

icl_timer* icl_init_timer(icl_time_flag time_flag);
void icl_start_timer(icl_timer* timer);
void icl_restart_timer(icl_timer* timer);
double icl_stop_timer(icl_timer* timer);
void icl_release_timer(icl_timer* timer);

//#ifdef __cplusplus 
//#ifndef _WIN32 
//}
//#endif
//#endif
