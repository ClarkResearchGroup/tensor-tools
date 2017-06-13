#ifndef My_Timer
#define My_Timer

#include <sys/time.h>
#include <iostream>

struct timeval tv;
double get_clock() {
	struct timeval tv; int ok;
	ok = gettimeofday(&tv, 0);
	if (ok<0) { std::printf("gettimeofday error");}
	return (tv.tv_sec * 1.0 + tv.tv_usec * 1.0E-6);
}

#endif
