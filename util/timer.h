#ifndef BASIC_TIMER
#define BASIC_TIMER

#include "types_and_headers.h"

// Timer
struct timeval _tval;
double get_clock() {
	struct timeval _tval; int ok;
	ok = gettimeofday(&_tval, 0);
	if (ok<0) { std::printf("gettimeofday error");}
	return (_tval.tv_sec * 1.0 + _tval.tv_usec * 1.0E-6);
}

#endif
