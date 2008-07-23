/*
 * Tü, 2006/01/29 HSGG
 */

#ifndef MYTIME_H
#define MYTIME_H

#include <time.h>

double timediff(const clock_t ta, const clock_t tb);

clock_t msg_starting();
void msg_finished(const clock_t ta);

#endif
