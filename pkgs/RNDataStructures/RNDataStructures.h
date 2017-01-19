/* Include file for GAPS data structures module */

#ifndef __RN__DATA_STRUCTURES__H__
#define __RN__DATA_STRUCTURES__H__

// define a custom assertion
#define rn_assertion(__x__) \
    if (!(__x__)) { fprintf(stderr, "Assertion error %s at line %d in file %s\n", #__x__, __LINE__, __FILE__); exit(-1); }


/* data structure include files */

#include "R3Graphics/R3Graphics.h"
#include "RNMisc.h"
#include "RNMinBinaryHeap.h"


/* Initialization functions */

int RNInitDataStructures(void);
void RNStopDataStructures(void);


#endif







