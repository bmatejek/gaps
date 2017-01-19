/* Source file for GAPS data structures module */



/* Include files */

#include "RNBasics/RNBasics.h"



/* Private variables */

static int RNdata_structures_active_count = 0;



int RNInitDataStructures(void)
{
    // Check whether are already initialized 
    if ((RNdata_structures_active_count++) > 0) return TRUE;

    // Return OK status 
    return TRUE;
}



void RNStopDataStructures(void)
{
    // Check whether have been initialized 
    if ((--RNdata_structures_active_count) > 0) return;
}




