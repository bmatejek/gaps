// Source file for the neuron machine learning class



// include files

#include "NeuronML.h"



// private variables

static int neuronml_active_count = 0;



int NeuronMLInit(void)
{
  // check whether already initialized
  if ((neuronml_active_count++) > 0) return TRUE;

  // initialize dependencies
  if (!R3InitShapes()) return FALSE;

  // return OK status
  return TRUE;
}



void NeuronMLStop(void)
{
  // check whether already initialized
  if ((--neuronml_active_count) > 0) return;

  // stop dependencies
  RNStopDataStructures();
}
