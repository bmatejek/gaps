// Include file for neuron machine learning class

#ifndef __NEURON_ML__H__
#define __NEURON_ML__H__



// external dependencies

#include <vector>



// dependency includes

#include "RNDataStructures/RNDataStructures.h"



// class declarations

class NeuronFeature;
class NeuronClassifier;



// NeuronML includes

#include "NeuronFeature.h"
// TODO add NeuronClassifier
//#include "NeuronClassifier.h"
//#include "NeuronRFClassifier.h"



// initialization functions

int NeuronMLInit(void);
void NeuronMLStop(void);



#endif
