// Include file for the neuron feature class

#ifndef __NEURON_FEATURE__H__
#define __NEURON_FEATURE__H__


////////////////////////////////////////////////////////////////////////
// CLASS DEFINITION
////////////////////////////////////////////////////////////////////////

class NeuronFeature {
public:
    //////////////////////////////////
    //// CONSTRUCTORS/DESTRUCTORS ////
    //////////////////////////////////
    
    NeuronFeature(void);
    virtual ~NeuronFeature(void);
    
    
    ////////////////////////////
    //// PROPERTY FUNCTIONS ////
    ////////////////////////////
    
    virtual int CreateFeatures(void);
    
    
    ///////////////////////
    //// I/O FUNCTIONS ////
    ///////////////////////
    
    int ReadFeature(void);
    int WriteFeature(void);
    
    
protected:
    friend class NeuronClassifier;
    char feature_name[512];
    std::vector<RNScalar> attributes;
};

#endif