// Include file for minimum binary heap

#ifndef __RN__MIN__BINARY__HEAP__H__
#define __RN__MIN__BINARY__HEAP__H__



/////////////////////////////////////////////////////////////////////
// CLASS DEFINITION
/////////////////////////////////////////////////////////////////////

template <class PtrType>
class RNMinBinaryHeap {
public:
    //////////////////////////////////
	//// CONSTRUCTORS/DESTRUCTORS ////
    //////////////////////////////////

    RNMinBinaryHeap(PtrType base, RNScalar *value_ptr, int nentries);
    ~RNMinBinaryHeap();


    /////////////////////////////
	//// ATTRIBUTE FUNCTIONS ////
    /////////////////////////////

    int Size(void) const;
    int MinIndex(void) const;
    PtrType MinKey(void) const;
	PtrType KeyOf(int i) const;
    

    ////////////////////////////
    //// PROPERTY FUNCTIONS ////
    ////////////////////////////

    int IsEmpty(void) const;
    int Contains(int i) const;

    
    ////////////////////////////////
    //// MANIPULATION FUNCTIONS ////
    ////////////////////////////////
	
	// insert functions
    void Insert(int i, PtrType key);

    // key manipulation functions
	void ChangeKey(int i, PtrType key);
	void DecreaseKey(int i, PtrType key);
	void IncreaseKey(int i, PtrType key);
	
	// deletion functions
	PtrType DeleteMin(void);
    void Delete(int i);


    ////////////////////////////////////////////////////////////////////////
    // INTERNAL STUFF BELOW HERE
    ////////////////////////////////////////////////////////////////////////
   
private:
	// value property functions
	RNScalar Value(int i) const;
    int Compare(int i, int j) const;
    int Greater(int i, int j) const;

	// heap manipulation helper functions
    void Exch(int i, int j);
    void Swim(int k);
    void Sink(int k);

private:
    // instance variables
    int NMAX;
    int N;
    int *pq;
    int *qp;
    PtrType *keys;
    int value_offset;
};



// include files
#include "RNMinBinaryHeap.cpp"



#endif