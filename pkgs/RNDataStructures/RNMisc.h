// Include file for miscellaneous useful functions

#ifndef __RN__MISC__H__
#define __RN__MISC__H__


#include "RNDataStructures.h"
#include <H5Cpp.h>


// function list
void RNProgressBar(int index, int nindices);
void RNDeflateIntegerArray(int *entires, int nentries);
void RNBestFitLine(RNScalar *x, RNScalar *y, int n, RNScalar& alpha, RNScalar& beta, RNScalar& RSquared);
R3Grid **RNReadH5File(const char *h5_filename, const char *dataset_name);
R3CharGrid *RNReadH5ImageFile(const char *h5_filename, const char *dataset_name);
int RNWriteH5File(R3Grid **grids, int ngrids, const char *h5_filename, const char *dataset_name, RNBoolean isInteger);

#endif
