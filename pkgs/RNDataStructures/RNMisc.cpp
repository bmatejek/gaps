#include "RNDataStructures.h"



void
RNProgressBar(int index, int nindices)
{
   rn_assertion((0 <= index) && (index < nindices));

   if (index == 0) { printf("0%%"); fflush(stdout); return; }
   else if (index == nindices - 1) { printf("100%%"); fflush(stdout); return; }

   RNScalar percentage = (100 * (RNScalar)index) / nindices;
   RNScalar prev_percentage = (100 * (RNScalar)(index - 1)) / nindices;

   int percent = (int)(percentage);
   int prev_percent = (int)(prev_percentage);

   if (percent == prev_percent) return;

   if (percent % 10 == 0) {
      printf("%d%%", percent);
      fflush(stdout);
   }
   else if (percent % 2 == 0) {
      printf(".");
      fflush(stdout);
   }
}



void
RNDeflateIntegerArray(int *entries, int nentries)
{
   rn_assertion(entries != NULL);

   // find the maximum entry in entries
   int max_entry = -1;
   for (int ie = 0; ie < nentries; ++ie) {
      rn_assertion(entries[ie] >= 0);
      if (max_entry < entries[ie]) {
         max_entry = entries[ie];
      }
   }

   // find the entries that exist
   RNBoolean *entry_exists = new RNBoolean[max_entry + 1];
   for (int ie = 0; ie < max_entry + 1; ++ie) {
      entry_exists[ie] = FALSE;
   }

   // see if this entry exists
   for (int ie = 0; ie < nentries; ++ie) {
      int entry = entries[ie];
      entry_exists[entry] = TRUE;
   }

   // create mapping
   int *mapping = new int[max_entry + 1];
   int num_skipped = 0;
   for (int ie = 0; ie < max_entry + 1; ++ie) {
      if (entry_exists[ie]) mapping[ie] = ie - num_skipped;
      else { num_skipped++; mapping[ie] = -1; }
   }

   // reset entries' values
   for (int ie = 0; ie < nentries; ++ie) {
      entries[ie] = mapping[entries[ie]];
   }

   // free memory
   delete[] entry_exists;
   delete[] mapping;
}



void 
RNBestFitLine(RNScalar *x, RNScalar *y, int n, RNScalar& alpha, RNScalar& beta, RNScalar& rsquared)
{
   // just checking
   rn_assertion(x != NULL);
   rn_assertion(y != NULL);
   rn_assertion(n > 0);

   // calcluate the average of x and y
   RNScalar avgx = 0.0;
   RNScalar avgy = 0.0;
   for (int i = 0; i < n; ++i) {
      avgx += x[i];
      avgy += y[i];
   }
   avgx /= n;
   avgy /= n;
   
   // calculate beta
   RNScalar beta_numerator = 0.0;
   RNScalar beta_denominator = 0.0;
   for (int i = 0; i < n; ++i) {
      beta_numerator += (x[i] - avgx) * (y[i] - avgy);
      beta_denominator += (x[i] - avgx) * (x[i] - avgx);
   }

   beta = beta_numerator / beta_denominator;
   
   // calculate alpha from beta
   alpha = avgy - beta * avgx;
   
   // calculate the coefficient of determination
   RNScalar avgxy = 0.0;
   RNScalar avgxsquared = 0.0;
   RNScalar avgysquared = 0.0;
   for (int i = 0; i < n; ++i) {
      avgxy += x[i] * y[i];
      avgxsquared += x[i] * x[i];
      avgysquared += y[i] * y[i];
   }
   avgxy /= n;
   avgxsquared /= n;
   avgysquared /= n;

   // rxy
   rsquared = (avgxy - avgx * avgy) / sqrt((avgxsquared - avgx * avgx) * (avgysquared - avgy * avgy));
   
   // coefficient of determination
   rsquared = rsquared * rsquared;
}



R3Grid **RNReadH5File(const char *h5_filename, const char *dataset_name) 
{
    try {
        // open a read only file
        H5::H5File file(h5_filename, H5F_ACC_RDONLY);
        // read the dataset
        H5::DataSet dataset = file.openDataSet(dataset_name);
        // get the dataspace
        H5::DataSpace dataspace = dataset.getSpace();
        // get the number of dimensions in the dataspace
        int rank = dataspace.getSimpleExtentNdims();
        hsize_t dims[rank];
        int ndims = dataspace.getSimpleExtentDims(dims, NULL);

        int affinity_order = 0; // c, z, y, x
        if (ndims == 4 && dims[3] < 5) affinity_order = 1; // z, y, x, c

        long ngrids, zres, yres, xres;
        if (ndims == 2) {
          ngrids = 1;
          zres = 1;
          yres = dims[0];
          xres = dims[1];
        }
        else if (ndims == 3)
        {
            ngrids = 1;
            zres = dims[0];
            yres = dims[1];
            xres = dims[2];
        }
        else if (ndims == 4)
        {
            if (!affinity_order) {
                ngrids = dims[0];
                zres = dims[1];
                yres = dims[2];
                xres = dims[3];
            }
            else {
                zres = dims[0];
                yres = dims[1];
                xres = dims[2];
                ngrids = dims[3];
            }
        }
        else
        {
            printf("Unsupported number of dimensions: %d.\n", ndims);
            exit(-1);
        }

        R3Grid **grids = new R3Grid *[ngrids];
        for (int ig = 0; ig < ngrids; ++ig)
        {
            grids[ig] = new R3Grid(xres, yres, zres);
        }

        // get the class type
        H5T_class_t class_type = dataset.getTypeClass();

        if (class_type == H5T_INTEGER) {
            // read in the dataset
            int *data = new int[ngrids * zres * yres * xres];
            H5::DataSpace mem_space(rank, dims);
            dataset.read(data, H5::PredType::NATIVE_INT, mem_space, dataspace);

            int ii = 0;
            for (int ig = 0; ig < ngrids; ++ig) {
                for (int iz = 0; iz < zres; ++iz) {
                    for (int iy = 0; iy < yres; ++iy) {
                        for (int ix = 0; ix < xres; ++ix, ++ii) {
                            grids[ig]->SetGridValue(ix, iy, iz, data[ii]);
                        }
                    }
                }
            }
            
            // make sure not affinity order
            rn_assertion(!affinity_order);

            // free data
            delete data;
        }
        else if (class_type == H5T_FLOAT) {
            // read in the dataset
            float *data = new float[ngrids * zres * yres * xres];
            H5::DataSpace mem_space(rank, dims);
            dataset.read(data, H5::PredType::NATIVE_FLOAT, mem_space, dataspace);

            if (!affinity_order) {
                long ii = 0;
                for (long ig = 0; ig < ngrids; ++ig) {
                    for (long iz = 0; iz < zres; ++iz) {
                        for (long iy = 0; iy < yres; ++iy) {
                            for (long ix = 0; ix < xres; ++ix, ++ii) {
                                grids[ig]->SetGridValue(ix, iy, iz, data[ii]);
                            }
                        }
                    }
                }
            }
            else {
                long ii = 0;
                for (long iz = 0; iz < zres; ++iz) {
                    for (long iy = 0; iy < yres; ++iy) {
                        for (long ix = 0; ix < xres; ++ix) {
                            for (long ig = 0; ig < ngrids; ++ig, ++ii) {
                                grids[ig]->SetGridValue(ix, iy, iz, data[ii]);
                            }
                        }
                    }
                }
            }
            
            // free data
            delete data;
        }
        else {
            fprintf(stderr, "Unrecognized class type: %d\n", class_type);
            return NULL;
        }

        // return grids
        return grids;
    }
    catch (H5::FileIException error) {
        error.printError();
        return NULL;
    }
    catch (H5::DataSetIException error) {
        error.printError();
        return NULL;
    }
    catch (H5::DataSpaceIException error) {
        error.printError();
        return NULL;
    }
    catch (H5::DataTypeIException error) {
        error.printError();
        return NULL;
    }
    
    // should never reach here
    return NULL;
}



int RNWriteH5File(R3Grid **grids, int ngrids, const char *h5_filename, const char *dataset_name, RNBoolean isInteger)
{
    try {
        H5::Exception::dontPrint();

        // create a file
        H5::H5File file = H5::H5File(h5_filename, H5F_ACC_TRUNC);

        int rank;
        hsize_t *dims;
        if (ngrids == 1) {
            rank = 3;
            dims = new hsize_t[rank];
            dims[0] = grids[0]->ZResolution();
            dims[1] = grids[0]->YResolution();
            dims[2] = grids[0]->XResolution();
        }
        else {
            rank = 4;
            dims = new hsize_t[rank];
            dims[0] = ngrids;
            dims[1] = grids[0]->ZResolution();
            dims[2] = grids[0]->YResolution();
            dims[3] = grids[0]->XResolution();
        }

        H5::DataSpace dataspace(rank, dims);
        
        H5::DataSet dataset;
        if (isInteger) dataset = file.createDataSet(dataset_name, H5::PredType::NATIVE_INT, dataspace);
        else dataset = file.createDataSet(dataset_name, H5::PredType::NATIVE_FLOAT, dataspace);

        if (isInteger) {
            int *data = new int[ngrids * dims[0] * dims[1] * dims[2]];
            int ii = 0;
            for (int ig = 0; ig < ngrids; ++ig) {
                for (int iz = 0; iz < grids[ig]->ZResolution(); ++iz) {
                    for (int iy = 0; iy < grids[ig]->YResolution(); ++iy) {
                        for (int ix = 0; ix < grids[ig]->XResolution(); ++ix, ++ii) {
                            data[ii] = (int)(grids[ig]->GridValue(ix, iy, iz) + 0.5);
                        }
                    }
                }
            }
            dataset.write(data, H5::PredType::NATIVE_INT);
        
            delete[] data;
        }
        else {
            float *data = new float[ngrids * dims[0] * dims[1] * dims[2]];
            int ii = 0;
            for (int ig = 0; ig < ngrids; ++ig) {
                for (int iz = 0; iz < grids[ig]->ZResolution(); ++iz) {
                    for (int iy = 0; iy < grids[ig]->YResolution(); ++iy) {
                        for (int ix = 0; ix < grids[ig]->XResolution(); ++ix, ++ii) {
                            data[ii] = grids[ig]->GridValue(ix, iy, iz);
                        }
                    }
                }
            }
            dataset.write(data, H5::PredType::NATIVE_FLOAT);
        
            delete[] data;
        }

    }
    catch (H5::FileIException error) {
        error.printError();
        return 0;
    }
    catch (H5::DataSetIException error) {
        error.printError();
        return 0;
    }
    catch (H5::DataSpaceIException error) {
        error.printError();
        return 0;
    }

    return 1;
}
