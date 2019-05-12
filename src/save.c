#include "hdf5.h"

void save(const char * filename, const char * dsetname, const double * dset_data, const int * dimensions, int ndims)
{
  hid_t file_id, dataset_id, dataspace_id;
  hsize_t * dims;
  herr_t status;

  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  /* Create the data space for the dataset. */
  int i;
  for(i=0; i<ndims; i++)
    dims[i] = (hsize_t)dimensions[i];
  
  dataspace_id = H5Screate_simple(ndims, dims, NULL);

  /* Create the dataset. */
  dataset_id = H5Dcreate2(file_id, dsetname, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);

  /* End access to the dataset and release resources used by it. */
  status = H5Dclose(dataset_id);

  /* Terminate access to the data space. */ 
  status = H5Sclose(dataspace_id);

  /* Close the file. */
  status = H5Fclose(file_id);
}
