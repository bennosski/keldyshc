#include "hdf5.h"
#include "constants.h"
#include <string.h>
#include <unistd.h>
#include <stdio.h>

void dsave(const char * filename, const char * dsetname, const double * restrict dset_data, int len)
{
  hid_t file_id, dataset_id, dataspace_id;
  //hsize_t * dims;
  herr_t status;

  if( access(filename, F_OK) != -1)
  {
    file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
  }
  else
  {
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  }

  //file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  /* Create the data space for the dataset. */
  //int i;
  //for(i=0; i<ndims; i++)
  //  dims[i] = (hsize_t)dimensions[i];
  
  int ndims = 1;
  hsize_t dims[] = {(hsize_t)len};
  
  dataspace_id = H5Screate_simple(ndims, dims, NULL);

  /* Create the dataset. */
  dataset_id = H5Dcreate2(file_id, dsetname, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /* write the data */
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);

  /* End access to the dataset and release resources used by it. */
  status = H5Dclose(dataset_id);

  /* Terminate access to the data space. */ 
  status = H5Sclose(dataspace_id);

  /* Close the file. */
  status = H5Fclose(file_id);
}

void zsave(const char * filename, const char * dsetname, const cdouble * restrict dset_data, int len)
{

  hid_t file_id, dataset_id, dataspace_id, memtype, filetype;
  //hsize_t * dims;
  herr_t status;

  if( access(filename, F_OK) != -1)
  {
    file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
  }
  else
  {
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  }
 
  /* Create the data space for the dataset. */
  //int i;
  //for(i=0; i<ndims; i++)
  //  dims[i] = (hsize_t)dimensions[i];

  int ndims = 1;
  hsize_t dims[] = {(hsize_t)len};

  /* Create compound datatype for memory */

  memtype = H5Tcreate(H5T_COMPOUND, sizeof(cdouble));
  status = H5Tinsert(memtype, "real", 0, H5T_NATIVE_DOUBLE);
  status = H5Tinsert(memtype, "imag", sizeof(double), H5T_NATIVE_DOUBLE);

  /* Create the compound datatype for the file */
    
  dataspace_id = H5Screate_simple(ndims, dims, NULL);

  /* Create or load the dataset. */
  /*
  if(H5Lexists(file_id, dsetname, H5P_DEFAULT)==1)
  {
    printf("dataset exists %s\n", dsetname);
    dataset_id = H5Dopen2(file_id, dsetname, H5P_ACC_RDWR);
  }
  else
  {
    printf("dataset does not exist %s\n", dsetname);
    dataset_id = H5Dcreate(file_id, dsetname, memtype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }
  */

  dataset_id = H5Dcreate(file_id, dsetname, memtype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /* write the data */
  status = H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);

  /* End access to the dataset and release resources used by it. */
  status = H5Dclose(dataset_id);

  /* Terminate access to the data space. */ 
  status = H5Sclose(dataspace_id);

  /* CHECK: DO I NOT NEED TO FREE memtype????
     It creates an error when I try to: */
 
  status = H5Tclose(memtype);

  /* Close the file. */
  status = H5Fclose(file_id);  
}




