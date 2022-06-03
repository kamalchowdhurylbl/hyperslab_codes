/************************************************************
  
  This example shows how to write and read a hyperslab.  It 
  is derived from the h5_read.c and h5_write.c examples in 
  the "Introduction to HDF5".

 ************************************************************/
 
#include "hdf5.h"
#include "h5_async_lib.h"

#define FILE        "sds.h5"
#define DATASETNAME "IntArray" 
#define NX_SUB  3                      /* hyperslab dimensions */ 
#define NY_SUB  4 
#define NX 7                           /* output buffer dimensions */ 
#define NY 7 
#define NZ  3 
#define RANK         2
#define RANK_OUT     3

#define X     5                        /* dataset dimensions */
#define Y     6

int
main (void)
{
    hsize_t     dimsf[2];              /* dataset dimensions */
    int         data[X][Y];            /* data to write */

    /* 
     * Data  and output buffer initialization. 
     */
    hid_t       file, dataset;         /* handles */
    hid_t       dataspace;   
    hid_t       memspace; 
    hsize_t     dimsm[3];              /* memory space dimensions */
    hsize_t     dims_out[2];           /* dataset dimensions */      
    herr_t      status;                             

    int         data_out[NX][NY][NZ ]; /* output buffer */
    int data_out1[X][Y];
   
    hsize_t     count[2];              /* size of the hyperslab in the file */
    hsize_t    offset[2];             /* hyperslab offset in the file */
    hsize_t     count_out[3];          /* size of the hyperslab in memory */
    hsize_t    offset_out[3];         /* hyperslab offset in memory */
    int         i, j, k, status_n, rank;
    int print_dbg_msg=1;
    hbool_t op_failed;
    size_t  num_in_progress;
    hid_t   es_id = H5EScreate();

/*********************************************************  
   This writes data to the HDF5 file.  
 *********************************************************/  
 
    /* 
     * Data  and output buffer initialization. 
     */
    for (j = 0; j < X; j++) {
	for (i = 0; i < Y; i++)
	    data[j][i] = i + j;
    }     
    /*
     * 0 1 2 3 4 5 
     * 1 2 3 4 5 6
     * 2 3 4 5 6 7
     * 3 4 5 6 7 8
     * 4 5 6 7 8 9
     */

    /*
     * Create a new file using H5F_ACC_TRUNC access,
     * the default file creation properties, and the default file
     * access properties.
     */
    
    //file = H5Fcreate (FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);  
    file = H5Fcreate_async (FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT,es_id);  
     /* H5F_ACC_TRUNC -- If the file already exists, the file is opened with read-write access, 
        and new data will overwrite any existing data. 
        If the file does not exist, it is created and opened with read-write access.*/

    /*
     * Describe the size of the array and create the data space for fixed
     * size dataset. 
     */
    dimsf[0] = X;
    dimsf[1] = Y;
    dataspace = H5Screate_simple (RANK, dimsf, NULL); 

    /*
     * Create a new dataset within the file using defined dataspace and
     * default dataset creation properties.
     */
    dataset = H5Dcreate_async (file, DATASETNAME, H5T_STD_I32BE, dataspace,
                         H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT,es_id);   //H5T_STD_I32BE = 32-bit big-endian signed integers

    /*
     * Write the data to the dataset using default transfer properties.
     */
    status = H5Dwrite_async (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, data,es_id);  //H5T_NATIVE_INT = C-style int
   
    
   
    status = H5Dread_async (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                   H5P_DEFAULT, data_out1,es_id);

    
    status = H5ESwait(es_id, H5ES_WAIT_FOREVER, &num_in_progress, &op_failed);
    if (status < 0) {
        fprintf(stderr, "Error with H5ESwait\n");
        
    }
    if (print_dbg_msg)
        fprintf(stderr, "H5ESwait done\n");

    

     printf ("First Data:\n ");
    
    for (j = 0; j < X; j++) {
	for (i = 0; i < Y; i++) printf("%d ", data_out1[j][i]);
	printf("\n ");
    }
	printf("\n");
    /*
     * Close/release resources.
     */
    H5Sclose (dataspace);
   // H5Dclose (dataset);
    status = H5Dclose_async(dataset, es_id);
    if (status < 0) {
        fprintf(stderr, "Closing dataset failed\n");
        //ret = -1;
    }
   
   status = H5Fclose_async(file, es_id);
    if (status < 0) {
        fprintf(stderr, "Closing file failed\n");
        //ret = -1;
    } 
    //H5Fclose(file);
  status = H5ESwait(es_id, H5ES_WAIT_FOREVER, &num_in_progress, &op_failed);
    if (status < 0) {
        fprintf(stderr, "Error with H5ESwait\n");
        //ret = -1;
    } 

/*************************************************************  

  This reads the hyperslab from the sds.h5 file just 
  created, into a 2-dimensional plane of the 3-dimensional 
  array.

 ************************************************************/  

    for (j = 0; j < NX; j++) {
	for (i = 0; i < NY; i++) {
	    for (k = 0; k < NZ ; k++)
		data_out[j][i][k] = 0;
	}
    } 
 
    /*
     * Open the file and the dataset.
     */
    file = H5Fopen_async (FILE, H5F_ACC_RDONLY, H5P_DEFAULT,es_id); //H5F_ACC_RDONLY= An existing file is opened with read-only access. 
                                                        //If the file does not exist, H5Fopen fails. (Default)
    dataset = H5Dopen_async(file, DATASETNAME,H5P_DEFAULT,es_id);  //#define DATASETNAME "IntArray" 

    dataspace = H5Dget_space_async (dataset,es_id);    /* dataspace handle */
    rank      = H5Sget_simple_extent_ndims (dataspace);
    status_n  = H5Sget_simple_extent_dims (dataspace, dims_out, NULL);
    printf("\nRank: %d\nDimensions: %lu x %lu \n", rank,
	   (unsigned long)(dims_out[0]), (unsigned long)(dims_out[1]));

    /* 
     * Define hyperslab in the dataset. 
     */
    
    offset[0] = 1;
    offset[1] = 2; // offset=1x2
    count[0]  = NX_SUB;
    count[1]  = NY_SUB; //count=3x4
    status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);  
    //H5S_SELECT_SET=Replace the existing selection with the parameters from this call. Overlapping blocks are not supported with this operator.
    //H5Sselect_hyperslab selects a hyperslab region to 
    //add to the current selected region for a specified dataspace.
    //offset value from row 1, col 2


    /*
     * Define the memory dataspace.
     */
    dimsm[0] = NX;
    dimsm[1] = NY;
    dimsm[2] = NZ;  //dimsm =7x7x3
    memspace = H5Screate_simple (RANK_OUT, dimsm, NULL);   //RANK_OUT=3, RANK=2
    // 7x7x3 array 

    /* 
     * Define memory hyperslab. 
     */
    offset_out[0] = 3;
    //offset_out[0] = 0;
    offset_out[1] = 0;  // offset_out = 3 X 0 X 0
    offset_out[2] = 0;
    count_out[0]  = NX_SUB; //count_out=3 X 4 X 1 
    count_out[1]  = NY_SUB;
    count_out[2]  = 1;
    //count_out[2]  = 0;
    status = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset_out, NULL, 
                                  count_out, NULL);
                                  //H5Sselect_hyperslab selects a hyperslab region to add to the current selected region 
                                  //for the dataspace specified by space_id.

    /*
     * Read data from hyperslab in the file into the hyperslab in 
     * memory and display.
     */
    status = H5Dread_async (dataset, H5T_NATIVE_INT, memspace, dataspace,
                      H5P_DEFAULT, data_out,es_id);

    status = H5ESwait(es_id, H5ES_WAIT_FOREVER, &num_in_progress, &op_failed);
    if (status < 0) {
        fprintf(stderr, "Error with H5ESwait\n");
        
    }
    if (print_dbg_msg)
        fprintf(stderr, "H5ESwait done\n");


    printf ("Data:\n ");
    for (j = 0; j < NX; j++) {
	for (i = 0; i < NY; i++) printf("%d ", data_out[j][i][0]);
	printf("\n ");
    }
	printf("\n");
    /*
     * 0 0 0 0 0 0 0
     * 0 0 0 0 0 0 0
     * 0 0 0 0 0 0 0
     * 3 4 5 6 0 0 0  
     * 4 5 6 7 0 0 0
     * 5 6 7 8 0 0 0
     * 0 0 0 0 0 0 0
     */
    
    
    //H5ESclose(es_id);
    
    /*
     * Close and release resources.
     */
    //H5Dclose(dataset);
    
    //H5Dclose_async(dataset,es_id);
    status = H5Dclose_async(dataset, es_id);
    if (status < 0) {
        fprintf(stderr, "Closing dataset failed\n");
        //ret = -1;
    }
    /* status = H5ESwait(es_id, H5ES_WAIT_FOREVER, &num_in_progress, &op_failed);
    if (status < 0) {
        fprintf(stderr, "Error with H5ESwait\n");
        //ret = -1;
    } */
    H5Sclose (dataspace);
    H5Sclose (memspace);

    //H5Fclose_async(file,es_id);
    status = H5Fclose_async(file, es_id);
    if (status < 0) {
        fprintf(stderr, "Closing file failed\n");
        //ret = -1;
    }
    
    status = H5ESwait(es_id, H5ES_WAIT_FOREVER, &num_in_progress, &op_failed);
    if (status < 0) {
        fprintf(stderr, "Error with H5ESwait\n");
        //ret = -1;
    } 
    /*
     * Open the file and the dataset.
     */
    file = H5Fopen_async (FILE, H5F_ACC_RDONLY, H5P_DEFAULT,es_id); //H5F_ACC_RDONLY= An existing file is opened with read-only access. 
                                                        //If the file does not exist, H5Fopen fails. (Default)
    dataset = H5Dopen_async(file, DATASETNAME,H5P_DEFAULT,es_id);  //#define DATASETNAME "IntArray" 

    dataspace = H5Dget_space_async (dataset,es_id);    /* dataspace handle */
    rank      = H5Sget_simple_extent_ndims (dataspace);
    status_n  = H5Sget_simple_extent_dims (dataspace, dims_out, NULL);
    printf("\nRank: %d\nDimensions: %lu x %lu \n", rank,
	   (unsigned long)(dims_out[0]), (unsigned long)(dims_out[1]));

    /* 
     * Define hyperslab in the dataset. 
     */
    
    offset[0] = 2;
    offset[1] = 3; // offset=2x3
    count[0]  = 3;
    count[1]  = 3; //count=3x3
    status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);  
    //H5S_SELECT_SET=Replace the existing selection with the parameters from this call. Overlapping blocks are not supported with this operator.
    //H5Sselect_hyperslab selects a hyperslab region to 
    //add to the current selected region for a specified dataspace.
    //offset value from row 1, col 2


    /*
     * Define the memory dataspace.
     */
    dimsm[0] = NX;
    dimsm[1] = NY;
    dimsm[2] = NZ;  //dimsm =7x7x3
    memspace = H5Screate_simple (RANK_OUT, dimsm, NULL);   //RANK_OUT=3, RANK=2
    // 7x7x3 array 

    /* 
     * Define memory hyperslab. 
     */
    offset_out[0] = 0;
    //offset_out[0] = 0;
    offset_out[1] = 4;  // offset_out = 0 X 4 X 0
    offset_out[2] = 0;
    count_out[0]  = 3; //count_out=3 X 3 X 1 
    count_out[1]  = 3;
    count_out[2]  = 1;
    //count_out[2]  = 0;
    status = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset_out, NULL, 
                                  count_out, NULL);
                                  
    
                                  
    /*
     * Read data from hyperslab in the file into the hyperslab in 
     * memory and display.
     */
    status = H5Dread_async (dataset, H5T_NATIVE_INT, memspace, dataspace,
                      H5P_DEFAULT, data_out,es_id);

    status = H5ESwait(es_id, H5ES_WAIT_FOREVER, &num_in_progress, &op_failed);
    if (status < 0) {
        fprintf(stderr, "Error with H5ESwait\n");
        
    }
    if (print_dbg_msg)
        fprintf(stderr, "H5ESwait done\n");


    printf ("Data with z=0:\n ");
    for (j = 0; j < NX; j++) {
	for (i = 0; i < NY; i++) printf("%d ", data_out[j][i][0]);
	printf("\n ");
    }
	printf("\n");
    /*
     * 0 0 0 0 5 6 7
     * 0 0 0 0 6 7 8
     * 0 0 0 0 7 8 9
     * 3 4 5 6 0 0 0  
     * 4 5 6 7 0 0 0
     * 5 6 7 8 0 0 0
     * 0 0 0 0 0 0 0
     */
    

    offset_out[0] = 0;
    //offset_out[0] = 0;
    offset_out[1] = 4;  // offset_out = 0 X 4 X 1
    offset_out[2] = 1;
    count_out[0]  = 3; //count_out=3 X 3 X 1 
    count_out[1]  = 3;
    count_out[2]  = 1;
    //count_out[2]  = 0;
    status = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset_out, NULL, 
                                  count_out, NULL);
                                  
    /*
     * Read data from hyperslab in the file into the hyperslab in 
     * memory and display.
     */
    status = H5Dread_async (dataset, H5T_NATIVE_INT, memspace, dataspace,
                      H5P_DEFAULT, data_out,es_id);

    status = H5ESwait(es_id, H5ES_WAIT_FOREVER, &num_in_progress, &op_failed);
    if (status < 0) {
        fprintf(stderr, "Error with H5ESwait\n");
        
    }
    if (print_dbg_msg)
        fprintf(stderr, "H5ESwait done\n");


   
    /*
     * 0 0 0 0 5 6 7
     * 0 0 0 0 6 7 8
     * 0 0 0 0 7 8 9
     * 0 0 0 0 0 0 0  
     * 0 0 0 0 0 0 0
     * 0 0 0 0 0 0 0
     * 0 0 0 0 0 0 0
     */
    printf ("Data with Z=1:\n ");
    for (j = 0; j < NX; j++) {
	for (i = 0; i < NY; i++) printf("%d ", data_out[j][i][1]);
	printf("\n ");
    }
	printf("\n");
    
    //H5ESclose(es_id);
    
    /*
     * Close and release resources.
     */
    //H5Dclose(dataset);
    
    //H5Dclose_async(dataset,es_id);
    status = H5Dclose_async(dataset, es_id);
    if (status < 0) {
        fprintf(stderr, "Closing dataset failed\n");
        //ret = -1;
    }
    /* status = H5ESwait(es_id, H5ES_WAIT_FOREVER, &num_in_progress, &op_failed);
    if (status < 0) {
        fprintf(stderr, "Error with H5ESwait\n");
        //ret = -1;
    } */
    H5Sclose (dataspace);
    H5Sclose (memspace);

    //H5Fclose_async(file,es_id);
    status = H5Fclose_async(file, es_id);
    if (status < 0) {
        fprintf(stderr, "Closing file failed\n");
        //ret = -1;
    }
    
    status = H5ESwait(es_id, H5ES_WAIT_FOREVER, &num_in_progress, &op_failed);
    if (status < 0) {
        fprintf(stderr, "Error with H5ESwait\n");
        //ret = -1;
    }
    status = H5ESclose(es_id);
    if (status < 0) {
        fprintf(stderr, "Can't close second event set\n");
        //ret = -1;
    }
    

}     
