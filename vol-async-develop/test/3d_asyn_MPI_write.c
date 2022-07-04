/************************************************************
  
  This example shows how to write and read a hyperslab.  It 
  is derived from the h5_read.c and h5_write.c examples in 
  the "Introduction to HDF5". It works on the 3 dimensional data.

 ************************************************************/
 
#include "hdf5.h"


#define FILE        "sds.h5"
#define DATASETNAME "IntArray" 

#define RANK 3   //since it will work on  3 dimensional data
#define RANK_OUT 3

#define X 6
#define Y 6
#define Z 3

#define X1 6
#define Y1 6
#define Z1 3

#define X2 6
#define Y2 6
#define Z2 3


int
main (int argc, char **argv)
{
    hsize_t     dimsf[3];              /* dataset dimensions */
    int         data[X][Y][Z];            /* data to write */

    
    /* 
     * Data  and output buffer initialization. 
     */
    hid_t       file, dataset;         /* handles */
    hid_t       dataspace;   
    hid_t       memspace; 
    hsize_t     dimsm[3];              /* memory space dimensions 1D*/
    hsize_t     dims_out[3];           /* dataset dimensions 1D */      
    herr_t      status;      
    hid_t	property_list_id_MPIO;                 /* property list identifier */   
    hid_t	data_transfer_propertylist;             
                    

    
    int data_out[X][Y][Z];   //data out 3d is 6x6x3
    int data_out1[X1][Y1][Z1];  //data out1  is 2x6x3  
    int data_out2[X2][Y2][Z2];  //data out1  is 3x6x3  
     
    hsize_t     count[3];              /* size of the hyperslab in the file */
    hsize_t    offset[3];             /* hyperslab offset in the file */
    hsize_t     count_out[3];          /* size of the hyperslab in memory */
    hsize_t    offset_out[3];

    int         i, j, k, status_n, rank;
    int print_dbg_msg=1;
    hbool_t op_failed;
    size_t  num_in_progress;
    hid_t   es_id = H5EScreate();
    int provided;
     /*
     * MPI variables
     */
    int mpi_size, mpi_rank;
    MPI_Comm comm  = MPI_COMM_WORLD;
    MPI_Info info  = MPI_INFO_NULL;

    /*
     * Initialize MPI
     */
    //MPI_Init(&argc, &argv);
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);  


    property_list_id_MPIO = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(property_list_id_MPIO, comm, info);


    data_transfer_propertylist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(data_transfer_propertylist, H5FD_MPIO_COLLECTIVE);
   

   

    //3d data
    int l=0;
    for (k = 0; k < Z; k++)
    for (i = 0; i < X; i++)
        for (j = 0; j < Y; j++)
                {data[i][j][k] = l;   //6x6x6
                l++;}
            
    if(mpi_rank==0){    
    for(k=0;k<Z;k++){
            printf ("\nFirst 3D Data Z=%d:\n ",k);
            
            
            for (i = 0; i < X; i++){ 
                for (j = 0; j < Y; j++) 
                    printf("%5d", data[i][j][k]);
            printf("\n ");
            }
        } 
    }
    //3d operations
    
    file = H5Fcreate_async (FILE, H5F_ACC_TRUNC, H5P_DEFAULT, property_list_id_MPIO,es_id);  
    /*
    * Describe the size of the array and create the data space for fixed
    * size dataset. 
    */
    

    dimsf[0]=X;
    dimsf[1]=Y;
    dimsf[2]=Z;
    dataspace = H5Screate_simple (RANK, dimsf, NULL); 
    
    /*
    * Create a new dataset within the file using defined dataspace and
    * default dataset creation properties.
    */
    dataset = H5Dcreate_async (file, DATASETNAME, H5T_STD_I32BE, dataspace, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT,es_id);   //H5T_STD_I32BE = 32-bit big-endian signed integers

    dimsm[0] = X1;
    dimsm[1] = Y1;
    dimsm[2] = Z1;

        
        
    memspace = H5Screate_simple (RANK_OUT, dimsm, NULL);   //RANK_OUT=3


    if(mpi_rank==0){

        
        
        offset[0] = 0;
        offset[1]=0;
        offset[2]=0;   // select 0x0x0
        count[0]  = 2;
        count[1]  = 2;
        count[2]  = 3;
        status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);  
    
        printf("MPI rank=%d Hyperslab operation on dataspace using offset  %llux%llux%llu and count %llux%llux%llu \n",mpi_rank,offset[0],offset[1],offset[2],count[0],count[1],count[2]);
        offset_out[0] = 0;  //offset=0x0x0
        offset_out[1] = 0;
        offset_out[2] = 0;
        count_out[0]  = 2; //count_out=2 X 2 x 3  
        count_out[1]  = 2;   
        count_out[2]  = 3;   
        status = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
        
        printf("-----------------  on memory space using offset  %llux%llux%llu and count %llux%llux%llu  \n",offset_out[0],offset_out[1],offset_out[2],count_out[0],count_out[1],count_out[2]);


        //status = H5Dread_async (dataset, H5T_NATIVE_INT, memspace, dataspace, H5P_DEFAULT, data_out1,es_id);
        status = H5Dwrite_async (dataset, H5T_NATIVE_INT, memspace, dataspace, data_transfer_propertylist, data,es_id);
        
        offset[0] = 0;
        offset[1]=2;
        offset[2]=0;   // select 0x2x0
        count[0]  = 2;
        count[1]  = 2;
        count[2]  = 3;
        status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);  
    
        printf("MPI rank=%d Hyperslab operation on dataspace using offset  %llux%llux%llu and count %llux%llux%llu \n",mpi_rank,offset[0],offset[1],offset[2],count[0],count[1],count[2]);
        offset_out[0] = 0;  //offset=0x2x0
        offset_out[1] = 2;
        offset_out[2] = 0;
        count_out[0]  = 2; //count_out=2 X 2 x 3  
        count_out[1]  = 2;   
        count_out[2]  = 3;   
        status = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
        
        printf("-----------------  on memory space using offset  %llux%llux%llu and count %llux%llux%llu  \n",offset_out[0],offset_out[1],offset_out[2],count_out[0],count_out[1],count_out[2]);


        //status = H5Dread_async (dataset, H5T_NATIVE_INT, memspace, dataspace, H5P_DEFAULT, data_out1,es_id);
        status = H5Dwrite_async (dataset, H5T_NATIVE_INT, memspace, dataspace, data_transfer_propertylist, data,es_id);
        
        
        offset[0] = 0;
        offset[1]=4;
        offset[2]=0;   // select 0x4x0
        count[0]  = 2;
        count[1]  = 2;
        count[2]  = 3;
        status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);  
    
        printf("MPI rank=%d Hyperslab operation on dataspace using offset  %llux%llux%llu and count %llux%llux%llu \n",mpi_rank,offset[0],offset[1],offset[2],count[0],count[1],count[2]);
        offset_out[0] = 0;  //offset=0x4x0
        offset_out[1] = 4;
        offset_out[2] = 0;
        count_out[0]  = 2; //count_out=2 X 2 x 3  
        count_out[1]  = 2;   
        count_out[2]  = 3;   
        status = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
        
        printf("-----------------  on memory space using offset  %llux%llux%llu and count %llux%llux%llu  \n",offset_out[0],offset_out[1],offset_out[2],count_out[0],count_out[1],count_out[2]);


        //status = H5Dread_async (dataset, H5T_NATIVE_INT, memspace, dataspace, H5P_DEFAULT, data_out1,es_id);
        status = H5Dwrite_async (dataset, H5T_NATIVE_INT, memspace, dataspace, data_transfer_propertylist, data,es_id);
        

    }
    if(mpi_rank==1){

        
        offset[0] = 2;
        offset[1] = 0;
        offset[2] = 0;    //offset=2x0x0
        count[0]  = 4;
        count[1]  = 2;  //count=4x2x3
        count[2]  = 3;
        status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);  
        printf("MPI rank=%d Hyperslab operation on dataspace using offset  %llux%llux%llu and count %llux%llux%llu \n",mpi_rank,offset[0],offset[1],offset[2],count[0],count[1],count[2]);

        offset_out[0] = 2;
        offset_out[1] = 0; //offset_out=2x0x0
        offset_out[2] = 0;
        count_out[0]  = 4;
        count_out[1]  = 2;  //count_out=4x2x3
        count_out[2]  = 3;
        status = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
        printf("-----------------  on memory space using offset  %llux%llux%llu and count %llux%llux%llu  \n",offset_out[0],offset_out[1],offset_out[2],count_out[0],count_out[1],count_out[2]);
        
        
    // status = H5Dread_async (dataset, H5T_NATIVE_INT, memspace, dataspace, H5P_DEFAULT, data_out2,es_id);
        status = H5Dwrite_async (dataset, H5T_NATIVE_INT, memspace, dataspace, data_transfer_propertylist, data,es_id);    
        
        
    //offset from 3x0x0 and count 3x2x3 
        offset[0] = 2;
        offset[1] = 2;
        offset[2] = 0;    //offset=2x2x0
        count[0]  = 4;
        count[1]  = 2;  //count=4x2x3
        count[2]  = 3;
        status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);  
        printf("MPI rank=%d Hyperslab operation on dataspace using offset  %llux%llux%llu and count %llux%llux%llu \n",mpi_rank,offset[0],offset[1],offset[2],count[0],count[1],count[2]);

        offset_out[0] = 2;
        offset_out[1] = 2; //offset_out=2x2x0
        offset_out[2] = 0;
        count_out[0]  = 4;
        count_out[1]  = 2;  //count_out=4x2x3
        count_out[2]  = 3;
        status = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
        printf("-----------------  on memory space using offset  %llux%llux%llu and count %llux%llux%llu  \n",offset_out[0],offset_out[1],offset_out[2],count_out[0],count_out[1],count_out[2]);
        
        
    // status = H5Dread_async (dataset, H5T_NATIVE_INT, memspace, dataspace, H5P_DEFAULT, data_out2,es_id);
        status = H5Dwrite_async (dataset, H5T_NATIVE_INT, memspace, dataspace, data_transfer_propertylist, data,es_id);    
        
        
        offset[0] = 2;
        offset[1] = 4;
        offset[2] = 0;    //offset=2x4x0
        count[0]  = 4;
        count[1]  = 2;  //count=4x2x3
        count[2]  = 3;
        status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);  
        printf("MPI rank=%d Hyperslab operation on dataspace using offset  %llux%llux%llu and count %llux%llux%llu \n",mpi_rank,offset[0],offset[1],offset[2],count[0],count[1],count[2]);

        offset_out[0] = 2;
        offset_out[1] = 4; //offset_out=2x4x0
        offset_out[2] = 0;
        count_out[0]  = 4;
        count_out[1]  = 2;  //count_out=4x2x3
        count_out[2]  = 3;
        status = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
        printf("-----------------  on memory space using offset  %llux%llux%llu and count %llux%llux%llu  \n",offset_out[0],offset_out[1],offset_out[2],count_out[0],count_out[1],count_out[2]);
        
        
    // status = H5Dread_async (dataset, H5T_NATIVE_INT, memspace, dataspace, H5P_DEFAULT, data_out2,es_id);
        status = H5Dwrite_async (dataset, H5T_NATIVE_INT, memspace, dataspace, data_transfer_propertylist, data,es_id);    
        
        
                    }
            
           
        
            
            

            
    status = H5ESwait(es_id, H5ES_WAIT_FOREVER, &num_in_progress, &op_failed);
    if (status < 0) {
        fprintf(stderr, "Error with H5ESwait\n");
        
    }
    if (print_dbg_msg)
        fprintf(stderr, "H5ESwait done\n");
    fflush(stdout);

    MPI_Barrier(comm);

    if(mpi_rank==0){
        status = H5Dread_async (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_out,es_id);
    
    
        status = H5ESwait(es_id, H5ES_WAIT_FOREVER, &num_in_progress, &op_failed);
        if (status < 0) {
            fprintf(stderr, "Error with H5ESwait\n");
            
        }
        printf("\nData out from the file\n");
        for(k=0;k<Z;k++){
            printf ("\n 3D Data Z=%d:\n ",k);
            
        
            for (i = 0; i < X; i++){ 
                for (j = 0; j < Y; j++) 
                    printf("%5d", data_out[i][j][k]);
            printf("\n ");
            }
        } 
    }
        /*
    * Close/release resources.
    */
    H5Sclose (dataspace);
    H5Sclose (memspace);
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
        
        
   

   
    status = H5ESclose(es_id);
    if (status < 0) {
        fprintf(stderr, "Can't close second event set\n");
        //ret = -1;
    }
    H5Pclose(property_list_id_MPIO);
    H5Pclose(data_transfer_propertylist);
    MPI_Barrier(comm);
    MPI_Finalize();

    return 0;

}     
