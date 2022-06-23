/************************************************************
  
  This example shows how to write and read a hyperslab.  It 
  is derived from the h5_read.c and h5_write.c examples in 
  the "Introduction to HDF5". It works on the 3 dimensional data.

 ************************************************************/
 
#include "hdf5.h"


#define FILE        "sds.h5"
#define DATASETNAME "IntArray" 

#define RANK 1   //since it will work on  2 dimensional data
#define RANK_OUT 1

#define X 40



#define X1 40




int
main (int argc, char **argv)
{
    hsize_t     dimsf[1];              /* dataset dimensions */
    int         data[X];            /* data to write */

    
    /* 
     * Data  and output buffer initialization. 
     */
    hid_t       file, dataset;         /* handles */
    hid_t       dataspace;   
    hid_t       memspace; 
    hsize_t     dimsm[1];              /* memory space dimensions 1D*/
    hsize_t     dims_out[1];           /* dataset dimensions 1D */      
    herr_t      status;      
    hid_t	property_list_id_MPIO;                 /* property list identifier */   
    hid_t	data_transfer_propertylist;             
                    

    
    int data_out[X];   //data out 3d is 6x6x3
    int data_out1[X1];  //data out1  is 2x6x3  
     
    hsize_t     count[1];              /* size of the hyperslab in the file */
    hsize_t    offset[1];             /* hyperslab offset in the file */
    hsize_t     count_out[1];          /* size of the hyperslab in memory */
    hsize_t    offset_out[1];

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
    
    for (i = 0; i < X; i++)
                {data[i] = l;   //6x6x6
                l++;}
            
    if(mpi_rank==0){    

        for (i = 0; i < X; i++)
            printf("%5d", data[i]);
        printf("\n");
            
        } 
    
    //2d operations
    
    file = H5Fcreate_async (FILE, H5F_ACC_TRUNC, H5P_DEFAULT, property_list_id_MPIO,es_id);  
    
    

    dimsf[0]=X;
    
    
    dataspace = H5Screate_simple (RANK, dimsf, NULL); 
    
   
    dataset = H5Dcreate_async (file, DATASETNAME, H5T_STD_I32BE, dataspace, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT,es_id);   //H5T_STD_I32BE = 32-bit big-endian signed integers

    dimsm[0] = X1;
        
    memspace = H5Screate_simple (RANK_OUT, dimsm, NULL);   //RANK_OUT=3


    if(mpi_rank==0){

        
        
        offset[0] = 0;
       
        count[0]  = 5;
        
        status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);  
    
        printf("MPI rank=%d Hyperslab operation on dataspace using offset  %llu and count %llu \n",mpi_rank,offset[0],count[0]);
        offset_out[0] = 0;  //offset=0x0
        count_out[0]  = 5; //count_out=2 X 2   
        
        status = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
        
        printf("-----------------  on memory space using offset  %llu and count %llu  \n",offset_out[0],count_out[0]);

        status = H5Dwrite_async (dataset, H5T_NATIVE_INT, memspace, dataspace, data_transfer_propertylist, data,es_id);    
        
        
        offset[0] = 5;
       
        count[0]  = 3;
        
        status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);  
    
        printf("MPI rank=%d Hyperslab operation on dataspace using offset  %llu and count %llu \n",mpi_rank,offset[0],count[0]);
        offset_out[0] = 5;  //offset=0x0
        count_out[0]  = 3; //count_out=2 X 2   
        
        status = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
        
        printf("-----------------  on memory space using offset  %llu and count %llu  \n",offset_out[0],count_out[0]);

        status = H5Dwrite_async (dataset, H5T_NATIVE_INT, memspace, dataspace, data_transfer_propertylist, data,es_id);    
        
        
        offset[0] = 8;
        count[0]  = 2;
        
        status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);  
    
        printf("MPI rank=%d Hyperslab operation on dataspace using offset  %llu and count %llu \n",mpi_rank,offset[0],count[0]);
        offset_out[0] = 8;  //offset=0x0
        count_out[0]  = 2; //count_out=2 X 2   
        
        status = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
        
        printf("-----------------  on memory space using offset  %llu and count %llu  \n",offset_out[0],count_out[0]);

        status = H5Dwrite_async (dataset, H5T_NATIVE_INT, memspace, dataspace, data_transfer_propertylist, data,es_id);    
        
        
        
        offset[0] = 10;
        count[0]  = 10;
        
        status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);  
    
        printf("MPI rank=%d Hyperslab operation on dataspace using offset  %llu and count %llu \n",mpi_rank,offset[0],count[0]);
        offset_out[0] = 10;  //offset=0x0
        count_out[0]  = 10; //count_out=2 X 2   
        
        status = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
        
        printf("-----------------  on memory space using offset  %llu and count %llu  \n",offset_out[0],count_out[0]);

        status = H5Dwrite_async (dataset, H5T_NATIVE_INT, memspace, dataspace, data_transfer_propertylist, data,es_id);    
        
        
          
        
        

    }
    if(mpi_rank==1){

        
        offset[0] = 20;
        
        count[0]  = 5;
        status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);  
        printf("MPI rank=%d Hyperslab operation on dataspace using offset  %llu and count %llu \n",mpi_rank,offset[0],count[0]);

        offset_out[0] = 20;
        
        count_out[0]  = 5;
        status = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
        printf("-----------------  on memory space using offset  %llu and count %llu \n",offset_out[0],count_out[0]);
        
        
        // status = H5Dread_async (dataset, H5T_NATIVE_INT, memspace, dataspace, H5P_DEFAULT, data_out2,es_id);
        status = H5Dwrite_async (dataset, H5T_NATIVE_INT, memspace, dataspace, data_transfer_propertylist, data,es_id);    
        
          
        
        offset[0] = 25;
        
        count[0]  = 3;
        status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);  
        printf("MPI rank=%d Hyperslab operation on dataspace using offset  %llu and count %llu \n",mpi_rank,offset[0],count[0]);


        offset_out[0] = 25;
        
        count_out[0]  = 3;
        status = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
        printf("-----------------  on memory space using offset  %llu and count %llu \n",offset_out[0],count_out[0]);
        
        
        // status = H5Dread_async (dataset, H5T_NATIVE_INT, memspace, dataspace, H5P_DEFAULT, data_out2,es_id);
        status = H5Dwrite_async (dataset, H5T_NATIVE_INT, memspace, dataspace, data_transfer_propertylist, data,es_id);    
        
          
        
        offset[0] = 28;
        
        count[0]  = 2;
        status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);  
        printf("MPI rank=%d Hyperslab operation on dataspace using offset  %llu and count %llu \n",mpi_rank,offset[0],count[0]);


        offset_out[0] = 28;
        
        count_out[0]  = 2;
        status = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
        printf("-----------------  on memory space using offset  %llu and count %llu \n",offset_out[0],count_out[0]);
        
        
        // status = H5Dread_async (dataset, H5T_NATIVE_INT, memspace, dataspace, H5P_DEFAULT, data_out2,es_id);
        status = H5Dwrite_async (dataset, H5T_NATIVE_INT, memspace, dataspace, data_transfer_propertylist, data,es_id);    
        
          
        
        offset[0] = 30;
        
        count[0]  = 10;
        status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);  
        printf("MPI rank=%d Hyperslab operation on dataspace using offset  %llu and count %llu \n",mpi_rank,offset[0],count[0]);


        offset_out[0] = 30;
        
        count_out[0]  = 10;
        status = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
        printf("-----------------  on memory space using offset  %llu and count %llu \n",offset_out[0],count_out[0]);
        
        
        // status = H5Dread_async (dataset, H5T_NATIVE_INT, memspace, dataspace, H5P_DEFAULT, data_out2,es_id);
        status = H5Dwrite_async (dataset, H5T_NATIVE_INT, memspace, dataspace, data_transfer_propertylist, data,es_id);    
        
          
        
           }
            
           
        
            
            

            
    status = H5ESwait(es_id, H5ES_WAIT_FOREVER, &num_in_progress, &op_failed);
    if (status < 0) {
        fprintf(stderr, "Error with H5ESwait\n");
        
    }
    if (print_dbg_msg)
        fprintf(stderr, "H5ESwait done\n");
    
        MPI_Barrier(comm);

    if(mpi_rank==0){
        status = H5Dread_async (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_out,es_id);

        printf("\nData out from the file\n");
        
            
        status = H5ESwait(es_id, H5ES_WAIT_FOREVER, &num_in_progress, &op_failed);
        if (status < 0) {
            fprintf(stderr, "Error with H5ESwait\n");
            
        }
        for (i = 0; i < X; i++)
                printf(" %4d ", data_out[i]);
        printf("\n");
        
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
