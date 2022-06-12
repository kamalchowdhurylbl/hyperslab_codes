/*  
 *  This example writes data to the HDF5 file by rows.
 *  Number of processes is assumed to be 1 or multiples of 2 (up to 8)
 */
 
#include "hdf5.h"
#include "stdlib.h"


#define H5FILE_NAME     "SDS_row.h5"
#define DATASETNAME 	"IntArray" 
#define NX     8                      /* dataset dimensions */
#define NY     5 
#define RANK   2

int
main (int argc, char **argv)
{
    /*
     * HDF5 APIs definitions
     */ 	
    hid_t       file_id, dset_id;         /* file and dataset identifiers */
    hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
    hsize_t     dimsf[2];                 /* dataset dimensions */
    int         *data,*data_out;                    /* pointer to data buffer to write */
    hsize_t	count[2];	          /* hyperslab selection parameters */
    hsize_t	offset[2];
    hid_t	plist_id;                 /* property list identifier */
    int         i;
    herr_t	status;
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
 
    /* 
     * Set up file access property list with parallel I/O access
     */
     plist_id = H5Pcreate(H5P_FILE_ACCESS);
     H5Pset_fapl_mpio(plist_id, comm, info);

    /*
     * Create a new file collectively and release property list identifier.
     */
    file_id = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);
   

    /*
     * Create the dataspace for the dataset.
     */
    dimsf[0] = NX;
    dimsf[1] = NY; //8x5  rank=2
    filespace = H5Screate_simple(RANK, dimsf, NULL); 

    /*
     * Create the dataset with default properties and close filespace.
     */
    dset_id = H5Dcreate(file_id, DATASETNAME, H5T_NATIVE_INT, filespace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

    /* 
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    count[0] = dimsf[0]/mpi_size;
    count[1] = dimsf[1];
    offset[0] = mpi_rank * count[0];
    offset[1] = 0;
    memspace = H5Screate_simple(RANK, count, NULL);

    printf("mpi_rank=%d \n count[0]=%ld count[1]=%ld offset[0]=%ld offset[1]=%ld \n",mpi_rank,count[0],count[1],offset[0],offset[1]);

    /*
     * Select hyperslab in the file.
     */
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    /*
     * Initialize data buffer 
     */
    data = (int *) malloc(sizeof(int)*count[0]*count[1]);
    data_out = (int *) malloc(sizeof(int)*count[0]*count[1]);
    for (i=0; i < count[0]*count[1]; i++) {
        data[i] = mpi_rank + 10;
    }
   /*
         10, 10, 10, 10, 10,
         10, 10, 10, 10, 10,
         11, 11, 11, 11, 11,
         11, 11, 11, 11, 11,
         12, 12, 12, 12, 12,
         12, 12, 12, 12, 12,
         13, 13, 13, 13, 13,
         13, 13, 13, 13, 13
   */
    /*
     * Create property list for collective dataset write.
     */
   plist_id = H5Pcreate(H5P_DATASET_XFER);
   H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    
    status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace,
		      plist_id, data);
    
    free(data);

    /*
     * Close/release resources.
     */
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
   H5Pclose(plist_id);
    H5Fclose(file_id);

     if(mpi_rank==0){

            file_id = H5Fopen(H5FILE_NAME, H5F_ACC_RDONLY, H5P_DEFAULT);

            /*
            * Open the dataset.
            */
            dset_id = H5Dopen(file_id,DATASETNAME,H5P_DEFAULT);

            /*
            * Read data back to the buffer matrix.
            */
            status = H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, data_out);  //H5Dread=Reads raw data from a dataset into a buffer.

            for (i=0; i < count[0]*count[1]; i++) {
                printf(" %d ",data_out[i]);

            }

            H5Dclose(dset_id);
            
            H5Fclose(file_id);
    }
 
    MPI_Finalize();

    return 0;
}     
