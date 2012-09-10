#include "atomic_gridding_kernels.cu"
#include "../std_gridding_kernels.cu"
#include "cuda_utils.hpp"
#include "gridding_gpu.hpp"
#include "cufft_config.hpp"
/** gridding3D_gpu
  * forward gridding from image to grid/k-space
  * TODO
  * NFFT
**/
void gridding3D_gpu(CufftType**	data,			//kspace data array 
					int			data_count,		//data count, samples per trajectory
					int			n_coils,		//number of coils 
					DType*		crds,			//
					DType*		imdata,			//
					int			imdata_count,	//			
					int			grid_width,		//
					DType*		kernel,			//
					int			kernel_count,	//
					int			kernel_width,	//
					int*		sectors,		//
					int			sector_count,	//
					int*		sector_centers,	//
					int			sector_width,	//
					int			im_width,		//
					DType		osr,			//
					const GriddingOutput gridding_out)
{
	showMemoryInfo();
	GriddingInfo* gi_host = initAndCopyGriddingInfo(sector_count,sector_width,kernel_width,kernel_count,grid_width,im_width,osr);

	//cuda mem allocation
	DType *imdata_d, *crds_d, *kernel_d;//, *temp_gdata_d;
	CufftType *gdata_d, *data_d;
	int* sector_centers_d, *sectors_d;
	if (DEBUG)
		printf("allocate and copy imdata of size %d...\n",2*imdata_count*n_coils);
	allocateAndCopyToDeviceMem<DType>(&imdata_d,imdata,2*imdata_count*n_coils);

	if (DEBUG)
		printf("allocate and copy gdata of size %d...\n",gi_host->grid_width_dim );

	allocateAndSetMem<CufftType>(&gdata_d, gi_host->grid_width_dim,0);

	if (DEBUG)
		printf("allocate and copy data of size %d...\n",data_count * n_coils);
	allocateDeviceMem<CufftType>(&data_d,data_count * n_coils);

	if (DEBUG)
		printf("allocate and copy coords of size %d...\n",3*data_count);
	allocateAndCopyToDeviceMem<DType>(&crds_d,crds,3*data_count);
	
	if (DEBUG)
		printf("allocate and copy kernel in const memory of size %d...\n",kernel_count);
	HANDLE_ERROR(cudaGetSymbolAddress((void**)&kernel_d, KERNEL));
	HANDLE_ERROR(cudaMemcpyToSymbol(KERNEL,(void*)kernel,kernel_count*sizeof(DType)));
	//allocateAndCopyToDeviceMem<DType>(&kernel_d,kernel,kernel_count);
	if (DEBUG)
		printf("allocate and copy sectors of size %d...\n",sector_count+1);
	allocateAndCopyToDeviceMem<int>(&sectors_d,sectors,sector_count+1);
	if (DEBUG)
		printf("allocate and copy sector_centers of size %d...\n",3*sector_count);
	allocateAndCopyToDeviceMem<int>(&sector_centers_d,sector_centers,3*sector_count);
	if (DEBUG)
		printf("sector pad width: %d\n",gi_host->sector_pad_width);
	
	//Inverse fft plan and execution
	cufftHandle fft_plan;
	if (DEBUG)
		printf("creating cufft plan with %d,%d,%d dimensions\n",gi_host->grid_width,gi_host->grid_width,gi_host->grid_width);
	cufftResult res = cufftPlan3d(&fft_plan, gi_host->grid_width,gi_host->grid_width,gi_host->grid_width, CufftTransformType) ;
	if (res != CUFFT_SUCCESS) 
		printf("error on CUFFT Plan creation!!! %d\n",res);
	int err;

	//iterate over coils and compute result
	for (int coil_it = 0; coil_it < n_coils; coil_it++)
	{
		int data_coil_offset = coil_it * data_count;
		int im_coil_offset = 2 * coil_it * imdata_count;//gi_host->width_dim;
		//reset temp array
		//cudaMemset(temp_gdata_d,0, sizeof(DType)*temp_grid_count);
		cudaMemset(data_d,0, sizeof(CufftType)*data_count);
		cudaMemset(gdata_d,0, sizeof(CufftType)*gi_host->grid_width_dim);
		
		if (DEBUG && (cudaThreadSynchronize() != cudaSuccess))
			printf("error at thread synchronization 1: %s\n",cudaGetErrorString(cudaGetLastError()));
		// apodization Correction
		performForwardDeapodization(imdata_d + im_coil_offset,gi_host);
		
		if (DEBUG && (cudaThreadSynchronize() != cudaSuccess))
			printf("error at thread synchronization 2: %s\n",cudaGetErrorString(cudaGetLastError()));
		// resize by oversampling factor and zero pad
		performPadding(imdata_d + im_coil_offset,gdata_d,gi_host);
		
		if (DEBUG && (cudaThreadSynchronize() != cudaSuccess))
			printf("error at thread synchronization 3: %s\n",cudaGetErrorString(cudaGetLastError()));
		// shift image to get correct zero frequency position
		performFFTShift(gdata_d,FORWARD,gi_host->grid_width);
		
		if (DEBUG && (cudaThreadSynchronize() != cudaSuccess))
			printf("error at thread synchronization 4: %s\n",cudaGetErrorString(cudaGetLastError()));
		// eventually free imdata_d
		// Forward FFT to kspace domain
		if (err=pt2CufftExec(fft_plan, gdata_d, gdata_d, CUFFT_FORWARD) != CUFFT_SUCCESS)
		{
			printf("cufft has failed with err %i \n",err);
			showMemoryInfo(true);
		}
		
 		if (DEBUG && (cudaThreadSynchronize() != cudaSuccess))
			printf("error at thread synchronization 5: %s\n",cudaGetErrorString(cudaGetLastError()));
		performFFTShift(gdata_d,FORWARD,gi_host->grid_width);
		
		if (DEBUG && (cudaThreadSynchronize() != cudaSuccess))
			printf("error at thread synchronization 6: %s\n",cudaGetErrorString(cudaGetLastError()));
		// convolution and resampling to non-standard trajectory
		performForwardConvolution(data_d,crds_d,gdata_d,kernel_d,sectors_d,sector_centers_d,gi_host);
		//check for errors
		if (DEBUG && (cudaThreadSynchronize() != cudaSuccess))
			printf("error at thread synchronization 7: %s\n",cudaGetErrorString(cudaGetLastError()));
		//get result
	  	copyFromDevice<CufftType>(data_d, *data + data_coil_offset,data_count);
	}//iterate over coils
	cufftDestroy(fft_plan);
	// Destroy the cuFFT plan.
	if (DEBUG && (cudaThreadSynchronize() != cudaSuccess))
		printf("error at thread synchronization 8: %s\n",cudaGetErrorString(cudaGetLastError()));
	freeTotalDeviceMemory(data_d,crds_d,gdata_d,imdata_d,sectors_d,sector_centers_d,NULL);//NULL as stop
	
	if (DEBUG && (cudaThreadSynchronize() != cudaSuccess))
		printf("error at thread synchronization 9: %s\n",cudaGetErrorString(cudaGetLastError()));
  free(gi_host);
}

/** gridding3D_gpu
  * adjoint gridding from k-space to grid
  * TODO
  * NFFT^H
**/
void gridding3D_gpu_adj(DType*		data,			//kspace data array 
						int			data_count,		//data count, samples per trajectory
						int			n_coils,		//number of coils 
						DType*		crds,			//
						CufftType**	imdata,			//
						int			imdata_count,	//			
						int			grid_width,		//
						DType*		kernel,			//
						int			kernel_count,	//
						int			kernel_width,	//
						int*		sectors,		//
						int			sector_count,	//
						int*		sector_centers,	//
						int			sector_width,	//
						int			im_width,		//
						DType		osr,			//
						const GriddingOutput gridding_out)
{
	assert(sectors != NULL);
	
	showMemoryInfo();

	//split and run sectors into blocks
	//and each data point to one thread inside this block 
	GriddingInfo* gi_host = initAndCopyGriddingInfo(sector_count,sector_width,kernel_width,kernel_count,grid_width,im_width,osr);
	
	DType* data_d, *kernel_d;
  DType* crds_d;
	CufftType *gdata_d, *imdata_d;
	int* sector_centers_d, *sectors_d;

	if (DEBUG)
		printf("allocate and copy imdata of size %d...\n",imdata_count);
	allocateAndCopyToDeviceMem<CufftType>(&imdata_d,*imdata,imdata_count);//Konvention!!!

	if (DEBUG)
		printf("allocate and copy gdata of size %d...\n",gi_host->grid_width_dim);
	allocateDeviceMem<CufftType>(&gdata_d,gi_host->grid_width_dim);

	if (DEBUG)
		printf("allocate and copy data of size %d...\n",2*data_count*n_coils);
	allocateAndCopyToDeviceMem<DType>(&data_d,data,2*data_count*n_coils);

	if (DEBUG)
		printf("allocate and copy coords of size %d...\n",3*data_count);
	allocateAndCopyToDeviceMem<DType>(&crds_d,crds,3*data_count);
	
	if (DEBUG)
		printf("allocate and copy kernel in const memory of size %d...\n",kernel_count);
	HANDLE_ERROR(cudaGetSymbolAddress((void**)&kernel_d, "KERNEL"));
	HANDLE_ERROR(cudaMemcpyToSymbol("KERNEL",kernel,kernel_count*sizeof(DType)));
	//allocateAndCopyToDeviceMem<DType>(&kernel_d,kernel,kernel_count);
	if (DEBUG)
		printf("allocate and copy sectors of size %d...\n",sector_count+1);
	allocateAndCopyToDeviceMem<int>(&sectors_d,sectors,sector_count+1);
	if (DEBUG)
		printf("allocate and copy sector_centers of size %d...\n",3*sector_count);
	allocateAndCopyToDeviceMem<int>(&sector_centers_d,sector_centers,3*sector_count);
	if (DEBUG)
		printf("sector pad width: %d\n",gi_host->sector_pad_width);
	
	//Inverse fft plan and execution
	cufftHandle fft_plan;
	if (DEBUG)
		printf("creating cufft plan with %d,%d,%d dimensions\n",gi_host->grid_width,gi_host->grid_width,gi_host->grid_width);
	cufftResult res = cufftPlan3d(&fft_plan, gi_host->grid_width,gi_host->grid_width,gi_host->grid_width, CufftTransformType) ;
	if (res != CUFFT_SUCCESS) 
		printf("error on CUFFT Plan creation!!! %d\n",res);
	int err;

	//iterate over coils and compute result
	for (int coil_it = 0; coil_it < n_coils; coil_it++)
	{
		int data_coil_offset = 2 * coil_it * data_count;
		int im_coil_offset = coil_it * imdata_count;//gi_host->width_dim;
		cudaMemset(gdata_d,0, sizeof(CufftType)*gi_host->grid_width_dim);
		
		if (DEBUG && (cudaThreadSynchronize() != cudaSuccess))
			printf("error at adj thread synchronization 1: %s\n",cudaGetErrorString(cudaGetLastError()));
		performConvolution(data_d+data_coil_offset,crds_d,gdata_d,kernel_d,sectors_d,sector_centers_d,NULL,gi_host);

		if (DEBUG && (cudaThreadSynchronize() != cudaSuccess))
			printf("error at adj  thread synchronization 2: %s\n",cudaGetErrorString(cudaGetLastError()));
		if (gridding_out == CONVOLUTION)
		{
			if (DEBUG)
				printf("stopping output after CONVOLUTION step\n");
			//get output
			copyFromDevice<CufftType>(gdata_d,*imdata,gi_host->grid_width_dim);
			if (DEBUG)
				printf("test value at point zero: %f\n",(*imdata)[0].x);
			freeTotalDeviceMemory(data_d,crds_d,imdata_d,gdata_d,sectors_d,sector_centers_d,NULL);//NULL as stop token

			free(gi_host);
			/* Destroy the cuFFT plan. */
			cufftDestroy(fft_plan);
			return;
		}

		if (DEBUG && (cudaThreadSynchronize() != cudaSuccess))
			printf("error at adj  thread synchronization 3: %s\n",cudaGetErrorString(cudaGetLastError()));
		performFFTShift(gdata_d,INVERSE,gi_host->grid_width);
		//Inverse FFT
		if (err=pt2CufftExec(fft_plan, gdata_d, gdata_d, CUFFT_INVERSE) != CUFFT_SUCCESS)
		{
			printf("cufft has failed at adj with err %i \n",err);
			showMemoryInfo(true);
		}
	
		if (DEBUG && (cudaThreadSynchronize() != cudaSuccess))
			printf("error at adj  thread synchronization 4: %s\n",cudaGetErrorString(cudaGetLastError()));
		if (gridding_out == FFT)
		{
			if (DEBUG)
				printf("stopping output after FFT step\n");
			//get output
			copyFromDevice<CufftType>(gdata_d,*imdata,gi_host->grid_width_dim);
			
			//free memory
			if (cufftDestroy(fft_plan) != CUFFT_SUCCESS)
				printf("error on destroying cufft plan\n");
			freeTotalDeviceMemory(data_d,crds_d,imdata_d,gdata_d,sectors_d,sector_centers_d,NULL);//NULL as stop token
			free(gi_host);
			// Destroy the cuFFT plan.
			return;
		}

		if (DEBUG && (cudaThreadSynchronize() != cudaSuccess))
			printf("error at adj  thread synchronization 5: %s\n",cudaGetErrorString(cudaGetLastError()));
		performFFTShift(gdata_d,INVERSE,gi_host->grid_width);

		if (DEBUG && (cudaThreadSynchronize() != cudaSuccess))
			printf("error at adj  thread synchronization 6: %s\n",cudaGetErrorString(cudaGetLastError()));
		performCrop(gdata_d,imdata_d,gi_host);
		
		if (DEBUG && (cudaThreadSynchronize() != cudaSuccess))
			printf("error at adj  thread synchronization 7: %s\n",cudaGetErrorString(cudaGetLastError()));
		performDeapodization(imdata_d,gi_host);
		//check for errors
		if (DEBUG && (cudaThreadSynchronize() != cudaSuccess))
			printf("error: at adj  thread synchronization 8: %s\n",cudaGetErrorString(cudaGetLastError()));
		
		//get result
		copyFromDevice<CufftType>(imdata_d,*imdata+im_coil_offset,imdata_count);
	}//iterate over coils

	if (DEBUG && (cudaThreadSynchronize() != cudaSuccess))
		printf("error: at adj  thread synchronization 9: %s\n",cudaGetErrorString(cudaGetLastError()));
	cufftDestroy(fft_plan);
	// Destroy the cuFFT plan.

	if (DEBUG && (cudaThreadSynchronize() != cudaSuccess))
		printf("error: at adj  thread synchronization 10: %s\n",cudaGetErrorString(cudaGetLastError()));
	freeTotalDeviceMemory(data_d,crds_d,gdata_d,imdata_d,sectors_d,sector_centers_d,NULL);//NULL as stop
	free(gi_host);
}
