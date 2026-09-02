#include <stdio.h>
#include <unistd.h>

// These run on TBB worker threads, with sibling workers still inside CUDA calls on
// other GPUs.  exit() would run the destructors of the global thrust::device_vectors
// underneath them, so the real message is followed by a "CUDA free failed:
// cudaErrorCudartUnloading" cascade and a SIGABRT that buries it -- exactly what the
// 14of22 seed-buffer overflow looked like in the field.  _exit() skips that; stderr
// is flushed explicitly because _exit() does not.
static inline void die_after(int code) {
    fflush(stderr);
    _exit(code);
}

// wrap of cudaSetDevice error checking in one place.  
static inline void check_cuda_setDevice(int device_id, const char* tag) {
    cudaError_t err = cudaSetDevice(device_id);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaSetDevice failed for device %d in %s failed with error \" %s \" \n", device_id, tag, cudaGetErrorString(err));
        die_after(11);
    }
}

// wrap of cudaMalloc error checking in one place.  
static inline void check_cuda_malloc(void** buf, size_t bytes, const char* tag) {
    cudaError_t err = cudaMalloc(buf, bytes);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMalloc of %lu bytes for %s failed with error \" %s \" \n", bytes, tag, cudaGetErrorString(err));
        die_after(12);
    }
}
	 
// wrap of cudaMemcpy error checking in one place.  
static inline void check_cuda_memcpy(void* dst_buf, void* src_buf, size_t bytes, cudaMemcpyKind kind, const char* tag) {
    cudaError_t err = cudaMemcpy(dst_buf, src_buf, bytes, kind);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaMemcpy of %lu bytes for %s failed with error \" %s \" \n", bytes, tag, cudaGetErrorString(err));
        die_after(13);
    }
}
	 
// wrap of cudaFree error checking in one place.  
static inline void check_cuda_free(void* buf, const char* tag) {
    cudaError_t err = cudaFree(buf);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: cudaFree for %s failed with error \" %s \" \n", tag, cudaGetErrorString(err));
        die_after(14);
    }
}
