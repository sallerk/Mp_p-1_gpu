// Copyright (C) Mihai Preda (tinycl.h declarations), shim by this project.
//
// Resolves the OpenCL entry points from the driver's OpenCL.dll at run time, so
// the program needs no OpenCL SDK, no headers and no import library -- the same
// trick mersenne_tf uses. PRPLL's src/tinycl.h already declares the whole API as
// plain extern "C" with no calling convention, which is exactly right on x64
// Windows (there is only one calling convention there), so we can simply define
// those symbols here and forward each one on first use.
//
// Resolution is lazy and per-function: an entry point the driver does not export
// (clSVMAlloc on an OpenCL 1.2 stack, say) only becomes an error if it is
// actually called.

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>

#include <cstdio>
#include <cstdlib>

#include "tinycl.h"

namespace {

HMODULE clLib() {
  static HMODULE lib = [] {
    HMODULE h = LoadLibraryA("OpenCL.dll");
    if (!h) {
      fprintf(stderr,
              "Could not load OpenCL.dll.\n"
              "It ships with your GPU driver; install or update the driver.\n");
      exit(1);
    }
    return h;
  }();
  return lib;
}

void* clSym(const char* name) {
  void* p = (void*) GetProcAddress(clLib(), name);
  if (!p) {
    fprintf(stderr, "OpenCL.dll does not export %s -- driver too old?\n", name);
    exit(1);
  }
  return p;
}

} // namespace

// Defines an extern "C" function that forwards to the driver's implementation.
// PARAMS and ARGS are parenthesised so their commas survive macro expansion.
#define FWD(RET, NAME, PARAMS, ARGS)                        \
  RET NAME PARAMS {                                         \
    using Fn = RET (*) PARAMS;                              \
    static Fn fn = (Fn) clSym(#NAME);                       \
    return fn ARGS;                                         \
  }

extern "C" {

FWD(unsigned, clGetPlatformIDs,
    (unsigned a, cl_platform_id* b, unsigned* c), (a, b, c))
FWD(int, clGetPlatformInfo,
    (cl_platform_id a, cl_device_info b, size_t c, void* d, size_t* e), (a, b, c, d, e))
FWD(int, clGetDeviceIDs,
    (cl_platform_id a, cl_device_type b, unsigned c, cl_device_id* d, unsigned* e), (a, b, c, d, e))
FWD(int, clGetDeviceInfo,
    (cl_device_id a, cl_device_info b, size_t c, void* d, size_t* e), (a, b, c, d, e))

FWD(cl_context, clCreateContext,
    (const intptr_t* a, unsigned b, const cl_device_id* c,
     void (*d)(const char*, const void*, size_t, void*), void* e, int* f), (a, b, c, d, e, f))
FWD(int, clReleaseContext, (cl_context a), (a))

FWD(cl_command_queue, clCreateCommandQueueWithProperties,
    (cl_context a, cl_device_id b, const cl_queue_properties* c, int* d), (a, b, c, d))
FWD(int, clReleaseCommandQueue, (cl_command_queue a), (a))
FWD(int, clGetCommandQueueInfo,
    (cl_command_queue a, cl_command_queue_info b, size_t c, void* d, size_t* e), (a, b, c, d, e))

FWD(cl_program, clCreateProgramWithSource,
    (cl_context a, unsigned b, const char** c, const size_t* d, int* e), (a, b, c, d, e))
FWD(cl_program, clCreateProgramWithBinary,
    (cl_context a, unsigned b, const cl_device_id* c, const size_t* d,
     const unsigned char** e, int* f, int* g), (a, b, c, d, e, f, g))
FWD(int, clBuildProgram,
    (cl_program a, unsigned b, const cl_device_id* c, const char* d,
     void (*e)(cl_program, void*), void* f), (a, b, c, d, e, f))
FWD(int, clCompileProgram,
    (cl_program a, unsigned b, const cl_device_id* c, const char* d,
     unsigned e, const cl_program* f, const char* const* g,
     void (*h)(cl_program, void*), void* i), (a, b, c, d, e, f, g, h, i))
FWD(cl_program, clLinkProgram,
    (cl_context a, unsigned b, const cl_device_id* c, const char* d,
     unsigned e, const cl_program* f,
     void (*g)(cl_program, void*), void* h, int* i), (a, b, c, d, e, f, g, h, i))
FWD(int, clGetProgramBuildInfo,
    (cl_program a, cl_device_id b, cl_program_build_info c, size_t d, void* e, size_t* f), (a, b, c, d, e, f))
FWD(int, clGetProgramInfo,
    (cl_program a, cl_program_info b, size_t c, void* d, size_t* e), (a, b, c, d, e))
FWD(int, clReleaseProgram, (cl_program a), (a))

FWD(cl_kernel, clCreateKernel, (cl_program a, const char* b, int* c), (a, b, c))
FWD(int, clReleaseKernel, (cl_kernel a), (a))
FWD(int, clSetKernelArg, (cl_kernel a, unsigned b, size_t c, const void* d), (a, b, c, d))
FWD(int, clSetKernelArgSVMPointer, (cl_kernel a, unsigned b, const void* c), (a, b, c))
FWD(int, clGetKernelInfo,
    (cl_kernel a, cl_kernel_info b, size_t c, void* d, size_t* e), (a, b, c, d, e))
FWD(int, clGetKernelArgInfo,
    (cl_kernel a, unsigned b, cl_kernel_arg_info c, size_t d, void* e, size_t* f), (a, b, c, d, e, f))
FWD(int, clGetKernelWorkGroupInfo,
    (cl_kernel a, cl_device_id b, cl_kernel_work_group_info c, size_t d, void* e, size_t* f), (a, b, c, d, e, f))

FWD(cl_mem, clCreateBuffer,
    (cl_context a, cl_mem_flags b, size_t c, void* d, int* e), (a, b, c, d, e))
FWD(int, clReleaseMemObject, (cl_mem a), (a))

FWD(int, clEnqueueNDRangeKernel,
    (cl_command_queue a, cl_kernel b, unsigned c, const size_t* d, const size_t* e,
     const size_t* f, unsigned g, const cl_event* h, cl_event* i), (a, b, c, d, e, f, g, h, i))
FWD(int, clEnqueueReadBuffer,
    (cl_command_queue a, cl_mem b, cl_bool c, size_t d, size_t e, void* f,
     unsigned g, const cl_event* h, cl_event* i), (a, b, c, d, e, f, g, h, i))
FWD(int, clEnqueueWriteBuffer,
    (cl_command_queue a, cl_mem b, cl_bool c, size_t d, size_t e, const void* f,
     unsigned g, const cl_event* h, cl_event* i), (a, b, c, d, e, f, g, h, i))
FWD(int, clEnqueueCopyBuffer,
    (cl_command_queue a, cl_mem b, cl_mem c, size_t d, size_t e, size_t f,
     unsigned g, const cl_event* h, cl_event* i), (a, b, c, d, e, f, g, h, i))
FWD(int, clEnqueueFillBuffer,
    (cl_command_queue a, cl_mem b, const void* c, size_t d, size_t e, size_t f,
     unsigned g, const cl_event* h, cl_event* i), (a, b, c, d, e, f, g, h, i))
FWD(int, clEnqueueMarkerWithWaitList,
    (cl_command_queue a, unsigned b, const cl_event* c, cl_event* d), (a, b, c, d))

FWD(int, clFlush, (cl_command_queue a), (a))
FWD(int, clFinish, (cl_command_queue a), (a))

FWD(int, clReleaseEvent, (cl_event a), (a))
FWD(int, clWaitForEvents, (unsigned a, const cl_event* b), (a, b))
FWD(int, clGetEventInfo,
    (cl_event a, cl_event_info b, size_t c, void* d, size_t* e), (a, b, c, d, e))
FWD(int, clGetEventProfilingInfo,
    (cl_event a, cl_profiling_info b, size_t c, void* d, size_t* e), (a, b, c, d, e))

FWD(void*, clSVMAlloc, (cl_context a, cl_svm_mem_flags b, size_t c, unsigned d), (a, b, c, d))
FWD(void, clSVMFree, (cl_context a, void* b), (a, b))

} // extern "C"
