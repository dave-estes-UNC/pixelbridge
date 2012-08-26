#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/time.h>
#include <OpenGL/OpenGL.h>

#include "ClNddiDisplay.h"

using namespace nddi;

// public

ClNddiDisplay::ClNddiDisplay(std::vector<unsigned int> frameVolumeDimensionalSizes,
                             int inputVectorSize) {
	ClNddiDisplay(frameVolumeDimensionalSizes, 320, 240, inputVectorSize);
}

ClNddiDisplay::ClNddiDisplay(std::vector<unsigned int> frameVolumeDimensionalSizes,
                             int displayWidth, int displayHeight,
                             int inputVectorSize) {
	
	frameVolumeDimensionalSizes_ = frameVolumeDimensionalSizes;
	displayWidth_ = displayWidth;
	displayHeight_ = displayHeight;
	
    // Create the CostModel
    costModel = new CostModel();

	// Setup Input Vector
	clInputVector_ = new ClInputVector(costModel, inputVectorSize);
    // TODO(CDE): NULL this out
    inputVector_ = (InputVector*)clInputVector_;
	
	// Setup coefficient plane with zeroed coefficient matrices
    clCoefficientPlane_ = new ClCoefficientPlane(costModel, displayWidth_, displayHeight_, CM_WIDTH, CM_HEIGHT);
    // TODO(CDE): NULL this out
	coefficientPlane_ = (ClCoefficientPlane*)clCoefficientPlane_;
	
	// Setup framevolume and initialize to black
    clFrameVolume_ = new ClFrameVolume(costModel, frameVolumeDimensionalSizes);
    // TODO(CDE): NULL this out
	frameVolume_ = (FrameVolume*)clFrameVolume_;
	
    // We won't be using the reference components or frameBuffer_, so make them null.
    //inputVector_ = NULL;
	//coefficientPlane_ = NULL;
	//frameVolume_ = NULL;
    frameBuffer_ = NULL;

    // Initialize GL Texture
    InitializeGl();
	
	// Initialize CL
	InitializeCl();
}

ClNddiDisplay::~ClNddiDisplay() {
	
	delete(clInputVector_);
	delete(clFrameVolume_);
	delete(clCoefficientPlane_);

    Cleanup(false);
}

void ClNddiDisplay::Cleanup(bool shouldExit)
{
    if (clQueue_ != 0)
        clReleaseCommandQueue(clQueue_);
    
    if (clKernel_ != 0)
        clReleaseKernel(clKernel_);
    
    if (clProgram_ != 0)
        clReleaseProgram(clProgram_);
    
    if (clContext_ != 0)
        clReleaseContext(clContext_);
    
	if( clKernel_ != 0 ) 
		clReleaseKernel(clKernel_);
    
    if (clFrameVolumeDims_ != 0)
        clReleaseMemObject(clFrameVolumeDims_);

    if (clFrameBuffer_ != 0)
        clReleaseMemObject(clFrameBuffer_);

	if( texture_ != 0 ) 
	{
		glBindBuffer(GL_TEXTURE_RECTANGLE_ARB, texture_ );
		glDeleteBuffers(1, &texture_);
	}
    if (shouldExit)
        exit(0);
}

// Private

void ClNddiDisplay::InitializeGl() {

    texture_ = NULL;
    
    // Allocate texture as a CL-friendly image
    glGenTextures( 1, &texture_ );
    glTexEnvi( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );
    glBindTexture(GL_TEXTURE_RECTANGLE_ARB, texture_);
    glTexImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, GL_RGBA32F_ARB, displayWidth_,
                 displayHeight_, 0, GL_LUMINANCE, GL_FLOAT, NULL );
    
    GLenum err = glGetError();
    if (err) {
        std::cout << "Error setting up GL Texture...or previous GL command." << std::endl;
        Cleanup(true);
    }
}

void ClNddiDisplay::InitializeCl() {

	int err; // error code returned from api calls
	
	// Get an ID for the device
	err = clGetDeviceIDs(NULL, CL_DEVICE_TYPE_GPU, 1, &clDeviceId_, NULL);
	//err = clGetDeviceIDs(NULL, CL_DEVICE_TYPE_CPU, 1, &clDeviceId_, NULL);
	if (err != CL_SUCCESS) {
		std::cout << "Failed to get device IDs." << std::endl;
		Cleanup(true);
	}
    
    // Query the extensions supported
    size_t extensionsSize;
    err = clGetDeviceInfo(clDeviceId_, CL_DEVICE_EXTENSIONS, 0, NULL, &extensionsSize);
	if (err != CL_SUCCESS) {
		std::cout << "Failed to find out size of extensions string." << std::endl;
		Cleanup(true);
	}
    char* extensions = (char*)malloc(extensionsSize);
    err = clGetDeviceInfo(clDeviceId_, CL_DEVICE_EXTENSIONS, extensionsSize, extensions, &extensionsSize);
	if (err != CL_SUCCESS) {
		std::cout << "Failed to query the extensions." << std::endl;
		Cleanup(true);
	}
    std::cout << "Extensions: " << extensions << std::endl;
    free(extensions);
    
    CGLContextObj kCGLContext = CGLGetCurrentContext();
    CGLShareGroupObj kCGLShareGroup = CGLGetShareGroup(kCGLContext);
    cl_context_properties clContextProperties[] =
    {
        CL_CONTEXT_PROPERTY_USE_CGL_SHAREGROUP_APPLE,
        (cl_context_properties)kCGLShareGroup,
        0
    };

	// Create a context
	clContext_ = clCreateContext(clContextProperties, 1, &clDeviceId_, NULL, NULL, &err);
	if (!clContext_) {
		std::cout << "Failed to create context." << std::endl;
		Cleanup(true);
	}
	
	// Create a command queue
	clQueue_ = clCreateCommandQueue(clContext_, clDeviceId_, 0, &err);
	if (!clQueue_) {
		std::cout << "Failed to create command queue." << std::endl;
		Cleanup(true);
	}

	// Read program file
    const char* kernelFileName = "cl/computePixel.cl";
    std::ifstream kernelFile;
    kernelFile.open(kernelFileName, std::ifstream::in);
    // TODO(CDE): Figure this crap out! Why does getenv("PWD") fail when running release version in xcode?
    if (!kernelFile.is_open()) { kernelFile.open("/Users/cdestes/School/UNC/Research/pixelbridge/xcode/DerivedData/pixelbridge/Build/Products/Release/cl/computePixel.cl", std::ios::in); }
    if (!kernelFile.is_open())
    {
        std::cerr << "Failed to open program file: " << getenv("PWD") << kernelFileName << "." << std::endl;
        Cleanup(true);
    }
    std::ostringstream oss;
    oss << kernelFile.rdbuf();
    std::string srcStdStr = oss.str();
    const char *srcStr = srcStdStr.c_str();

	// Create the compute program from the source buffer
	clProgram_ = clCreateProgramWithSource(clContext_, 1, (const char **)&srcStr, NULL, &err);
	if (!clProgram_) {
		std::cout << "Failed to create program." << std::endl;
		Cleanup(true);
	}
	
	// Build the program executable
	err = clBuildProgram(clProgram_, 0, NULL, NULL, NULL, NULL);
	if (err != CL_SUCCESS) {
		size_t len;
		char buffer[2048];
		
		std::cout << "Error: Failed to build program executable\n" << std::endl;
		clGetProgramBuildInfo(clProgram_, clDeviceId_, CL_PROGRAM_BUILD_LOG,
							  sizeof(buffer), buffer, &len);
		printf("%s\n", buffer);
		Cleanup(true);
	}
	
	// Create the compute kernel in the program we wish to run
	clKernel_ = clCreateKernel(clProgram_, "computePixel", &err);
	if (!clKernel_ || err != CL_SUCCESS) {
		std::cout << "Failed to create compute kernel." << std::endl;
		Cleanup(true);
	}
	
    // Initialize the CL NDDI components
    cl_mem* inputVectorBuffer = clInputVector_->initializeCl(clContext_, clQueue_);
    cl_mem* coefficientPlaneBuffer = clCoefficientPlane_->initializeCl(clContext_, clQueue_);
    cl_mem* frameVolumeBuffer = clFrameVolume_->initializeCl(clContext_, clQueue_);
    if (!inputVectorBuffer || !coefficientPlaneBuffer || !frameVolumeBuffer) {
        std::cout << "Failed to do the CL initialization of the NDDI components." << std::endl;
		Cleanup(true);
    }
    
	// Create the frame volume dimensions buffer
    cl_uint fv_dims[CM_HEIGHT];
    for (int i = 0; i < CM_HEIGHT; i++) {
        fv_dims[i] = frameVolumeDimensionalSizes_[i];
    }
	clFrameVolumeDims_ = clCreateBuffer(clContext_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                        sizeof(cl_uint) * CM_HEIGHT, fv_dims, NULL);
    
    // Create the framebuffer output as a GL texture, so be directly used by the GPU to render the results
	clFrameBuffer_ = clCreateFromGLTexture2D(clContext_, CL_MEM_READ_WRITE, GL_TEXTURE_RECTANGLE_ARB, 0, texture_, NULL );
	if (!clFrameVolumeDims_ || !clFrameBuffer_) {
		std::cout << "Failed to create input/output memory arrays." << std::endl;
		Cleanup(true);
	}

	cl_uint cm_dims[2] = { CM_WIDTH, CM_HEIGHT };
    cl_uint display_dims[2] = { displayWidth_, displayHeight_ };

    // Set the arguments to our compute kernel
    err  = clSetKernelArg(clKernel_, 0, sizeof(cl_mem), inputVectorBuffer);
    err |= clSetKernelArg(clKernel_, 1, sizeof(cl_mem), coefficientPlaneBuffer);
    err |= clSetKernelArg(clKernel_, 2, sizeof(cl_mem), frameVolumeBuffer);
    err |= clSetKernelArg(clKernel_, 3, sizeof(cl_mem), &clFrameVolumeDims_);
    err |= clSetKernelArg(clKernel_, 4, sizeof(cl_mem), &clFrameBuffer_);
    err |= clSetKernelArg(clKernel_, 5, sizeof(cl_uint), &cm_dims[0]);
    err |= clSetKernelArg(clKernel_, 6, sizeof(cl_uint), &cm_dims[1]);
    err |= clSetKernelArg(clKernel_, 7, sizeof(cl_uint), &display_dims[0]);
    err |= clSetKernelArg(clKernel_, 8, sizeof(cl_uint), &display_dims[1]);
    if (err != CL_SUCCESS) {
        std::cout << "Failed to set kernel arguments." << std::endl;
        Cleanup(true);
    }
        
    // Initialize global and local workgroup size
    global[0] = displayWidth_; global[1] = displayHeight_;
	local[0] = 1; local[1] = 1;
    
    // Get the maximum work-group size for executing the kernel on the device
    err = clGetKernelWorkGroupInfo(clKernel_, clDeviceId_, CL_KERNEL_WORK_GROUP_SIZE,
                                   sizeof(size_t), &local, NULL);
    if (err != CL_SUCCESS) {
        std::cout << "Failed to get kernel workgroup info." << err << std::endl;
        Cleanup(true);
    }
    
    // Ignoring query above and hard coding workgroup size based on experimental data
    local[0]=128; local[1]=4;
    
    // Adjust global based on the new workgroup size in local
    global[0] = displayWidth_ / local[0] * local[0];
    if (global[0] < displayWidth_)
        global[0] += local[0];
    global[1] = displayHeight_ / local[1] * local[1];
    if (global[1] < displayHeight_)
        global[1] += local[1];
    cout << "Global Size: " << global[0] << " " << global[1] << " and Local Workgroup Size: " << local[0] << " " << local[1] << endl;
}


// Turn on for the timing data in the Render method below
#define OUTPUT_RENDER_TIMING_DATA

void ClNddiDisplay::Render() {
    
#ifdef OUTPUT_RENDER_TIMING_DATA
	timeval startTime, endTime; // Used for timing data
	gettimeofday(&startTime, NULL);
#endif
	
	int err;            // Holds return value of CL calls
	unsigned int count; // Number of pixels to compute
	
	// Initialize count
	count = displayWidth_ * displayHeight_;

    // Make sure GL's done and acquire the GL Object
    glFinish();
	err = clEnqueueAcquireGLObjects(clQueue_, 1, &clFrameBuffer_, 0, NULL, NULL);
	if (err) {
		std::cout << "Failed to enqueue acquire GL Objects command." << std::endl;
		Cleanup(true);
	}
	
	// Execute the kernel over the entire range of the data set
	err = clEnqueueNDRangeKernel(clQueue_, clKernel_, 2, NULL, global, local,
								 0, NULL, NULL);
	if (err) {
		std::cout << "Failed to enqueue ND range kernel command." << std::endl;
		Cleanup(true);
	}
	
	// Release the GL Object and wait for CL command queue to empty
	err = clEnqueueReleaseGLObjects(clQueue_, 1, &clFrameBuffer_, 0, NULL, NULL);
	if (err) {
		std::cout << "Failed to enqueue release GL Objects command." << std::endl;
		Cleanup(true);
	}
	clFinish(clQueue_);
	
#ifdef OUTPUT_RENDER_TIMING_DATA
	gettimeofday(&endTime, NULL);
	printf("Render Statistics:\n  Size: %dx%d - FPS: %f\n",
		   displayWidth_,
		   displayHeight_,
		   1.0f / ((double)(endTime.tv_sec * 1000000
							+ endTime.tv_usec
							- startTime.tv_sec * 1000000
							- startTime.tv_usec) / 1000000.0f)
		   );
#endif
}

void ClNddiDisplay::PutPixel(Pixel p, std::vector<unsigned int> location) {
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * (1 + frameVolumeDimensionalSizes_.size()));
    
    // Set the single pixel
	clFrameVolume_->PutPixel(p, location);
	
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void ClNddiDisplay::CopyPixelStrip(Pixel* p, std::vector<unsigned int> start, std::vector<unsigned int> end) {
	
	int dimensionToCopyAlong;
	bool dimensionFound = false;
	
    // Find the dimension to copy along
	for (int i = 0; !dimensionFound && (i < frameVolumeDimensionalSizes_.size()); i++) {
		if (start[i] != end[i]) {
			dimensionToCopyAlong = i;
			dimensionFound = true;
		}
	}
	
    // Register transmission cost now that we know the length of the strip sent
    costModel->registerTransmissionCharge(4 * ((end[dimensionToCopyAlong] - start[dimensionToCopyAlong] + 1) + 2 * frameVolumeDimensionalSizes_.size()));
    
    // Copy the pixels
    clFrameVolume_->CopyPixelStrip(p, start, end);
	
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void ClNddiDisplay::CopyPixels(Pixel* p, std::vector<unsigned int> start, std::vector<unsigned int> end) {
	
    // Register transmission cost first
    int pixelsToCopy = 1;
    for (int i = 0; i < start.size(); i++) {
        pixelsToCopy *= end[i] - start[i] + 1;
    }
    costModel->registerTransmissionCharge(4 * (pixelsToCopy + 2 * frameVolumeDimensionalSizes_.size()));
    
    // Copy pixels
    clFrameVolume_->CopyPixels(p, start, end);
    
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void ClNddiDisplay::FillPixel(Pixel p, std::vector<unsigned int> start, std::vector<unsigned int> end) {
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * (1 + 2 * frameVolumeDimensionalSizes_.size()));
    
    // Fill pixels
    clFrameVolume_->FillPixel(p, start, end);
    
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void ClNddiDisplay::CopyFrameVolume(std::vector<unsigned int> start, std::vector<unsigned int> end, std::vector<unsigned int> dest) {
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * (3 * frameVolumeDimensionalSizes_.size()));
    
    // Copy pixels
    clFrameVolume_->CopyFrameVolume(start, end, dest);
    
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void ClNddiDisplay::UpdateInputVector(std::vector<int> input) {
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * input.size());
	
    // Update the input vector
	clInputVector_->UpdateInputVector(input);
	
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void ClNddiDisplay::PutCoefficientMatrix(std::vector< std::vector<int> > coefficientMatrix,
                                           std::vector<unsigned int> location) {
	
    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * (CM_WIDTH * CM_SIZE + frameVolumeDimensionalSizes_.size()));
    
    // Update the coefficient matrix
    clCoefficientPlane_->PutCoefficientMatrix(coefficientMatrix, location);
	
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void ClNddiDisplay::FillCoefficientMatrix(std::vector< std::vector<int> > coefficientMatrix,
                                            std::vector<unsigned int> start,
                                            std::vector<unsigned int> end) {
    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * (CM_WIDTH * CM_SIZE + 2 * 2));
    
    // Fill the coefficient matrices
    clCoefficientPlane_->FillCoefficientMatrix(coefficientMatrix, start, end);
	
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}

void ClNddiDisplay::FillCoefficient(int coefficient,
                                      int row, int col,
                                      std::vector<unsigned int> start,
                                      std::vector<unsigned int> end) {
    // Register transmission cost first
    costModel->registerTransmissionCharge(4 * (3 + 2 * 2));
	
    // Fill the coefficient matrices
    clCoefficientPlane_->FillCoefficient(coefficient, row, col, start, end);
    
#ifndef SUPRESS_EXCESS_RENDERING
	Render();
#endif
}
