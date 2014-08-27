#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/time.h>

#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#else
#include <GL/glx.h>
// glx.h includes X11/Xlib.h which defines these macros.
// Get rid of them so we can compile.
#undef DisplayWidth
#undef DisplayHeight
#endif

#include "PixelBridgeFeatures.h"
#include "ClNddiDisplay.h"

using namespace nddi;

// public

ClNddiDisplay::ClNddiDisplay(vector<unsigned int> &frameVolumeDimensionalSizes,
                             int numCoefficientPlanes, int inputVectorSize)
:   clFrameVolume_(NULL),
    clInputVector_(NULL),
    clKernelFillCoefficient_(0),
    maxCommandPacketSize_(0)
{
    ClNddiDisplay(frameVolumeDimensionalSizes, 320, 240, numCoefficientPlanes, inputVectorSize);
}

ClNddiDisplay::ClNddiDisplay(vector<unsigned int> &frameVolumeDimensionalSizes,
                             int displayWidth, int displayHeight,
                             int numCoefficientPlanes, int inputVectorSize)
{
    numPlanes_ = numCoefficientPlanes;
    frameVolumeDimensionalSizes_ = frameVolumeDimensionalSizes;
    displayWidth_ = displayWidth;
    displayHeight_ = displayHeight;
    quiet_ = true;

    // Create the CostModel
    costModel = new CostModel();

    // Setup Input Vector
    clInputVector_ = new ClInputVector(costModel, inputVectorSize);
    // TODO(CDE): NULL this out
    inputVector_ = (InputVector*)clInputVector_;

    // Setup coefficient plane with zeroed coefficient matrices
    clCoefficientPlane_ = new ClCoefficientPlanes(costModel, displayWidth_, displayHeight_, numCoefficientPlanes, CM_WIDTH, CM_HEIGHT);
    // TODO(CDE): NULL this out
    coefficientPlanes_ = (ClCoefficientPlanes*)clCoefficientPlane_;

    // Setup framevolume and initialize to black
    clFrameVolume_ = new ClFrameVolume(costModel, frameVolumeDimensionalSizes);
    // TODO(CDE): NULL this out
    frameVolume_ = (FrameVolume*)clFrameVolume_;

    // We won't be using the base components or frameBuffer_, so make them null.
    //inputVector_ = NULL;
    //coefficientPlanes_ = NULL;
    //frameVolume_ = NULL;
    frameBuffer_ = NULL;

    // Set the maximum command packet size based on a 3D frame volume and enough 8x8 tiles to fill the display
    // TODO(CDE): Revisit this maximum when I start using the command packet to update the coefficient plane
    maxCommandPacketSize_ = (4 * 2) + (displayWidth_ / 8 + 1) * (displayHeight / 8 + 1) * (4 * (8 * 8 + 3));

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

    if (clKernelComputePixel_ != 0)
        clReleaseKernel(clKernelComputePixel_);

    if (clProgramComputePixel_ != 0)
        clReleaseProgram(clProgramComputePixel_);

    if (clContext_ != 0)
        clReleaseContext(clContext_);

    if (clFrameVolumeDims_ != 0)
        clReleaseMemObject(clFrameVolumeDims_);

    if (clFrameBuffer_ != 0)
        clReleaseMemObject(clFrameBuffer_);

    if (clCommandPacket_ != 0)
        clReleaseMemObject(clCommandPacket_);

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

    texture_ = 0;

    // Allocate texture as a CL-friendly image
    glGenTextures( 1, &texture_ );
    glTexEnvi( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );
    glBindTexture(GL_TEXTURE_RECTANGLE_ARB, texture_);
    glTexImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, GL_RGBA32F_ARB, displayWidth_,
                 displayHeight_, 0, GL_LUMINANCE, GL_FLOAT, NULL );

    GLenum err = glGetError();
    if (err) {
        cout << "Error setting up GL Texture...or previous GL command." << endl;
        Cleanup(true);
    }
}

void ClNddiDisplay::InitializeCl() {

    cl_uint  numPlatforms;
    cl_int   err;

    // Get platform ID
    err = clGetPlatformIDs(1, &clPlatformId_, &numPlatforms);

    // Get an ID for the device
    err = clGetDeviceIDs(clPlatformId_, CL_DEVICE_TYPE_GPU, 1, &clDeviceId_, NULL);
    //err = clGetDeviceIDs(NULL, CL_DEVICE_TYPE_CPU, 1, &clDeviceId_, NULL);
    if (err != CL_SUCCESS) {
        cout << "Failed to get device IDs." << endl;
        Cleanup(true);
    }

    // Query the extensions supported
    size_t extensionsSize;
    err = clGetDeviceInfo(clDeviceId_, CL_DEVICE_EXTENSIONS, 0, NULL, &extensionsSize);
    if (err != CL_SUCCESS) {
        cout << "Failed to find out size of extensions string." << endl;
        Cleanup(true);
    }
    char* extensions = (char*)malloc(extensionsSize);
    err = clGetDeviceInfo(clDeviceId_, CL_DEVICE_EXTENSIONS, extensionsSize, extensions, &extensionsSize);
    if (err != CL_SUCCESS) {
        cout << "Failed to query the extensions." << endl;
        Cleanup(true);
    }
    cout << "Extensions: " << extensions << endl;
    free(extensions);

    // Create the context
    cl_context_properties clContextProperties[] =
    {
    CL_CONTEXT_PLATFORM, (cl_context_properties)clPlatformId_,
#ifdef __APPLE__
#ifndef NO_CL
#error OpenCL Not supported on Mac OS X Perfectly yet.
#endif
#else
    CL_GL_CONTEXT_KHR, (cl_context_properties)glXGetCurrentContext(),
    CL_GLX_DISPLAY_KHR, (cl_context_properties)glXGetCurrentDisplay(),
#endif
        0
    };

    // Create a context
    clContext_ = clCreateContext(clContextProperties, 1, &clDeviceId_, NULL, NULL, &err);
    if (!clContext_) {
        cout << "Failed to create context." << endl;
        Cleanup(true);
    }

    // Create a command queue
#ifdef CL_PROFILING_ENABLED
    clQueue_ = clCreateCommandQueue(clContext_, clDeviceId_, CL_QUEUE_PROFILING_ENABLE, &err);
#else
    clQueue_ = clCreateCommandQueue(clContext_, clDeviceId_, 0, &err);
#endif
    if (!clQueue_) {
        cout << "Failed to create command queue." << endl;
        Cleanup(true);
    }

    // Load kernels
    // TODO(CDE): Figure this crap out! Why does getenv("PWD") fail when running release version in xcode?
    char path[] = "/home/cdestes/Work/pixelbridge/src";

    char computePixelFileName[] = "computePixel.cl";
    char computePixelName[] = "computePixel";
    LoadKernel(path, computePixelFileName, computePixelName, &clProgramComputePixel_, &clKernelComputePixel_);

    char fillCoefficientFileName[] = "fillCoefficient.cl";
    char fillCoefficientName[] = "fillCoefficient";
    LoadKernel(path, fillCoefficientFileName, fillCoefficientName, &clProgramFillCoefficient_, &clKernelFillCoefficient_);

    // Initialize the CL NDDI components
    cl_mem inputVectorBuffer = clInputVector_->initializeCl(clContext_, clQueue_);
    cl_mem coefficientPlaneBuffer = clCoefficientPlane_->initializeCl(clContext_, clQueue_);
    cl_mem frameVolumeBuffer = clFrameVolume_->initializeCl(clContext_, clQueue_);
    if (!inputVectorBuffer || !coefficientPlaneBuffer || !frameVolumeBuffer) {
        cout << "Failed to do the CL initialization of the NDDI components." << endl;
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
    clFrameBuffer_ = clCreateFromGLTexture(clContext_, CL_MEM_READ_WRITE, GL_TEXTURE_RECTANGLE_ARB, 0, texture_, NULL );
    if (!clFrameVolumeDims_ || !clFrameBuffer_) {
        cout << "Failed to create input/output memory arrays." << endl;
        Cleanup(true);
    }

    // Create the command packet buffer
    clCommandPacket_ = clCreateBuffer(clContext_, CL_MEM_READ_ONLY,
                                      maxCommandPacketSize_, NULL, NULL);
    clFrameVolume_->setCommandPacket(clCommandPacket_, maxCommandPacketSize_);

    cl_uint cm_dims[2] = { CM_WIDTH, CM_HEIGHT };
    cl_uint display_dims[2] = { displayWidth_, displayHeight_ };

    // Set the arguments to our computePixel kernel
    err  = clSetKernelArg(clKernelComputePixel_, 0, sizeof(cl_mem), &inputVectorBuffer);
    err |= clSetKernelArg(clKernelComputePixel_, 1, sizeof(cl_mem), &coefficientPlaneBuffer);
    err |= clSetKernelArg(clKernelComputePixel_, 2, sizeof(cl_mem), &frameVolumeBuffer);
    err |= clSetKernelArg(clKernelComputePixel_, 3, sizeof(cl_mem), &clFrameVolumeDims_);
    err |= clSetKernelArg(clKernelComputePixel_, 4, sizeof(cl_mem), &clFrameBuffer_);
    err |= clSetKernelArg(clKernelComputePixel_, 5, sizeof(cl_uint), &cm_dims[0]);
    err |= clSetKernelArg(clKernelComputePixel_, 6, sizeof(cl_uint), &cm_dims[1]);
    err |= clSetKernelArg(clKernelComputePixel_, 7, sizeof(cl_uint), &display_dims[0]);
    err |= clSetKernelArg(clKernelComputePixel_, 8, sizeof(cl_uint), &display_dims[1]);
    if (err != CL_SUCCESS) {
        cout << "Failed to set computePixel kernel arguments." << endl;
        Cleanup(true);
    }

    // Initialize global and local workgroup size for computePixel kernel
    globalComputePixel_[0] = displayWidth_; globalComputePixel_[1] = displayHeight_;
    localComputePixel_[0] = 1; localComputePixel_[1] = 1;

    // Get the maximum work-group size for executing the kernel on the device
    err = clGetKernelWorkGroupInfo(clKernelComputePixel_, clDeviceId_, CL_KERNEL_WORK_GROUP_SIZE,
                                   sizeof(size_t), &localComputePixel_, NULL);
    if (err != CL_SUCCESS) {
        cout << "Failed to get kernel workgroup info." << err << endl;
        Cleanup(true);
    }

    // Adjust global based on the new workgroup size in local
    globalComputePixel_[0] = displayWidth_ / localComputePixel_[0] * localComputePixel_[0];
    if (globalComputePixel_[0] < displayWidth_)
        globalComputePixel_[0] += localComputePixel_[0];
    globalComputePixel_[1] = displayHeight_ / localComputePixel_[1] * localComputePixel_[1];
    if (globalComputePixel_[1] < displayHeight_)
        globalComputePixel_[1] += localComputePixel_[1];
    cout << "computePixel Global Size: " << globalComputePixel_[0] << " " << globalComputePixel_[1] << " and Local Workgroup Size: " << localComputePixel_[0] << " " << localComputePixel_[1] << endl;

    // Set the arguments to our computePixel kernel
    err  = clSetKernelArg(clKernelFillCoefficient_, 0, sizeof(cl_mem), &clCommandPacket_);
    err |= clSetKernelArg(clKernelFillCoefficient_, 1, sizeof(cl_mem), &coefficientPlaneBuffer);
    err |= clSetKernelArg(clKernelFillCoefficient_, 2, sizeof(cl_uint), &cm_dims[0]);
    err |= clSetKernelArg(clKernelFillCoefficient_, 3, sizeof(cl_uint), &cm_dims[1]);
    err |= clSetKernelArg(clKernelFillCoefficient_, 4, sizeof(cl_uint), &display_dims[0]);
    err |= clSetKernelArg(clKernelFillCoefficient_, 5, sizeof(cl_uint), &display_dims[1]);
    if (err != CL_SUCCESS) {
        cout << "Failed to set fillCoefficient kernel arguments." << endl;
        Cleanup(true);
    }

    // Initialize local workgroup size for fillCoefficient kernel
    // Note: The global size is recalculated every time before enqueing the kernel
    localFillCoefficient_[0] = 1; localFillCoefficient_[1] = 1;

    // Get the maximum work-group size for executing the kernel on the device
    err = clGetKernelWorkGroupInfo(clKernelFillCoefficient_, clDeviceId_, CL_KERNEL_WORK_GROUP_SIZE,
                                   sizeof(size_t), &localFillCoefficient_, NULL);
    if (err != CL_SUCCESS) {
        cout << "Failed to get kernel workgroup info." << err << endl;
        Cleanup(true);
    }
    cout << "fillCoefficient Local Workgroup Size: " << localFillCoefficient_[0] << " " << localFillCoefficient_[1] << endl;
}

void ClNddiDisplay::LoadKernel(char *path, char *file, char *name, cl_program *program, cl_kernel *kernel) {

    cl_int   err;
    char filename[256];

    // Read feature header file
    sprintf(filename, "%s/%s", path, "PixelBridgeFeatures.h");
    ifstream headerFile;
    headerFile.open(filename, ifstream::in);
    if (!headerFile.is_open()) {
        headerFile.open(filename, ios::in);
    }
    if (!headerFile.is_open())
    {
        cerr << "Failed to open PixelBridgeFeatures.h." << endl;
        Cleanup(true);
    }

    // Read program file
    sprintf(filename, "%s/%s", path, file);
    ifstream kernelFile;
    kernelFile.open(filename, ifstream::in);
    if (!kernelFile.is_open()) {
        kernelFile.open(filename, ios::in);
    }
    if (!kernelFile.is_open())
    {
        cerr << "Failed to open program file: " << filename << "." << endl;
        Cleanup(true);
    }

    // Create the combined stream for both the header and kernel files.
    ostringstream oss;
    oss << headerFile.rdbuf();
    oss << kernelFile.rdbuf();
    string srcStdStr = oss.str();
    const char *srcStr = srcStdStr.c_str();

    // Create the compute program from the source buffer
    *program = clCreateProgramWithSource(clContext_, 1, (const char **)&srcStr, NULL, &err);
    if (!*program) {
        cout << "Failed to create program." << endl;
        Cleanup(true);
    }

    // Build the program executable
    err = clBuildProgram(*program, 0, NULL, NULL, NULL, NULL);
    if (err != CL_SUCCESS) {
        size_t len;
        char buffer[2048];

        cout << "Error: Failed to build program executable\n" << endl;
        clGetProgramBuildInfo(*program, clDeviceId_, CL_PROGRAM_BUILD_LOG,
                              sizeof(buffer), buffer, &len);
        printf("%s\n", buffer);
        Cleanup(true);
    }

    // Create the compute kernel in the program we wish to run
    *kernel = clCreateKernel(*program, name, &err);
    if (!*kernel || err != CL_SUCCESS) {
        cout << "Failed to create kernel." << endl;
        Cleanup(true);
    }
}


void ClNddiDisplay::Render() {

    timeval startTime, endTime; // Used for timing data
    if (!quiet_)
        gettimeofday(&startTime, NULL);

    int err;            // Holds return value of CL calls
    unsigned int count; // Number of pixels to compute

    // Initialize count
    count = displayWidth_ * displayHeight_;

    // Make sure GL's done and acquire the GL Object
    glFinish();
    err = clEnqueueAcquireGLObjects(clQueue_, 1, &clFrameBuffer_, 0, NULL, NULL);
    if (err) {
        cout << "Failed to enqueue acquire GL Objects command " << err << endl;
        Cleanup(true);
    }

    // Execute the kernel over the entire range of the data set
    err = clEnqueueNDRangeKernel(clQueue_, clKernelComputePixel_, 2, NULL, globalComputePixel_, localComputePixel_,
                                 0, NULL, NULL);
    if (err) {
        cout << "Failed to enqueue ND range kernel command " << err << endl;
        Cleanup(true);
    }

    // Release the GL Object and wait for CL command queue to empty
    err = clEnqueueReleaseGLObjects(clQueue_, 1, &clFrameBuffer_, 0, NULL, NULL);
    if (err) {
        cout << "Failed to enqueue release GL Objects command " << err << endl;
        Cleanup(true);
    }
    clFinish(clQueue_);

    if (!quiet_) {
        gettimeofday(&endTime, NULL);
        printf("Render Statistics:\n  Size: %dx%d - FPS: %f\n",
               displayWidth_,
               displayHeight_,
               1.0f / ((double)(endTime.tv_sec * 1000000
                                + endTime.tv_usec
                                - startTime.tv_sec * 1000000
                                - startTime.tv_usec) / 1000000.0f)
               );
    }
}

void ClNddiDisplay::PutPixel(Pixel p, vector<unsigned int> &location) {

    // Register transmission cost first
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_PIXELS(1) +         // One Pixel
                                          CALC_BYTES_FOR_FV_COORD_TUPLES(1), // One Coordinate Tuple
                                          0);

    // Set the single pixel
    clFrameVolume_->PutPixel(p, location);

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void ClNddiDisplay::CopyPixelStrip(Pixel* p, vector<unsigned int> &start, vector<unsigned int> &end) {

    int dimensionToCopyAlong;
    bool dimensionFound = false;

    // Find the dimension to copy along
    for (int i = 0; !dimensionFound && (i < frameVolumeDimensionalSizes_.size()); i++) {
        if (start[i] != end[i]) {
            dimensionToCopyAlong = i;
            dimensionFound = true;
        }
    }
    int pixelsToCopy = end[dimensionToCopyAlong] - start[dimensionToCopyAlong] + 1;

    // Register transmission cost now that we know the length of the strip sent
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_PIXELS(pixelsToCopy) +    // A strip of pixels
                                          CALC_BYTES_FOR_FV_COORD_TUPLES(2),       // Two Coordinate Tuples
                                          0);

    // Copy the pixels
    clFrameVolume_->CopyPixelStrip(p, start, end);

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void ClNddiDisplay::CopyPixels(Pixel* p, vector<unsigned int> &start, vector<unsigned int> &end) {

    // Register transmission cost first
    int pixelsToCopy = 1;
    for (int i = 0; i < start.size(); i++) {
        pixelsToCopy *= end[i] - start[i] + 1;
    }
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_PIXELS(pixelsToCopy) +    // Range of pixels
                                          CALC_BYTES_FOR_FV_COORD_TUPLES(2),       // Two Coordinate Tuples
                                          0);

    // Copy pixels
    clFrameVolume_->CopyPixels(p, start, end);

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void ClNddiDisplay::CopyPixelTiles(vector<Pixel*> &p, vector<vector<unsigned int> > &starts, vector<unsigned int> &size) {

    size_t tile_count = p.size();

    // Ensure parameter vectors' sizes match
    assert(starts.size() == tile_count);
    assert(starts[0].size() == frameVolumeDimensionalSizes_.size());
    assert(size.size() == 2);

    // Register transmission cost first
    int pixelsToCopy = 1;
    for (int i = 0; i < size.size(); i++) {
        pixelsToCopy *= size[i];
    }
    pixelsToCopy *= starts.size();
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_PIXELS(pixelsToCopy) +            // t tiles of x by y tiles of pixels
                                          CALC_BYTES_FOR_FV_COORD_TUPLES(tile_count + 1) + // t start coordinate tuples + 1 tuple for tile size dimensions
                                          CALC_BYTES_FOR_TILE_COORD_DOUBLES(1),            // 1 X by Y tile dimension double
                                          0);

    // Copy pixels (copies to host array, sets up packet and sends it to the device
    clFrameVolume_->CopyPixelTiles(p, starts, size);

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void ClNddiDisplay::FillPixel(Pixel p, vector<unsigned int> &start, vector<unsigned int> &end) {

    // Register transmission cost first
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_PIXELS(1) +         // One Pixel
                                          CALC_BYTES_FOR_FV_COORD_TUPLES(2), // Two Coordinate Tuples
                                          0);

    // Fill pixels
    clFrameVolume_->FillPixel(p, start, end);

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void ClNddiDisplay::CopyFrameVolume(vector<unsigned int> &start, vector<unsigned int> &end, vector<unsigned int> &dest) {

    // Register transmission cost first
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_FV_COORD_TUPLES(3), // Three Coordinate Tuples
                                          0);

    // Copy pixels
    clFrameVolume_->CopyFrameVolume(start, end, dest);

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void ClNddiDisplay::UpdateInputVector(vector<int> &input) {

    assert(input.size() == inputVector_->getSize() - 2);

    // Register transmission cost first
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_IV_UPDATE(), // Input Vector
                                          0);

    // Update the input vector
    clInputVector_->UpdateInputVector(input);

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void ClNddiDisplay::PutCoefficientMatrix(vector< vector<int> > &coefficientMatrix,
                                           vector<unsigned int> &location) {

    // Register transmission cost first
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_CMS(1) +             // One coefficient matrix
                                          CALC_BYTES_FOR_CP_COORD_TRIPLES(1), // One Coefficient Plane Coordinate triple
                                          0);

    // Update the coefficient matrix
    clCoefficientPlane_->PutCoefficientMatrix(coefficientMatrix, location);

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void ClNddiDisplay::FillCoefficientMatrix(vector< vector<int> > &coefficientMatrix,
                                            vector<unsigned int> &start,
                                            vector<unsigned int> &end) {
    // Register transmission cost first
    costModel->registerTransmissionCharge(CALC_BYTES_FOR_CMS(1) +             // One coefficient matrix
                                          CALC_BYTES_FOR_CP_COORD_TRIPLES(2), // Two Coefficient Plane Coordinate triples
                                          0);

    // Fill the coefficient matrices
    clCoefficientPlane_->FillCoefficientMatrix(coefficientMatrix, start, end);

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void ClNddiDisplay::FillCoefficient(int coefficient,
                                      int row, int col,
                                      vector<unsigned int> &start,
                                      vector<unsigned int> &end) {
    assert(row >= 0 && row < CM_HEIGHT);
    assert(col >= 0 && col < CM_WIDTH);
    assert(start.size() == 3);
    assert(end.size() == 3);

    // Register transmission cost first
    costModel->registerTransmissionCharge(BYTES_PER_COEFF * 1 +                // One coefficient
                                          CALC_BYTES_FOR_CM_COORD_DOUBLES(1) + // One Coefficient Matrix Coordinate double
                                          CALC_BYTES_FOR_CP_COORD_TRIPLES(2),  // Two Coefficient Plane Coordinate triples
                                          0);

    // Fill the coefficient matrices
    clCoefficientPlane_->FillCoefficient(coefficient, row, col, start, end);

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

// TODO(CDE): Consider moving this to ClCoefficientPlane even though it executes a kernel
void ClNddiDisplay::FillCoefficientTiles(vector<int> &coefficients,
                                         vector<vector<unsigned int> > &positions,
                                         vector<vector<unsigned int> > &starts,
                                         vector<unsigned int> &size) {

    cl_int   err;
    static coefficient_update_t *packet = NULL;

    // Set the number of instances for kernel
    size_t tile_count = coefficients.size();
    assert(positions.size() == tile_count);
    assert(starts.size() == tile_count);

    // Register transmission cost first
    costModel->registerTransmissionCharge(BYTES_PER_COEFF * tile_count +                 // t coefficients
                                          CALC_BYTES_FOR_CM_COORD_DOUBLES(tile_count) +  // t Coefficient Matrix Coordinate doubles
                                          CALC_BYTES_FOR_CP_COORD_TRIPLES(tile_count) +  // t Coefficient Plane Coordinate triples
                                          CALC_BYTES_FOR_TILE_COORD_DOUBLES(1),          // 1 X by Y tile dimension double
                                          0);

    // Build the packet
    if (packet) free(packet);
    packet = (coefficient_update_t*)malloc(sizeof(coefficient_update_t) * tile_count);
    for (int i = 0; i < tile_count; i++) {
        packet[i].coefficient = coefficients[i];
        packet[i].posCol = positions[i][0];
        packet[i].posRow = positions[i][1];
        packet[i].startX = starts[i][0];
        packet[i].startY = starts[i][1];
        packet[i].sizeW = size[0];
        packet[i].sizeH = size[1];
    }

    // Enqueue command to write the packet
    err = clEnqueueWriteBuffer(clQueue_, clCommandPacket_, CL_FALSE,
                               0, sizeof(coefficient_update_t) * tile_count, packet,
                               0, NULL, NULL); // TODO(CDE): Set the event param when I move it to ClCoefficientPlane
    if (err != CL_SUCCESS) {
        cout << __FUNCTION__ << " - Failed to create enqueue write buffer command." << err << endl;
    }

    // Set the num kernel arg
    err = clSetKernelArg(clKernelFillCoefficient_, 6, sizeof(cl_uint), &tile_count);
    if (err != CL_SUCCESS) {
        cout << "Failed to set fillCoefficient kernel argument." << endl;
        Cleanup(true);
    }

    // Adjust global based on the new workgroup size in local
    globalFillCoefficient_[0] = tile_count / localFillCoefficient_[0] * localFillCoefficient_[0];
    if (globalFillCoefficient_[0] < tile_count)
        globalFillCoefficient_[0] += localFillCoefficient_[0];
    globalFillCoefficient_[1] = localFillCoefficient_[1];
    //cout << "fillCoefficient Global Size: " << globalFillCoefficient_[0] << " " << globalFillCoefficient_[1] << " and Local Workgroup Size: " << localFillCoefficient_[0] << " " << localFillCoefficient_[1] << endl;

    // Enqueue kernel
    err = clEnqueueNDRangeKernel(clQueue_, clKernelFillCoefficient_, 1, NULL, globalFillCoefficient_, localFillCoefficient_,
                                 0, NULL, NULL);
    if (err) {
        cout << __FUNCTION__ << " - Failed to enqueue ND range kernel command " << err << endl;
        Cleanup(true);
    }

#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = true;
#else
    Render();
#endif
}

void ClNddiDisplay::FillScaler(Scaler scaler,
                               vector<unsigned int> &start,
                               vector<unsigned int> &end) {
    // TODO(CDE): Implement #MultiPlaneCL
}


void ClNddiDisplay::FillScalerTiles(vector<uint64_t> &scalers,
                                    vector<vector<unsigned int> > &starts,
                                    vector<unsigned int> &size) {
    // TODO(CDE): Implement #MultiPlaneCL
}
void ClNddiDisplay::FillScalerTileStack(vector<uint64_t> &scalers,
                                        vector<unsigned int> &start,
                                        vector<unsigned int> &size) {
    // TODO(CDE): Implement #MultiPlaneCL
}
