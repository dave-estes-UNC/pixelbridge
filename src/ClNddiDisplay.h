#ifndef CL_NDDI_DISPLAY_H
#define CL_NDDI_DISPLAY_H


#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#include <CL/cl_gl.h>
#endif

#include "GlNddiDisplay.h"
#include "NDimensionalDisplayInterfaceExtended.h"
#include "ClInputVector.h"
#include "ClCoefficientPlane.h"
#include "ClFrameVolume.h"

// Must match struct in fillCoefficient.cl
typedef struct {
    int  coefficient;
    uint posCol;
    uint posRow;
    uint startX;
    uint startY;
    uint sizeW;
    uint sizeH;
} coefficient_update_t;

/**
 * This is a version of the GlNddiDisplay that will use special OpenCL-aware
 * NDDI components and will do the rendering with OpenCL kernels.
 */
class ClNddiDisplay : public GlNddiDisplay, public NDimensionalDisplayInterfaceExtended {

public:
    ClNddiDisplay(vector<unsigned int> &frameVolumeDimensionalSizes,
                  int inputVectorSize);
    ClNddiDisplay(vector<unsigned int> &frameVolumeDimensionalSizes,
                  int displayWidth, int displayHeight,
                  int inputVectorSize);
    ~ClNddiDisplay();
    void PutPixel(Pixel p, vector<unsigned int> &location);
    void CopyPixelStrip(Pixel* p, vector<unsigned int> &start, vector<unsigned int> &end);
    void CopyPixels(Pixel* p, vector<unsigned int> &start, vector<unsigned int> &end);
    void CopyPixelTiles(vector<Pixel*> &p, vector<vector<unsigned int> > &starts, vector<unsigned int> &size);
    void FillPixel(Pixel p, vector<unsigned int> &start, vector<unsigned int> &end);
    void CopyFrameVolume(vector<unsigned int> &start, vector<unsigned int> &end, vector<unsigned int> &dest);
    void UpdateInputVector(vector<int> &input);
    void PutCoefficientMatrix(vector< vector<int> > &coefficientMatrix, vector<unsigned int> &location);
    void FillCoefficientMatrix(vector< vector<int> > &coefficientMatrix, vector<unsigned int> &start, vector<unsigned int> &end);
    void FillCoefficient(int coefficient, int row, int col, vector<unsigned int> &start, vector<unsigned int> &end);
    void FillCoefficientTiles(vector<int> &coefficients, vector<vector<unsigned int> > &positions, vector<vector<unsigned int> > &starts, vector<unsigned int> &size);

    // To satisfy the NDimensionalDisplayInterfaceExtended interface
    void CopyFrameVolume(vector<unsigned int> &start, vector<unsigned int> &end, vector<unsigned int> &dest, bool blend) {};

private:
    void Render();

    void InitializeGl();
    void InitializeCl();
    void LoadKernel(char *path, char *file, char *name, cl_program *program, cl_kernel *kernel);

    void Cleanup(bool shouldExit);

    ClInputVector       *clInputVector_;
    ClCoefficientPlane  *clCoefficientPlane_;
    ClFrameVolume       *clFrameVolume_;

    cl_platform_id       clPlatformId_;
    cl_device_id         clDeviceId_;
    cl_context           clContext_;
    cl_command_queue     clQueue_;
    cl_program           clProgramComputePixel_;
    cl_kernel            clKernelComputePixel_;
    cl_program           clProgramFillCoefficient_;
    cl_kernel            clKernelFillCoefficient_;
    cl_mem               clFrameVolumeDims_;
    cl_mem               clFrameBuffer_;
    cl_mem               clCommandPacket_;

    size_t               globalComputePixel_[2];
	size_t               localComputePixel_[2];
    size_t               globalFillCoefficient_[2];
	size_t               localFillCoefficient_[2];

	size_t               maxCommandPacketSize_;
};

#endif // CL_NDDI_DISPLAY_H
