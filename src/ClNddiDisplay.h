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

/**
 * This is a version of the GlNddiDisplay that will use special OpenCL-aware
 * NDDI components and will do the rendering with OpenCL kernels.
 */
class ClNddiDisplay : public GlNddiDisplay, public NDimensionalDisplayInterfaceExtended {

public:
    ClNddiDisplay(std::vector<unsigned int> frameVolumeDimensionalSizes,
                  int inputVectorSize);
    ClNddiDisplay(std::vector<unsigned int> frameVolumeDimensionalSizes,
                  int displayWidth, int displayHeight,
                  int inputVectorSize);
    ~ClNddiDisplay();
    void PutPixel(Pixel p, std::vector<unsigned int> location);
    void CopyPixelStrip(Pixel* p, std::vector<unsigned int> start, std::vector<unsigned int> end);
    void CopyPixels(Pixel* p, std::vector<unsigned int> start, std::vector<unsigned int> end);
    void FillPixel(Pixel p, std::vector<unsigned int> start, std::vector<unsigned int> end);
    void CopyFrameVolume(std::vector<unsigned int> start, std::vector<unsigned int> end, std::vector<unsigned int> dest);
    void UpdateInputVector(std::vector<int> input);
    void PutCoefficientMatrix(std::vector< std::vector<int> > coefficientMatrix, std::vector<unsigned int> location);
    void FillCoefficientMatrix(std::vector< std::vector<int> > coefficientMatrix, std::vector<unsigned int> start, std::vector<unsigned int> end);
    void FillCoefficient(int coefficient, int row, int col, std::vector<unsigned int> start, std::vector<unsigned int> end);

    // To satisfy the NDimensionalDisplayInterfaceExtended interface
    void CopyFrameVolume(std::vector<unsigned int> start, std::vector<unsigned int> end, std::vector<unsigned int> dest, bool blend) {};
    void CopyPixelTiles(Pixel* p, std::vector<std::vector<unsigned int> > starts, std::vector<unsigned int> size);

private:
    void Render();

    void InitializeGl();
    void InitializeCl();

    void Cleanup(bool shouldExit);

    ClInputVector       * clInputVector_;
    ClCoefficientPlane  * clCoefficientPlane_;
    ClFrameVolume       * clFrameVolume_;

    cl_platform_id    clPlatformId_;
    cl_device_id      clDeviceId_;
    cl_context        clContext_;
    cl_command_queue  clQueue_;
    cl_program        clProgramComputePixel_;
    cl_kernel         clKernel_;
    cl_mem            clFrameVolumeDims_;
    cl_mem            clFrameBuffer_;
    cl_mem            clCommandPacket_;

    size_t global[2];
	size_t local[2];

	size_t            maxCommandPacketSize_;
};

#endif // CL_NDDI_DISPLAY_H
