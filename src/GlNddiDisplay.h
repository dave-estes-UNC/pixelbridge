#ifndef GL_NDDI_DISPLAY_H
#define GL_NDDI_DISPLAY_H

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/freeglut.h>
#endif

#include "BaseNddiDisplay.h"

using namespace nddi;

/**
 * This class adds a public method that returns a reference to the frame buffer
 * so that a GLUT-based application can render it to the window.
 */
class GlNddiDisplay : public nddi::BaseNddiDisplay {

public:
    GlNddiDisplay() {}
    GlNddiDisplay(std::vector<unsigned int> &frameVolumeDimensionalSizes,
                  int numCoefficientPlanes, int inputVectorSize);
    GlNddiDisplay(std::vector<unsigned int> &frameVolumeDimensionalSizes,
                  int displayWidth, int displayHeight,
                  int numCoefficientPlanes, int inputVectorSize);
    ~GlNddiDisplay();

    /**
     * Used by GLUT-based application to get access to the frame buffer in order to display it.
     *
     * @return Texture holding a rendered frame.
     */
    GLuint GetFrameBufferTex();

    /**
     * Use to register in bulk all of the rendering cost. This will updates the counts correctly,
     * but any more detailed information is lost.
     */
    Pixel* GetFrameBuffer();
    
    /**
     * Triggers a simulted render that only records the cost estimated cost involved.
     */
    void SimulateRender();

    void SetPixelByteSignMode(SignMode mode);
    void SetFullScaler(uint16_t scaler);
    uint16_t GetFullScaler() { return fullScaler; }
    
private:
    void Render();
    nddi::Pixel ComputePixel(unsigned int x, unsigned int y);
#ifndef NO_OMP
    nddi::Pixel ComputePixel(unsigned int x, unsigned int y, int* iv, nddi::Pixel* fv);
#endif
    void RegisterBulkRenderCost();
    

protected:
    SignMode  pixelSignMode_;
    uint16_t  fullScaler;
    size_t    accumulatorShifter_;
    GLuint    texture_;
    Pixel    *frameBuffer_;
};

#endif // GL_NDDI_DISPLAY_H
