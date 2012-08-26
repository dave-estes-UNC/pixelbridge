#ifndef GL_NDDI_DISPLAY_H
#define GL_NDDI_DISPLAY_H

#include <GLUT/glut.h>

#include "BaseNddiDisplay.h"

/**
 * This class adds a public method that returns a reference to the frame buffer
 * so that a GLUT-based application can render it to the window.
 */
class GlNddiDisplay : public nddi::BaseNddiDisplay {
    
public:
    GlNddiDisplay() {}
    GlNddiDisplay(std::vector<unsigned int> frameVolumeDimensionalSizes,
                  int inputVectorSize);
    GlNddiDisplay(std::vector<unsigned int> frameVolumeDimensionalSizes,
                  int displayWidth, int displayHeight,
                  int inputVectorSize);
    ~GlNddiDisplay();
    
    /**
     * Used by GLUT-based application to get access to the frame buffer in order to display it.
     *
     * @return Texture holding a rendered frame.
     */
    GLuint GetFrameBuffer();
    
private:
    void Render();
    nddi::Pixel ComputePixel(unsigned int x, unsigned int y);
    nddi::Pixel ComputePixel(unsigned int x, unsigned int y, int* iv, nddi::Pixel* fv);

protected:
    GLuint texture_;
};

#endif // GL_NDDI_DISPLAY_H
