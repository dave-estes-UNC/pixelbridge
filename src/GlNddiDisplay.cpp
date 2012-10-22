#include <iostream>
#include <stdio.h>
#include <sys/time.h>

#include "PixelBridgeFeatures.h"
#include "GlNddiDisplay.h"

using namespace nddi;

// public

GlNddiDisplay::GlNddiDisplay(std::vector<unsigned int> frameVolumeDimensionalSizes,
                             int inputVectorSize) {
	GlNddiDisplay(frameVolumeDimensionalSizes, 320, 240, inputVectorSize);
}

GlNddiDisplay::GlNddiDisplay(std::vector<unsigned int> frameVolumeDimensionalSizes,
                             int displayWidth, int displayHeight,
                             int inputVectorSize) {

	frameVolumeDimensionalSizes_ = frameVolumeDimensionalSizes;
	displayWidth_ = displayWidth;
	displayHeight_ = displayHeight;

    // Create the CostModel
    costModel = new CostModel();

    // Setup Input Vector
	inputVector_ = new InputVector(costModel, inputVectorSize);

	// Setup framevolume and initialize to black
	frameVolume_ = new FrameVolume(costModel, frameVolumeDimensionalSizes);

	// Setup coefficient plane with zeroed coefficient matrices
	coefficientPlane_ = new CoefficientPlane(costModel, displayWidth_, displayHeight_, CM_WIDTH, CM_HEIGHT);

    // allocate a texture name
    glGenTextures( 1, &texture_ );

	// Setup framebuffer and initialize to black
	frameBuffer_ = (Pixel*)malloc(sizeof(Pixel) * displayWidth_ * displayHeight_);
	memset(frameBuffer_, 0x00, sizeof(Pixel) * displayWidth_ * displayHeight_);
}

// TODO(CDE): Why is the destructor for GlNddiDisplay being called when we're using a ClNddiDisplay?
GlNddiDisplay::~GlNddiDisplay() {

	delete(inputVector_);
	delete(frameVolume_);
	delete(coefficientPlane_);

    if (frameBuffer_)
        free(frameBuffer_);

    glDeleteTextures(1, &texture_);
}

// Private

// Turn on for the timing data in the Render method below
#define OUTPUT_RENDER_TIMING_DATA

void GlNddiDisplay::Render() {

#ifdef OUTPUT_RENDER_TIMING_DATA
	timeval startTime, endTime; // Used for timing data
	gettimeofday(&startTime, NULL);
#endif

    // Even though the InputVector can be invoked concurrently, it's really slow. So
    // we'll use a local copy instead and update the cost model in bulk later.
#ifndef NO_OMP
    int   * iv = inputVector_->data();
    Pixel * fv = frameVolume_->data();
#pragma omp parallel for
#endif // !NO_OMP
	for (int j = 0; j < displayHeight_; j++) {
		for (int i = 0; i < displayWidth_; i++) {
#ifndef NO_OMP
			frameBuffer_[j * displayWidth_ + i].packed = ComputePixel(i, j, iv, fv).packed;
#else
			frameBuffer_[j * displayWidth_ + i].packed = ComputePixel(i, j).packed;
#endif
		}
	}

    // Update the cost model for the in bulk now if we are using OpenMP since we bypassed the traditional
    // getters for input vector, frame volume, and coefficient matrix.
#ifndef NO_OMP
    costModel->registerBulkMemoryCharge(INPUT_VECTOR_COMPONENT,
                                        displayWidth_ * displayHeight_ * CM_HEIGHT * (CM_WIDTH - 2),
                                        READ_ACCESS,
                                        NULL,
                                        displayWidth_ * displayHeight_ * CM_HEIGHT * (CM_WIDTH - 2) * 4L,
                                        0);
    costModel->registerBulkMemoryCharge(COEFFICIENT_PLANE_COMPONENT,
                                        displayWidth_ * displayHeight_ * CM_HEIGHT * CM_WIDTH,
                                        READ_ACCESS,
                                        NULL,
                                        displayWidth_ * displayHeight_ * CM_HEIGHT * CM_WIDTH * 4L,
                                        0);
    costModel->registerBulkMemoryCharge(FRAME_VOLUME_COMPONENT,
                                        displayWidth_ * displayHeight_,
                                        READ_ACCESS,
                                        NULL,
                                        displayWidth_ * displayHeight_ * 4L,
                                        0);
    costModel->registerPixelMappingCharge(displayWidth_ * displayHeight_);
#endif

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

Pixel GlNddiDisplay::ComputePixel(unsigned int x, unsigned int y) {

	// Grab the coefficient matrix
    CoefficientMatrix * matrix = coefficientPlane_->getCoefficientMatrix(x, y);

	// Compute the position vector for the proper pixel in the frame volume.
	std::vector<unsigned int> fvPosition;
	// Matrix multiply the input vector by the coefficient matrix
	for (int j = 0; j < CM_HEIGHT; j++) {
        // Initialize to zero
		fvPosition.push_back(0);
        // No need to read the x and y from the input vector, just multiply directly.
        fvPosition[j] += matrix->getCoefficient(0, j) * x;
        fvPosition[j] += matrix->getCoefficient(1, j) * y;
        // Then multiply the remainder of the input vector
		for (int i = 2; i < CM_WIDTH; i++) {
			fvPosition[j] += matrix->getCoefficient(i, j) * inputVector_->getValue(i);
		}
	}

    costModel->registerPixelMappingCharge(1);

	return frameVolume_->getPixel(fvPosition);
}

Pixel GlNddiDisplay::ComputePixel(unsigned int x, unsigned int y, int* iv, Pixel* fv) {

    // Grab the coefficient matrix
    int * cm = coefficientPlane_->getCoefficientMatrix(x, y)->data();

	// Compute the position vector for the proper pixel in the frame volume.
	std::vector<unsigned int> fvPosition;
	// Matrix multiply the input vector by the coefficient matrix
	for (int j = 0; j < CM_HEIGHT; j++) {
        // Initialize to zero
		fvPosition.push_back(0);
        // No need to read the x and y from the input vector, just multiply directly.
        fvPosition[j] += cm[j * CM_WIDTH + 0] * x;
        fvPosition[j] += cm[j * CM_WIDTH + 1] * y;
        // Then multiply the remainder of the input vector
		for (int i = 2; i < CM_WIDTH; i++) {
			fvPosition[j] += cm[j * CM_WIDTH + i] * iv[i];
		}
	}

    // Compute the offset and grab the pixel directly from the frame volume
    unsigned int offset = 0;
    unsigned int multiplier = 1;

    for (int i = 0; i < fvPosition.size(); i++) {
        offset += fvPosition[i] * multiplier;
        multiplier *= frameVolumeDimensionalSizes_[i];
    }

    return fv[offset];
}

GLuint GlNddiDisplay::GetFrameBuffer() {

#ifdef SUPRESS_EXCESS_RENDERING
	Render();
#endif

// TODO(CDE): Temporarily putting this here until GlNddiDisplay and ClNddiDisplay
//            are using the exact same kind of GL textures
#ifdef NO_CL
    // select our current texture
    glBindTexture( GL_TEXTURE_2D, texture_ );

    // select modulate to mix texture with color for shading
    glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );

    // when texture area is small, bilinear filter the closest mipmap
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
					GL_LINEAR_MIPMAP_NEAREST );
    // when texture area is large, bilinear filter the first mipmap
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );

    // if wrap is true, the texture wraps over at the edges (repeat)
    //       ... false, the texture ends at the edges (clamp)
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,
					GL_CLAMP );
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T,
					GL_CLAMP );

    // build our texture mipmaps
    gluBuild2DMipmaps( GL_TEXTURE_2D, 3, displayWidth_, displayHeight_,
					  GL_RGBA, GL_UNSIGNED_BYTE, frameBuffer_ );
#endif

	return texture_;
}
