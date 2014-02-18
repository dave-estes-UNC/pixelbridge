#include <iostream>
#include <stdio.h>
#include <sys/time.h>
#include <cmath>

#include "PixelBridgeFeatures.h"
#include "GlNddiDisplay.h"

inline uint8_t CLAMP_SIGNED_BYTE(int32_t i) {
    uint32_t ret;

    if (i < 0) {
        ret = 0;
    } else if (i > 0xff) {
        ret = 0xff;
    } else {
        ret = i;
    }

    return (uint8_t)ret;
}

inline uint8_t CLAMP_UNSIGNED_BYTE(uint32_t i) {
    uint32_t ret;

    if (i > 0xff) {
        ret = 0xff;
    } else {
        ret = i;
    }

    return ret;
}

inline uint8_t TRUNCATE_BYTE(int32_t i) {
    uint32_t ret;

    ret = i & 0xff;

    return (uint8_t)ret;
}

// public

GlNddiDisplay::GlNddiDisplay(vector<unsigned int> &frameVolumeDimensionalSizes,
                             int numCoefficientPlanes, int inputVectorSize) {
    texture_ = 0;
    GlNddiDisplay(frameVolumeDimensionalSizes, 320, 240, numCoefficientPlanes, inputVectorSize);
}

GlNddiDisplay::GlNddiDisplay(vector<unsigned int> &frameVolumeDimensionalSizes,
                             int displayWidth, int displayHeight,
                             int numCoefficientPlanes, int inputVectorSize) {

    numPlanes_ = numCoefficientPlanes;
    frameVolumeDimensionalSizes_ = frameVolumeDimensionalSizes;
    displayWidth_ = displayWidth;
    displayHeight_ = displayHeight;
    pixelSignMode_ = UNSIGNED_MODE;
    quiet_ = true;

    // Create the CostModel
    costModel = new CostModel();

    // Setup Input Vector
    inputVector_ = new InputVector(costModel, inputVectorSize);

    // Setup framevolume and initialize to black
    frameVolume_ = new FrameVolume(costModel, frameVolumeDimensionalSizes);

    // Setup coefficient plane with zeroed coefficient matrices
    coefficientPlane_ = new CoefficientPlane(costModel, displayWidth_, displayHeight_, numCoefficientPlanes, CM_WIDTH, CM_HEIGHT);

    // allocate a texture name
    glGenTextures( 1, &texture_ );

    // Setup framebuffer and initialize to black
    frameBuffer_ = (Pixel*)malloc(sizeof(Pixel) * displayWidth_ * displayHeight_);
    memset(frameBuffer_, 0x00, sizeof(Pixel) * displayWidth_ * displayHeight_);

    // Set the full scaler and the accumulator
    SetFullScaler(DEFAULT_FULL_SCALER);
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

void GlNddiDisplay::Render() {

    timeval startTime, endTime; // Used for timing data
    if (!quiet_)
        gettimeofday(&startTime, NULL);

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
    RegisterBulkRenderCost();
#endif

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
#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = false;
#endif

}

Pixel GlNddiDisplay::ComputePixel(unsigned int x, unsigned int y) {

    int32_t    rAccumulator = 0, gAccumulator = 0, bAccumulator = 0;
    Pixel      q;

    q.packed = 0;

    // Accumulate color channels for the pixels chosen by each plane
    for (unsigned int p = 0; p < numPlanes_; p++) {

        // Grab the scaler for this location
        Scaler scaler = coefficientPlane_->getScaler(x, y, p);
#ifdef SKIP_COMPUTE_WHEN_SCALER_ZERO
        if (scaler.packed == 0) continue;
#endif

        // Grab the coefficient matrix
        CoefficientMatrix * matrix = coefficientPlane_->getCoefficientMatrix(x, y, p);
        assert(matrix);

        // Compute the position vector for the proper pixel in the frame volume.
        vector<unsigned int> fvPosition;
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
        q = frameVolume_->getPixel(fvPosition);
        if (pixelSignMode_ == UNSIGNED_MODE) {
            rAccumulator += (uint8_t)q.r * scaler.r;
            gAccumulator += (uint8_t)q.g * scaler.g;
            bAccumulator += (uint8_t)q.b * scaler.b;
        } else {
            rAccumulator += (int8_t)q.r * scaler.r;
            gAccumulator += (int8_t)q.g * scaler.g;
            bAccumulator += (int8_t)q.b * scaler.b;
        }
    }

    // Note: This shift operation will be absolutely necessary when this is implemented
    //       in hardware to avoid the division operation.
    if (pixelSignMode_ == UNSIGNED_MODE) {
        q.r = CLAMP_UNSIGNED_BYTE(rAccumulator >> accumulatorShifter_);
        q.g = CLAMP_UNSIGNED_BYTE(gAccumulator >> accumulatorShifter_);
        q.b = CLAMP_UNSIGNED_BYTE(bAccumulator >> accumulatorShifter_);
    } else {
        q.r = CLAMP_SIGNED_BYTE(rAccumulator >> accumulatorShifter_);
        q.g = CLAMP_SIGNED_BYTE(gAccumulator >> accumulatorShifter_);
        q.b = CLAMP_SIGNED_BYTE(bAccumulator >> accumulatorShifter_);
    }
    q.a = 255;

    costModel->registerPixelMappingCharge(1);

    return q;
}

#ifndef NO_OMP
// Duplicate version of the proper ComputePixel which avoids locking
Pixel GlNddiDisplay::ComputePixel(unsigned int x, unsigned int y, int* iv, Pixel* fv) {

    int32_t    rAccumulator = 0, gAccumulator = 0, bAccumulator = 0;
    Pixel      q;

    q.packed = 0;

    // Accumulate color channels for the pixels chosen by each plane
    for (unsigned int p = 0; p < numPlanes_; p++) {

        // Grab the scaler for this location.
#ifdef NARROW_DATA_STORES
        int16_t * scalerData = coefficientPlane_->data();
#else
        int * scalerData = coefficientPlane_->data();
#endif
        Scaler scaler;
        scaler.r = scalerData[(p * displayWidth_ * displayHeight_ + y * displayWidth_ + x) * 3 + 0];
        scaler.g = scalerData[(p * displayWidth_ * displayHeight_ + y * displayWidth_ + x) * 3 + 1];
        scaler.b = scalerData[(p * displayWidth_ * displayHeight_ + y * displayWidth_ + x) * 3 + 2];
#ifdef SKIP_COMPUTE_WHEN_SCALER_ZERO
        if (scaler.packed == 0) continue;
#endif

        // Grab the coefficient matrix
#ifdef NARROW_DATA_STORES
        int16_t * cm = coefficientPlane_->getCoefficientMatrix(x, y, p)->data();
#else
        int * cm = coefficientPlane_->getCoefficientMatrix(x, y, p)->data();
#endif

        // Compute the position vector for the proper pixel in the frame volume.
        vector<unsigned int> fvPosition;
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

        q = fv[offset];
        if (pixelSignMode_ == UNSIGNED_MODE) {
            rAccumulator += (uint8_t)q.r * scaler.r;
            gAccumulator += (uint8_t)q.g * scaler.g;
            bAccumulator += (uint8_t)q.b * scaler.b;
        } else {
            rAccumulator += (int8_t)q.r * scaler.r;
            gAccumulator += (int8_t)q.g * scaler.g;
            bAccumulator += (int8_t)q.b * scaler.b;
        }
    }

    // Note: This shift operation will be absolutely necessary when this is implemented
    //       in hardware to avoid the division operation.
    if (pixelSignMode_ == UNSIGNED_MODE) {
        q.r = CLAMP_UNSIGNED_BYTE(rAccumulator >> accumulatorShifter_);
        q.g = CLAMP_UNSIGNED_BYTE(gAccumulator >> accumulatorShifter_);
        q.b = CLAMP_UNSIGNED_BYTE(bAccumulator >> accumulatorShifter_);
    } else {
        q.r = CLAMP_SIGNED_BYTE(rAccumulator >> accumulatorShifter_);
        q.g = CLAMP_SIGNED_BYTE(gAccumulator >> accumulatorShifter_);
        q.b = CLAMP_SIGNED_BYTE(bAccumulator >> accumulatorShifter_);
    }
    q.a = 255;

    return q;
}
#endif

GLuint GlNddiDisplay::GetFrameBufferTex() {

#ifdef SUPRESS_EXCESS_RENDERING
    if (changed_)
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

void GlNddiDisplay::RegisterBulkRenderCost() {

    // For each pixel computed with all of the planes, the input vector (except x,y) is read
    costModel->registerBulkMemoryCharge(INPUT_VECTOR_COMPONENT,
                                        displayWidth_ * displayHeight_ * numPlanes_ * CM_HEIGHT * (CM_WIDTH - 2),
                                        READ_ACCESS,
                                        NULL,
                                        displayWidth_ * displayHeight_ * numPlanes_ * CM_HEIGHT * (CM_WIDTH - 2) * BYTES_PER_IV_VALUE,
                                        0);
    // For each pixel computed with all of the planes, the coefficient and scaler is read
    costModel->registerBulkMemoryCharge(COEFFICIENT_PLANE_COMPONENT,
                                        displayWidth_ * displayHeight_ * numPlanes_ * (CM_HEIGHT * CM_WIDTH) +
                                        displayWidth_ * displayHeight_ * numPlanes_ * (1),
                                        READ_ACCESS,
                                        NULL,
                                        displayWidth_ * displayHeight_ * numPlanes_ * (CM_HEIGHT * CM_WIDTH) * BYTES_PER_COEFF +
                                        displayWidth_ * displayHeight_ * numPlanes_ * (1) * BYTES_PER_SCALER,
                                        0);
    // For each pixel computed with all of the planes, a pixel sample is pulled from the frame volume
    costModel->registerBulkMemoryCharge(FRAME_VOLUME_COMPONENT,
                                        displayWidth_ * displayHeight_ * numPlanes_,
                                        READ_ACCESS,
                                        NULL,
                                        displayWidth_ * displayHeight_ * numPlanes_ * BYTES_PER_PIXEL,
                                        0);
    costModel->registerPixelMappingCharge(displayWidth_ * displayHeight_);
}

void GlNddiDisplay::SimulateRender() {
#ifdef SUPRESS_EXCESS_RENDERING
    if (changed_) {
        RegisterBulkRenderCost();
        changed_ = false;
    }
#endif
}

Pixel* GlNddiDisplay::GetFrameBuffer() {

#ifdef SUPRESS_EXCESS_RENDERING
    if (changed_)
        Render();
#endif

    return frameBuffer_;
}

void GlNddiDisplay::SetPixelByteSignMode(SignMode mode) {
    assert(mode == UNSIGNED_MODE || mode == SIGNED_MODE);
    pixelSignMode_ = mode;
}

void GlNddiDisplay::SetFullScaler(uint16_t scaler) {
    if (scaler & (scaler- 1)) {
        cout << "ERROR: THE FULL_SCALER specified is not a power of two." << endl;
    } else {
        // Set the full scaler
        fullScaler = scaler;

        // Initiliaze the shifter used during accumulation
        double s = log2((double)fullScaler);
        accumulatorShifter_ = int(s);
    }
}
