/*
 * Configuration.h
 *
 *  Created on: Feb 17, 2014
 *      Author: cdestes
 */

#ifndef CONFIGURATION_H_
#define CONFIGURATION_H_

typedef enum {
    COUNT,  // Just count (aka Perfect Pixel Latching)
    SIMPLE, // Simple Framebuffer
    FLAT,   // Tiled, but not cached
    CACHE,  // Tiled and cached
    DCT,    // 8x8 Macroblocks and DCT coefficients
    IT      // 4x4 Macroblocks and integer transform coefficients
} tiler_t;

typedef enum {
    NONE = 0,         // No blending
    FRAME_VOLUME,     // Uses frame volume blending operations
    TEMPORAL,         // Uses the input vector to switch rapidly between planes
    COEFFICIENT_PLANE // Uses multiple coefficient planes
} blend_t;

class Configuration {

public:

    tiler_t tiler;
    blend_t blend;
    size_t tileWidth;
    size_t tileHeight;
    size_t maxTiles;
    size_t sigBits;
    size_t quality;
    size_t startFrame;
    size_t maxFrames;
    size_t rewindStartFrame;
    size_t rewindFrames;
    size_t temporalFlipCountPerFrame;
    bool PSNR;
    bool verbose;
    bool headless;

public:

    Configuration() {
        tiler = CACHE;
        blend = NONE;
        tileWidth = 0;
        tileHeight = 0;
        maxTiles = 1000;
        sigBits = 8;
        quality = 4;
        startFrame = 0;
        maxFrames = 0;
        rewindStartFrame = 0;
        rewindFrames = 0;
        temporalFlipCountPerFrame = 4;
        PSNR = false;
        verbose = false;
        headless = false;
    }

};

extern Configuration globalConfiguration;

#endif /* CONFIGURATION_H_ */
