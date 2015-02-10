/*
 * Configuration.h
 *
 *  Created on: Feb 17, 2014
 *      Author: cdestes
 */

#ifndef CONFIGURATION_H_
#define CONFIGURATION_H_

#include <vector>

using namespace std;

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

/**
 * Multiscale Support allows the DCT Tiler to use scaled 8x8 super-macroblocks
 * at multiple scales. The frame will be first approximated by the larger
 * scales and then refined slowly with the lower scales. The scale_config_t pairs
 * designate each of the scales supported, started with the top-most planes (plane 0).
 *
 * Note: The last plane is still reserved for medium gray, and therefore is not available
 *       as part of this configuration.
 */
struct scale_config_t {
    size_t scale_multiplier;
    size_t first_plane_idx;
    size_t plane_count;
};


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
    vector<scale_config_t> dctScales;

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

        scale_config_t simple = {1, 0, 63};
        dctScales.push_back(simple);
    }

    void clearDctScales() {
        dctScales.clear();
    }

    void addDctScale(size_t mult, size_t first, size_t count) {
        scale_config_t scale;

        scale.scale_multiplier = mult;
        scale.first_plane_idx = first;
        scale.plane_count = count;

        dctScales.push_back(scale);
    }

};

extern Configuration globalConfiguration;

#endif /* CONFIGURATION_H_ */
