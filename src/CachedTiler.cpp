/*
 *  CachedTiler.cpp
 *  pixelbridge
 *
 *  Created by Dave Estes on 10/29/10.tile_map
 *  Copyright 2010 Dave Estes. All rights reserved.
 *
 */

#include <iostream>
#include <zlib.h>

#include "CachedTiler.h"


/**
 * The CachedTiler is created based on the dimensions of the NDDI display that's passed in. If those
 * dimensions change, then the CachedTiler should be destroyed and re-created.
 */
CachedTiler::CachedTiler (size_t display_width, size_t display_height,
                          size_t tile_width, size_t tile_height,
                          size_t max_tiles, size_t bits)
: tile_width_(tile_width),
  tile_height_(tile_height),
  max_tiles_(max_tiles),
  bits_(bits)
{
    quiet_ = globalConfiguration.headless || !globalConfiguration.verbose;

    // 3 dimensional matching the Tile Width x Height x max tiles
    vector<unsigned int> fvDimensions;
    fvDimensions.push_back(tile_width);
    fvDimensions.push_back(tile_height);
    fvDimensions.push_back(max_tiles);

#ifndef NO_CL
    display_ = new ClNddiDisplay(fvDimensions,                  // framevolume dimensional sizes
                                 display_width, display_height, // display size
                                 1,                             // number of coefficient planes in display
                                 3,                             // input vector size (x, y, and z)
                                 globalConfiguration.headless);
#else
    display_ = new GlNddiDisplay(fvDimensions,                  // framevolume dimensional sizes
                                 display_width, display_height, // display size
                                 1,                             // number of coefficient planes in display
                                 3,                             // input vector size (x, y, and z)
                                 globalConfiguration.headless);
#endif

    // Compute tile_map width
    tile_map_width_ = display_width / tile_width;
    if ((tile_map_width_ * tile_width) < display_width) { tile_map_width_++; }

    // Compute tile_map height
    tile_map_height_ = display_height / tile_height;
    if ((tile_map_height_ * tile_height) < display_height) { tile_map_height_++; }

    // Set up tile cache counters
    unchanged_tiles_ = cache_hits_ = cache_misses_ = age_counter_ = 0;

    // Set up the tile map, one column at a time
    display_map_.resize(tile_map_width_);
    for (int i = 0; i < tile_map_width_; i++) {
        for (int j = 0; j < tile_map_height_; j++) {
            display_map_[i].push_back(NULL);
        }
    }

    // Initialize Input Vector
    vector<int> iv;
    iv.push_back(1);
    display_->UpdateInputVector(iv);

    // Initialize Frame Volume
    nddi::Pixel p;
    p.r = p.g = p.b = p.a = 0xff;
    vector<unsigned int> start, end;
    start.push_back(0); start.push_back(0); start.push_back(0);
    end.push_back(tile_width_); end.push_back(tile_height_); end.push_back(max_tiles_);
    display_->FillPixel(p, start, end);

    // Initialize Coefficient Planes
    InitializeCoefficientPlanes();
}

CachedTiler::~CachedTiler()
{
    // Free all of the allocated tile_t structs
    map<unsigned long, tile_t*>::iterator it;
    for (it = cache_map_.begin(); it != cache_map_.end(); it++) {
        free(it->second);
    }
}

/**
 * Returns the Display created and initialized by the tiler.
 */
GlNddiDisplay* CachedTiler::GetDisplay() {
    return display_;
}

/**
 * Intializes the Coefficient Planes for this tiler.
 *
 * @return The cost of this operation, including all of the NDDI operations
 */
void CachedTiler::InitializeCoefficientPlanes() {

    // Setup the coefficient matrix to complete 3x3 identity initially
    vector< vector<int> > coeffs;
    coeffs.resize(3);
    coeffs[0].push_back(1); coeffs[0].push_back(0); coeffs[0].push_back(0);
    coeffs[1].push_back(0); coeffs[1].push_back(1); coeffs[1].push_back(0);
    coeffs[2].push_back(0); coeffs[2].push_back(0); coeffs[2].push_back(0);

    // Setup start and end points to (0,0) initially
    vector<unsigned int> start, end;
    start.push_back(0); start.push_back(0); start.push_back(0);
    end.push_back(0); end.push_back(0); end.push_back(0);

    for (int j = 0; j < tile_map_height_; j++) {
        for (int i = 0; i < tile_map_width_; i++) {
            coeffs[2][0] = -i * tile_width_;
            coeffs[2][1] = -j * tile_height_;
            start[0] = i * tile_width_; start[1] = j * tile_height_;
            end[0] = (i + 1) * tile_width_ - 1; end[1] = (j + 1) * tile_height_ - 1;
            if (end[0] >= display_->DisplayWidth()) { end[0] = display_->DisplayWidth() - 1; }
            if (end[1] >= display_->DisplayHeight()) { end[1] = display_->DisplayHeight() - 1; }
            display_->FillCoefficientMatrix(coeffs, start, end);
        }
    }

    // Turn off all planes and then set the 0 plane to full on.
    start[0] = 0; start[1] = 0; start[2] = 0;
    end[0] = display_->DisplayWidth() - 1;
    end[1] = display_->DisplayHeight() - 1;
    end[2] = display_->NumCoefficientPlanes() - 1;
    Scaler s;
    s.packed = 0;
    display_->FillScaler(s, start, end);
    end[2] = 0;
    s.r = s.g = s.b = display_->GetFullScaler();
    display_->FillScaler(s, start, end);
}


/**
 * Update the tile_map, tilecache, and then the NDDI display based on the frame that's passed in. The
 * frame is returned from the ffmpeg player as an RGB buffer. There is not Alpha channel.
 *
 * @param buffer Pointer to an RGB buffer
 * @param width The width of the RGB buffer
 * @param height The height of the RGB buffer
 */
void CachedTiler::UpdateDisplay(uint8_t* buffer, size_t width, size_t height)
{
    int                            unchanged = 0, hits = 0, misses = 0;
    unsigned char                  mask = 0xff << (8 - bits_);
    tile_t                        *tile;
    Pixel                         *tile_pixels = NULL;
    Pixel                         *tile_pixels_sig_bits = NULL;

    assert(width >= display_->DisplayWidth());
    assert(height >= display_->DisplayHeight());

    // Break up the passed in buffer into one tile at a time
    for (size_t j_tile_map = 0; j_tile_map < tile_map_height_; j_tile_map++) {
        for (size_t i_tile_map = 0; i_tile_map < tile_map_width_; i_tile_map++) {

            // Increment age counter
            age_counter_++;

            // Allocate tile pixel arrays if necessary. Sometimes they're re-used.
            if (!tile_pixels)
                tile_pixels = (Pixel*)malloc(tile_width_ * tile_height_ * sizeof(Pixel));
            if (!tile_pixels_sig_bits)
                tile_pixels_sig_bits = (Pixel*)malloc(tile_width_ * tile_height_ * sizeof(Pixel));

#ifndef NO_OMP
#pragma omp parallel for ordered
#endif
            for (size_t j_tile = 0; j_tile < tile_height_; j_tile++) {
#ifndef NO_OMP
#pragma omp ordered
#endif
                {
                    // Compute the offset into the RGB buffer for this row in this tile
                    int bufferOffset = 3 * ((j_tile_map * tile_height_ + j_tile) * width + (i_tile_map * tile_width_));

                    for (size_t i_tile = 0; i_tile < tile_width_; i_tile++) {

                        Pixel p, psb;

                        // Just use a black pixel if our tile is hanging off the edge of the buffer
                        if ((j_tile_map * tile_height_ + j_tile >= display_->DisplayHeight()) ||
                            (i_tile_map * tile_width_ + i_tile >= display_->DisplayWidth()) ) {
                            p.r = p.g = p.b = 0; p.a = 255;
                        } else {
                            p.r = buffer[bufferOffset++];
                            p.g = buffer[bufferOffset++];
                            p.b = buffer[bufferOffset++];
                            p.a = 0xff;
                        }
                        tile_pixels[j_tile * tile_width_ + i_tile].packed = p.packed;

                        psb.r = p.r & mask;
                        psb.g = p.g & mask;
                        psb.b = p.b & mask;
                        psb.a = p.a & mask;

                        tile_pixels_sig_bits[j_tile * tile_width_ + i_tile].packed = psb.packed;
                    }
                }
            }

            // Calculate the checksum
            unsigned long  checksum;
#if (CHECKSUM_CALCULATOR == TRIVIAL)
            checksum  = (unsigned long)tile_pixels_sig_bits[0].packed << 32;
            checksum |= (unsigned long)tile_pixels_sig_bits[tile_width_ * tile_height_ - 1].packed;
#else
            unsigned long crc = crc32(0L, Z_NULL, 0);
#if (CHECKSUM_CALCULATOR == CRC)
            checksum = crc32(crc, (unsigned char*)tile_pixels_sig_bits, tile_width_ * tile_height_ * sizeof(Pixel));
#elif (CHECKSUM_CALCULATOR == ADLER)
            checksum = adler32(crc, (unsigned char*)tile_pixels_sig_bits, tile_width_ * tile_height_ * sizeof(Pixel));
#endif
#endif

            // If the tile is already in the tile cache
            tile = IsTileInCache(checksum);
            if (tile) {
                // Remove the old age mapping
                age_map_.erase(tile->age);

                // Update the age and map it
                tile->age = age_counter_;
                age_map_.insert(pair<unsigned long, tile_t*>(tile->age, tile));

                // If the display doesn't already contain this tile
                if (display_map_[i_tile_map][j_tile_map] != tile) {
                    // Update the cache hit counter
                    hits++;
                    // Update the display map
                    display_map_[i_tile_map][j_tile_map] = tile;

#ifdef USE_COPY_PIXEL_TILES
                    // Push tile, updating coefficient plane only
                    PushTile(tile, i_tile_map, j_tile_map);
#else
                    // Update Coefficient Matrices
                    UpdateCoefficientMatrices(i_tile_map, j_tile_map, tile);
#endif
                } else {
                    // Update the unchanged counter
                    unchanged++;
                }

            // The tile is not already in the cache, so we'll need to add it
            } else {
                // Cache miss
                misses++;

                // If we have room in the tile cache
                if (cache_map_.size() < max_tiles_) {

                    // Allocate a new tile
                    tile = (tile_t*)malloc(sizeof(tile_t));
                    tile->checksum = checksum;

                    // This will be pushed at the end with this new "index". Since we're growing the cache still,
                    // we can use the cache index of this new element as the zIndex
                    tile->zIndex = cache_map_.size();

                    // Update the display map with the new tile
                    display_map_[i_tile_map][j_tile_map] = tile;

                    // Map the checksum to this tile
                    cache_map_.insert(pair<unsigned long, tile_t*>(tile->checksum, tile));

                    // Update the age and map it
                    tile->age = age_counter_;
                    age_map_.insert(pair<unsigned long, tile_t*>(tile->age, tile));

#ifdef USE_COPY_PIXEL_TILES
                    // Push tile
                    PushTile(tile, tile_pixels, i_tile_map, j_tile_map);

                    // Force new tile to be allocated. This one will be freed after it's copied
                    tile_pixels = NULL;
#else
                    // Push the pixels to the frame volume.
                    UpdateFrameVolume(tile_pixels, tile);

                    // Update Coefficient Matrices
                    UpdateCoefficientMatrices(i_tile_map, j_tile_map, tile);
#endif

                // We didn't have room in the tile cache, so update an existing tile
                } else {

                    bool needToUpdateCoefficients = false;

                    // Determine if we can re-use the zIndex from the previous tile, saving coefficient matrix updates
                    tile = display_map_[i_tile_map][j_tile_map];
                    if (IsTileInUse(tile)) {
                        // It was in use, so try to find an expired one.
                        tile = GetExpiredCacheTile();

                        // If we couldn't, then we're in trouble
                        assert(tile);

                        needToUpdateCoefficients = true;

                        // Update the tile map with the new tile
                        display_map_[i_tile_map][j_tile_map] = tile;
                    }

                    // Remove the old checksum and age mappings
                    cache_map_.erase(tile->checksum);
                    age_map_.erase(tile->age);

                    // Set the new checksum and age
                    tile->checksum = checksum;
                    tile->age = age_counter_;

                    // Add the new checksum and age mappings
                    cache_map_.insert(pair<unsigned long, tile_t*>(tile->checksum, tile));
                    age_map_.insert(pair<unsigned long, tile_t*>(tile->age, tile));

#ifdef USE_COPY_PIXEL_TILES
                    if (!needToUpdateCoefficients) {
                        // Push tile, updating FrameVolume only
                        PushTile(tile, tile_pixels);
                    } else {
                        // Push tile, FrameVolume and CoeffientPlane
                        PushTile(tile, tile_pixels, i_tile_map, j_tile_map);
                    }

                    // Force new tile to be allocated. This one will be freed after it's copied
                    tile_pixels = NULL;
#else
                    // Push the pixels to the frame volume.
                    UpdateFrameVolume(tile_pixels, tile);

                    // Update Coefficient Matrices if needed
                    if (needToUpdateCoefficients) {
                        UpdateCoefficientMatrices(i_tile_map, j_tile_map, tile);
                    }
#endif
                }
            } // for (int i_tile_map = 0; i_tile_map < tile_map_width_; i_tile_map++) {
        } // for (int j_tile_map = 0; j_tile_map < tile_map_height_; j_tile_map++) {
    }

#ifdef USE_COPY_PIXEL_TILES
    // If any tiles were updated
    if (misses > 0) {

        // Set the tile size parameter
        vector<unsigned int> size;
        size.push_back(tile_width_); size.push_back(tile_height_);

        // Update the Frame Volume by copying the tiles over
        if (tile_pixels_list.size() > 0) {
            display_->CopyPixelTiles(tile_pixels_list, tile_starts_list, size);
        }

        // Free the tile memory and empty the vector
        while (!tile_pixels_list.empty()) {
            free(tile_pixels_list.back());
            tile_pixels_list.pop_back();
        }
        tile_starts_list.clear();

        // Update the coefficient plane
        if (coefficients_list.size() > 0) {
            display_->FillCoefficientTiles(coefficients_list,
                                           coefficient_positions_list,
                                           coefficient_plane_starts_list,
                                           size);
        }

        // Free the coefficient memory and empty the vector
        coefficients_list.clear();
        coefficient_positions_list.clear();
        coefficient_plane_starts_list.clear();
    }
#endif

    // Free alloc'd memory
    if (tile_pixels)
        free(tile_pixels);
    if (tile_pixels_sig_bits)
        free(tile_pixels_sig_bits);

    // Report cache statistics
    unchanged_tiles_ += unchanged;
    cache_hits_ += hits;
    cache_misses_ += misses;
    if (!quiet_) {
        cout << "Cached Tiling Statistics:" << endl << "  unchanged: " << unchanged_tiles_ << " cache hits: " << cache_hits_ << " cache misses: " << cache_misses_ << " cache size: " << cache_map_.size() << endl;
    }
}

/**
 * Checks to see if the tile is already in the tile cache based solely on the
 * checksum.
 *
 * @param checksum The checksum of the tile to search for.
 * @return The pointer to the tile_t.
 */
tile_t* CachedTiler::IsTileInCache(unsigned long checksum) {

    map<unsigned long, tile_t*>::iterator it;

    it = cache_map_.find(checksum);
    if (it != cache_map_.end()) {
        return it->second;
    } else {
        return NULL;
    }
}

/**
 * Checks to see if the tile is being used in the tile map based solely on the
 * checksum.
 *
 * @param tile The tile whose age is used to determine if it's in use.
 * @return true if still in use
 * @return false if not in use
 */
bool CachedTiler::IsTileInUse(tile_t* tile) {

    // If the tile's age is greater than the current counter minus the total number
    // of tiles on the display, then it's likely in use.
    return (tile->age >= (age_counter_ - tile_map_width_ * tile_map_height_));
}

/**
 * Finds a candidate to be ejected from the cache.
 *
 * @return The tile to be ejected
 */
tile_t* CachedTiler::GetExpiredCacheTile() {
    pair<map<unsigned long, tile_t*>::iterator, map<unsigned long, tile_t*>::iterator> its;

    its = age_map_.equal_range(0);

    return its.first->second;
}

/**
 * Updates region of the Frame Volume corresponding to the tile's
 * zIndex.
 *
 * @param pixels The array of pixels to use in the update.
 * @param tile The zIndex of the tile will be used to choose the region in the frame volume.
 * @return The cost of the NDDI operations.
 */
#ifndef USE_COPY_PIXEL_TILES
void CachedTiler::UpdateFrameVolume(Pixel* pixels, tile_t *tile) {

    // Setup start and end points
    vector<unsigned int> start, end;
    start.push_back(0); start.push_back(0); start.push_back(tile->zIndex);
    end.push_back(tile_width_ - 1); end.push_back(tile_height_ - 1); end.push_back(tile->zIndex);

    display_->CopyPixels(pixels, start, end);
}


/**
 * Updates the coefficient matrices for the specified tile map location with
 * the tile provided.
 *
 * @param x X coordinate of the tile in the tile map.
 * @param y Y coordinate of the tile in the tile map.
 * @param tile The zIndex from this tile will be used to update the coefficient matrices.
 * @return The cost of the NDDI operations.
 */
void CachedTiler::UpdateCoefficientMatrices(size_t x, size_t y, tile_t *tile) {

    // Setup start and end points
    vector<unsigned int> start, end;
    start.push_back(x * tile_width_); start.push_back(y * tile_height_); start.push_back(0);
    end.push_back((x + 1) * tile_width_ - 1); end.push_back((y + 1) * tile_height_ - 1); end.push_back(0);
    if (end[0] >= display_->DisplayWidth()) { end[0] = display_->DisplayWidth() - 1; }
    if (end[1] >= display_->DisplayHeight()) { end[1] = display_->DisplayHeight() - 1; }

    display_->FillCoefficient(tile->zIndex, 2, 2, start, end);
}

#else

/**
 * Pushes the pixels for a tile and the coefficient updates into their respective lists.
 * The lists will later be sent in two bulk "packets" to the NDDI display.
 *
 * @param tile The zIndex of the tile will be used to choose the region in the frame volume.
 * @param pixels The array of pixels to use in the update.
 * @param i The new X coordinate of the tile in the tile map.
 * @param i The new Y coordinate of the tile in the tile map.
 */
void CachedTiler::PushTile(tile_t* tile, Pixel* pixels, size_t i, size_t j) {
    PushTile(tile, pixels);
    PushTile(tile, i, j);
}

/**
 * Pushes the pixels for a tile and the coefficient updates into their respective lists.
 * The lists will later be sent in two bulk "packets" to the NDDI display.
 *
 * @param tile The zIndex of the tile will be used to choose the region in the frame volume.
 * @param pixels The array of pixels to use in the update.
 */
void CachedTiler::PushTile(tile_t* tile, Pixel* pixels) {
    // Create and push the start coordinates
    vector<unsigned int> start;
    start.push_back(0); start.push_back(0); start.push_back(tile->zIndex);

#ifndef NO_OMP
#pragma omp critical
#endif
    {
        // Push the tile and starts
        tile_pixels_list.push_back(pixels);
        tile_starts_list.push_back(start);
    }
}

/**
 * Pushes the pixels for a tile and the coefficient updates into their respective lists.
 * The lists will later be sent in two bulk "packets" to the NDDI display.
 *
 * @param tile The zIndex of the tile will be used to choose the region in the frame volume.
 * @param i The X coordinate of the tile in the tile map.
 * @param j The Y coordinate of the tile in the tile map.
 */
void CachedTiler::PushTile(tile_t* tile, size_t i, size_t j) {
    // Create and push the position and start coordinates
    vector<unsigned int> position;
    position.push_back(2); position.push_back(2);

    vector<unsigned int> start;
    start.push_back(i * tile_width_); start.push_back(j * tile_height_); start.push_back(0);

#ifndef NO_OMP
#pragma omp critical
#endif
    {
        // Push the tile and starts
        coefficients_list.push_back(tile->zIndex);
        coefficient_positions_list.push_back(position);
        coefficient_plane_starts_list.push_back(start);
    }
}

#endif
