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

#include "PixelBridgeFeatures.h"
#include "CachedTiler.h"


/**
 * The CachedTiler is created based on the dimensions of the NDDI display that's passed in. If those
 * dimensions change, then the CachedTiler should be destroyed and re-created.
 */
CachedTiler::CachedTiler (GlNddiDisplay* display, size_t tile_width, size_t tile_height, size_t max_tiles, size_t bits, bool quiet)
: display_(display),
  tile_width_(tile_width),
  tile_height_(tile_height),
  max_tiles_(max_tiles),
  bits_(bits),
  quiet_(quiet)
{

	// Compute tile_map width
	tile_map_width_ = display_->DisplayWidth() / tile_width;
	if ((tile_map_width_ * tile_width) < display_->DisplayWidth()) { tile_map_width_++; }

	// Compute tile_map height
	tile_map_height_ = display_->DisplayHeight() / tile_height;
	if ((tile_map_height_ * tile_height) < display_->DisplayHeight()) { tile_map_height_++; }

	// Set up tile cache and counters
	tile_cache_.reserve(max_tiles_);
	unchanged_tiles_ = cache_hits_ = cache_misses_ = 0;

	// Blank tile for map
	tile_t tile;
	tile.checksum = 0;
	tile.zIndex = 0;

	// Set up the tile map, one column at a time
	tile_map_.resize(tile_map_width_);
	for (int i = 0; i < tile_map_width_; i++) {
		for (int j = 0; j < tile_map_height_; j++) {
			tile_map_[i].push_back(tile);
		}
	}
}

/**
 * Intializes the Coefficient Plane for this tiler.
 *
 * @return The cost of this operation, including all of the NDDI operations
 */
void CachedTiler::InitializeCoefficientPlane() {

	// Setup the coefficient matrix to complete 3x3 identity initially
	vector< vector<int> > coeffs;
    coeffs.resize(3);
    coeffs[0].push_back(1); coeffs[0].push_back(0); coeffs[0].push_back(0);
    coeffs[1].push_back(0); coeffs[1].push_back(1); coeffs[1].push_back(0);
    coeffs[2].push_back(0); coeffs[2].push_back(0); coeffs[2].push_back(0);

	// Setup start and end points to (0,0) initially
	vector<unsigned int> start, end;
	start.push_back(0); start.push_back(0);
	end.push_back(0); end.push_back(0);

	for (int j = 0; j < tile_map_height_; j++) {
		for (int i = 0; i < tile_map_width_; i++) {
			coeffs[2][0] = -i * tile_width_;
			coeffs[2][1] = -j * tile_height_;
			start[0] = i * tile_width_; start[1] = j * tile_height_;
			end[0] = (i + 1) * tile_width_ - 1; end[1] = (j + 1) * tile_height_ - 1;
			if (end[0] >= display_->DisplayWidth()) { end[0] = display_->DisplayWidth() - 1; }
			if (end[1] >= display_->DisplayHeight()) { end[1] = display_->DisplayHeight() - 1; }
			((ClNddiDisplay *)display_)->FillCoefficientMatrix(coeffs, start, end);
		}
	}
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
    Pixel                         *tile_pixels = NULL;
    Pixel                         *tile_pixels_sig_bits = NULL;

	// Break up the passed in buffer into one tile at a time
	for (int j_tile_map = 0; j_tile_map < tile_map_height_; j_tile_map++) {
		for (int i_tile_map = 0; i_tile_map < tile_map_width_; i_tile_map++) {
			// Initialize a tile
			tile_t tile;
			tile.checksum = 0L;
			tile.zIndex = 0;

			// Allocate tiles is necessary. Sometimes they're re-used.
            if (!tile_pixels)
            	tile_pixels = (Pixel*)malloc(tile_width_ * tile_height_ * sizeof(Pixel));
            if (!tile_pixels_sig_bits)
            	tile_pixels_sig_bits = (Pixel*)malloc(tile_width_ * tile_height_ * sizeof(Pixel));

#ifndef NO_OMP
#pragma omp parallel for ordered
#endif
			for (int j_tile = 0; j_tile < tile_height_; j_tile++) {
#ifndef NO_OMP
#pragma omp ordered
#endif
                {
                    // Compute the offset into the RGB buffer for this row in this tile
                    int bufferOffset = 3 * ((j_tile_map * tile_height_ + j_tile) * width + (i_tile_map * tile_width_));

                    for (int i_tile = 0; i_tile < tile_width_; i_tile++) {

                        Pixel p, psb;

                        // Just use a black pixel if our tile is hanging off the edge of the buffer
                        if ((j_tile_map * tile_height_ + j_tile >= height) ||
                            (i_tile_map * tile_width_ + i_tile >= width) ) {
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

			unsigned long crc = crc32(0L, Z_NULL, 0);
			tile.checksum = crc32(crc, (unsigned char*)tile_pixels_sig_bits, tile_width_ * tile_height_ * sizeof(Pixel));
			//tile.checksum = adler32(crc, (unsigned char*)tile_pixels_sig_bits, tile_width_ * tile_height_ * sizeof(Pixel));

			// If the tile is already in the tile cache
			int index = IsTileInCache(tile);
			if (index >= 0) {
				// Move it to the head of the tile cache
				tile = tile_cache_[index];
				tile_cache_.erase(tile_cache_.begin() + index);
				tile_cache_.insert(tile_cache_.begin(), tile);

				// If the tilemap doesn't already contain this tile
				if (tile_map_[i_tile_map][j_tile_map].checksum != tile.checksum) {
					// Update the cache hit counter
					hits++;
					// Update the tilemap
					tile_map_[i_tile_map][j_tile_map] = tile;

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
				if (tile_cache_.size() < max_tiles_) {
					// This will be pushed at the end with this new index
					size_t index = tile_cache_.size();

					// Since we're growing the cache still, we can use the cache index of this new element as the zIndex
					tile.zIndex = index;

					// Push the tile onto the head of the cache
					tile_cache_.insert(tile_cache_.begin(), tile);

					// Update the tile map with the new tile
					tile_map_[i_tile_map][j_tile_map] = tile;

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

				// We didn't have room in the tile cache
				} else {

					// Determine if we can re-use the zIndex from the previous tile, saving coefficient matrix updates
					tile_t previous_tile = tile_map_[i_tile_map][j_tile_map];
					tile_map_[i_tile_map][j_tile_map].checksum = 0;

					int index = IsTileInCache(previous_tile);

					// If the previous tile at this location is not being used anymore but it is still in the tile cache
					if (!IsTileInMap(previous_tile) && index != -1) {
						// Set the zIndex
						tile.zIndex = previous_tile.zIndex;

						// Update the tile cache be removing the previous one and inserting the new one at the front
						tile_cache_.erase(tile_cache_.begin() + index);
						tile_cache_.insert(tile_cache_.begin(), tile);

						// Update the tile map with the new tile. The zIndex is the same, so no need to update coefficient matrices
						tile_map_[i_tile_map][j_tile_map] = tile;

#ifdef USE_COPY_PIXEL_TILES
						// Push tile, updating FrameVolume only
						PushTile(tile, tile_pixels);

						// Force new tile to be allocated. This one will be freed after it's copied
						tile_pixels = NULL;
#else
						// Push the pixels to the frame volume.
						UpdateFrameVolume(tile_pixels, tile);
#endif

					// We couldn't re-use the previous zIndex, so we'll have to find another candidate in the cache
					} else {

						// Get the least recently used tile in the cache that is not currently used
						int index = GetExpiredCacheTile();

						if (index >= 0) {
							// Set the zIndex
							tile.zIndex = tile_cache_[index].zIndex;

							// Update the tile cache be removing the previous one and inserting the new one at the front
							tile_cache_.erase(tile_cache_.begin() + index);
							tile_cache_.insert(tile_cache_.begin(), tile);

							// Update the tile map with the new tile
							tile_map_[i_tile_map][j_tile_map] = tile;

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

						} else {
							// ERROR
							cout << "ERROR: Couldn't find a candidate to eject from the cache." << endl;
						}
					}
				}
			}
		}
	}

#ifdef USE_COPY_PIXEL_TILES
	// If any tiles were updated
	if (misses > 0) {

		// Set the tile size parameter
		vector<unsigned int> size;
		size.push_back(tile_width_); size.push_back(tile_height_);

		// Update the Frame Volume by copying the tiles over
		if (tile_pixels_list.size() > 0) {
			((ClNddiDisplay *)display_)->CopyPixelTiles(tile_pixels_list, tile_starts_list, size);
		}

		// Free the tile memory and empty the vector
		while (!tile_pixels_list.empty()) {
			free(tile_pixels_list.back());
			tile_pixels_list.pop_back();
		}
		tile_starts_list.clear();

		// Update the coefficient plane
		if (coefficients_list.size() > 0) {
			((ClNddiDisplay *)display_)->FillCoefficientTiles(coefficients_list,
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
		cout << "Cached Tiling Statistics:" << endl << "  unchanged: " << unchanged_tiles_ << " cache hits: " << cache_hits_ << " cache misses: " << cache_misses_ << " cache size: " << tile_cache_.size() << endl;
	}
}

/**
 * Checks to see if the tile is already in the tile cache based solely on the
 * checksum.
 *
 * @param tile The checksum in this tile is searched for in the tile cache based on a checksum match.
 * @return The index of the tile in the cache if it is found.
 * @return -1 if no match is found.
 */
int CachedTiler::IsTileInCache(tile_t tile) {

	for (int i = 0; i < tile_cache_.size(); i++) {
		if (tile_cache_[i].checksum == tile.checksum) {
			return i;
		}
	}

	return -1;
}

/**
 * Checks to see if the tile is being used in the tile map based solely on the
 * checksum.
 *
 * @param tile The checksum in this tile is searched for in the tile map based on a checksum match.
 * @return true if found at least once
 * @return false if not found
 */
bool CachedTiler::IsTileInMap(tile_t tile) {

	for (int i = 0; i < tile_map_width_; i++) {
		for (int j = 0; j < tile_map_height_; j++) {
			if (tile_map_[i][j].checksum == tile.checksum) {
				return true;
			}
		}
	}

	return false;
}

/**
 * Finds a candidate to be ejected from the cache.
 *
 * @return The index of the candidate tile
 * @return -1 if no candidate found
 */
int CachedTiler::GetExpiredCacheTile() {
	for (int i = tile_cache_.size() - 1; i >= 0; i--) {
		if (!IsTileInMap(tile_cache_[i])) {
			return i;
		}
	}
	return -1;
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
void CachedTiler::UpdateFrameVolume(Pixel* pixels, tile_t tile) {

	// Setup start and end points
	vector<unsigned int> start, end;
	start.push_back(0); start.push_back(0); start.push_back(tile.zIndex);
	end.push_back(tile_width_ - 1); end.push_back(tile_height_ - 1); end.push_back(tile.zIndex);

	((ClNddiDisplay *)display_)->CopyPixels(pixels, start, end);
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
void CachedTiler::UpdateCoefficientMatrices(size_t x, size_t y, tile_t tile) {

	// Setup start and end points
	vector<unsigned int> start, end;
	start.push_back(x * tile_width_); start.push_back(y * tile_height_);
	end.push_back((x + 1) * tile_width_ - 1); end.push_back((y + 1) * tile_height_ - 1);
	if (end[0] >= display_->DisplayWidth()) { end[0] = display_->DisplayWidth() - 1; }
	if (end[1] >= display_->DisplayHeight()) { end[1] = display_->DisplayHeight() - 1; }

	((ClNddiDisplay *)display_)->FillCoefficient(tile.zIndex, 2, 2, start, end);
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
void CachedTiler::PushTile(tile_t tile, Pixel* pixels, size_t i, size_t j) {
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
void CachedTiler::PushTile(tile_t tile, Pixel* pixels) {
	// Create and push the start coordinates
	vector<unsigned int> start;
	start.push_back(0); start.push_back(0); start.push_back(tile.zIndex);

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
void CachedTiler::PushTile(tile_t tile, size_t i, size_t j) {
	// Create and push the position and start coordinates
	vector<unsigned int> position;
	position.push_back(2); position.push_back(2);

	vector<unsigned int> start;
	start.push_back(i * tile_width_); start.push_back(j * tile_height_);

#ifndef NO_OMP
#pragma omp critical
#endif
	{
		// Push the tile and starts
		coefficients_list.push_back(tile.zIndex);
		coefficient_positions_list.push_back(position);
		coefficient_plane_starts_list.push_back(start);
	}
}

#endif
