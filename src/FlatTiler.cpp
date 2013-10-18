/*
 *  FlatTiler.cpp
 *  pixelbridge
 *
 *  Created by Dave Estes on 11/26/10.
 *  Copyright 2010 Dave Estes. All rights reserved.
 *
 */

#include <iostream>
#include <zlib.h>

#include "PixelBridgeFeatures.h"
#include "FlatTiler.h"


/**
 * The FlatTiler is created based on the dimensions of the NDDI display that's passed in. If those
 * dimensions change, then the FlatTiler should be destroyed and re-created.
 */
FlatTiler::FlatTiler (size_t display_width, size_t display_height,
                      size_t tile_width, size_t tile_height,
                      size_t bits, bool quiet)
: tile_width_(tile_width),
  tile_height_(tile_height),
  bits_(bits),
  quiet_(quiet)
{
	
    // 2 dimensional matching the Video Width x Height
    vector<unsigned int> fvDimensions;
    fvDimensions.push_back(display_width);
    fvDimensions.push_back(display_height);

#ifndef NO_CL
    display_ = new ClNddiDisplay(fvDimensions,                  // framevolume dimensional sizes
                                 display_width, display_height, // display size
                                 2); 						    // input vector size (x and y only)
#else
    display_ = new GlNddiDisplay(fvDimensions,                   // framevolume dimensional sizes
                                 display_width, display_height,  // display size
                                 2); 						     // input vector size (x and y only)
#endif
    
	// Compute tile_map width
	tile_map_width_ = display_->DisplayWidth() / tile_width;
	if ((tile_map_width_ * tile_width) < display_->DisplayWidth()) { tile_map_width_++; }
	
	// Compute tile_map height
	tile_map_height_ = display_->DisplayHeight() / tile_height;
	if ((tile_map_height_ * tile_height) < display_->DisplayHeight()) { tile_map_height_++; }
	
	// Set up the tile map, one column at a time with a checksum of zero.
	tile_map_.resize(tile_map_width_);
	for (int i = 0; i < tile_map_width_; i++) {
		for (int j = 0; j < tile_map_height_; j++) {
			tile_map_[i].push_back(0L);
		}
	}
    
	// Set tile update count to zero
	unchanged_tiles_ = tile_updates_ = 0;
    
    // Not Input Vector Initialization required

    // Initialize Frame Volume
    nddi::Pixel p;
    p.r = p.g = p.b = p.a = 0xff;
    vector<unsigned int> start, end;
    start.push_back(0); start.push_back(0);
    end.push_back(display_width - 1); end.push_back(display_height - 1);
    display_->FillPixel(p, start, end);
    
    // Initialize Coefficient Planes
    InitializeCoefficientPlanes();
}

/**
 * Returns the Display created and initialized by the tiler.
 */
GlNddiDisplay* FlatTiler::GetDisplay() {
    return display_;
}

/**
 * Initializes the Coefficient Planes for this tiler.
 */
void FlatTiler::InitializeCoefficientPlanes() {

    vector< vector<int> > coeffs;
    coeffs.resize(2);
    coeffs[0].push_back(1); coeffs[0].push_back(0);
    coeffs[1].push_back(0); coeffs[1].push_back(1);
    
	vector<unsigned int> start, end;
    start.push_back(0); start.push_back(0); start.push_back(0);
    end.push_back(display_->DisplayWidth() - 1); end.push_back(display_->DisplayHeight() - 1); end.push_back(0);

    display_->FillCoefficientMatrix(coeffs, start, end);

    // Turn off all planes and then set the 0 plane to full on.
    end[2] = NUM_COEFFICIENT_PLANES - 1;
    display_->FillScaler(0, start, end);
    end[2] = 0;
    display_->FillScaler(NUM_COEFFICIENT_PLANES, start, end);
}

/**
 * Update the tile_map, tilecache, and then the NDDI display based on the frame that's passed in. The
 * frame is returned from the ffmpeg player as an RGB buffer. There is not Alpha channel.
 * 
 * @param buffer Pointer to an RGB buffer
 * @param width The width of the RGB buffer
 * @param height The height of the RGB buffer
 */
void FlatTiler::UpdateDisplay(uint8_t* buffer, size_t width, size_t height)
{
	int                            unchanged = 0;
	int                            updates = 0;
	unsigned char                  mask = 0xff << (8 - bits_);
    Pixel                         *tile_pixels = NULL;
    Pixel                         *tile_pixels_sig_bits = NULL;
#ifdef USE_COPY_PIXEL_TILES
	vector<Pixel *>                tiles;
	vector<vector<unsigned int> >  starts;
#endif
	
	// Break up the passed in buffer into one tile at a time
#ifndef NO_OMP
#pragma omp parallel for ordered
#endif // !NO_OMP
	for (int j_tile_map = 0; j_tile_map < tile_map_height_; j_tile_map++) {

#ifndef NO_OMP
#pragma omp ordered
#endif
        {
            for (int i_tile_map = 0; i_tile_map < tile_map_width_; i_tile_map++) {
                // Initialize a tile's checksum
                unsigned long tile_checksum = 0L;
                
                // Use locals for this tile's width and height in case they need to be adjust at the edges
                int tw = tile_width_, th = tile_height_;
                if (i_tile_map == (tile_map_width_ - 1)) {
                    tw -= tile_map_width_ * tile_width_ - width;
                }
                if (j_tile_map == (tile_map_height_ - 1)) {
                    th -= tile_map_height_ * tile_height_ - height;
                }
                
                // Allocate tiles is necessary. Sometimes they're re-used.
                if (!tile_pixels)
                	tile_pixels = (Pixel*)malloc(tw * th * sizeof(Pixel));
                if (!tile_pixels_sig_bits)
                	tile_pixels_sig_bits = (Pixel*)malloc(tw * th * sizeof(Pixel));
                
                // Build the tile's pixel array while computing the checksum in an alternative tile which only
                // holds the significant bits
                for (int j_tile = 0; j_tile < th; j_tile++) {
                    // Compute the offset into the RGB buffer for this row in this tile
                    int bufferOffset = 3 * ((j_tile_map * tile_height_ + j_tile) * width + (i_tile_map * tile_width_));
                    
                    for (int i_tile = 0; i_tile < tw; i_tile++) {
                        Pixel p, psb;
                        
                        p.r = buffer[bufferOffset++];
                        p.g = buffer[bufferOffset++];
                        p.b = buffer[bufferOffset++];
                        p.a = 0xff;
                        tile_pixels[j_tile * tw + i_tile].packed = p.packed;

                        psb.r = p.r & mask;
                        psb.g = p.g & mask;
                        psb.b = p.b & mask;
                        psb.a = p.a & mask;
                        tile_pixels_sig_bits[j_tile * tw + i_tile].packed = psb.packed;
                    }
                }
                
#if (CHECKSUM_CALCULATOR == TRIVIAL)
                tile_checksum  = (unsigned long)tile_pixels_sig_bits[0].packed << 32;
                tile_checksum |= (unsigned long)tile_pixels_sig_bits[tw * th - 1].packed;
#else
                unsigned long crc = crc32(0L, Z_NULL, 0);
#if (CHECKSUM_CALCULATOR == CRC)
                tile_checksum = crc32(crc, (unsigned char*)tile_pixels_sig_bits, tw * th * sizeof(Pixel));
#elif (CHECKSUM_CALCULATOR == ADLER)
                tile_checksum = adler32(crc, (unsigned char*)tile_pixels_sig_bits, tw * th * sizeof(Pixel));
#endif
#endif
                
                // If the checksum in the tile map doesn't match, then update the frame volume
                if (tile_map_[i_tile_map][j_tile_map] != tile_checksum) {
                    tile_map_[i_tile_map][j_tile_map] = tile_checksum;
#ifdef USE_COPY_PIXEL_TILES
#ifndef NO_OMP
#pragma omp critical
#endif
                    {
                    	// Push the tile
                    	tiles.push_back(tile_pixels);

                    	// Force new tile to be allocated. This one will be freed after it's copied
                    	tile_pixels = NULL;

                    	// Create and push the start coordinates
                    	vector<unsigned int> start;
                    	start.push_back(i_tile_map * tile_width_); start.push_back(j_tile_map * tile_height_);
                    	starts.push_back(start);
                    }

#else
                    UpdateFrameVolume(tile_pixels, i_tile_map, j_tile_map);
#endif

#ifndef NO_OMP
#pragma omp atomic
#endif
                    updates++;
                } else {
#ifndef NO_OMP
#pragma omp atomic
#endif
                    unchanged++;
                }
            }
        }
    }

#ifdef USE_COPY_PIXEL_TILES
	// If any tiles were updated
	if (updates > 0) {

		// Update the Frame Volume by copying the tiles over
		vector<unsigned int> size;
		size.push_back(tile_width_); size.push_back(tile_height_);
		display_->CopyPixelTiles(tiles, starts, size);

		while (!tiles.empty()) {
			free(tiles.back());
			tiles.pop_back();
		}
	}
#endif

    // Free alloc'd memory
	if (tile_pixels)
		free(tile_pixels);
	if (tile_pixels_sig_bits)
		free(tile_pixels_sig_bits);

	// Report update statistics
	unchanged_tiles_ += unchanged;
	tile_updates_ += updates;

	if (!quiet_) {
		cout << "Flat Tiling Statistics:" << endl << "  unchanged tiles: " << unchanged_tiles_ << " tiles updated: " << tile_updates_ << endl;
	}
}

/**
 * Updates region of the Frame Volume corresponding to the tile's i and j location.
 *
 * @param pixels The array of pixels to use in the update.
 */
#ifndef USE_COPY_PIXEL_TILES
void FlatTiler::UpdateFrameVolume(Pixel* pixels, int i_map, int j_map) {
	
	// Setup start and end points
	vector<unsigned int> start, end;
	start.push_back(i_map * tile_width_); start.push_back(j_map * tile_height_);
	end.push_back(i_map * tile_width_ + tile_width_ - 1); end.push_back(j_map * tile_height_ + tile_height_ - 1);
	if (end[0] >= display_->DisplayWidth()) { end[0] = display_->DisplayWidth() - 1; }
	if (end[1] >= display_->DisplayHeight()) { end[1] = display_->DisplayHeight() - 1; }
    
    display_->CopyPixels(pixels, start, end);
}
#endif
