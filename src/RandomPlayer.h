//
//  RandomPlayer.h
//  pixelbridge
//
//  Created by Dave Estes on 8/26/13.
//  Copyright (c) 2013 Dave Estes. All rights reserved.
//

#ifndef pixelbridge_RandomPlayer_h
#define pixelbridge_RandomPlayer_h

#include "Player.h"

#define VIDEO_PIXEL_SIZE   3 // RGB
#define FIXED_WIDTH        64
#define FIXED_HEIGHT       64
#define FIXED_FRAME_COUNT  100


/**
 * This class will utilize an FFMPEG demuxer, decoder, and color convertor to decode
 * each frame of the provided video.
 */
class RandomPlayer : Player {
public:
	/**
	 * The constructor for the player.
	 */
	RandomPlayer();
    
	~RandomPlayer();
    
	/**
	 * Returns the width of the video.
	 *
	 * @returns The width of the video
	 */
	size_t width();
    
	/**
	 * Returns the height of the video.
	 *
	 * @returns The height of the video
	 */
	size_t height();
    
	/**
	 * Decodes a frame.
	 *
	 * @returns A pointer to the decoded frame.
	 */
	uint8_t* decodeFrame();

private:
    void createRandomFrame();

    
private:
	size_t width_;
	size_t height_;
    uint32_t         framesRemaining_;
    
	int              numBytes_;
	uint8_t*         buffer_;
};


#endif
