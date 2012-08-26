#ifndef REWINDER_H
#define REWINDER_H
/*
 *  Rewinder.h
 *  pixelbridge
 *
 *  Created by Dave Estes on 3/20/11.
 *  Copyright 2011 Dave Estes. All rights reserved.
 *
 */
#include <stdlib.h>
#include <string.h>

#include "FfmpegPlayer.h"

class Rewinder {
	
public:
	Rewinder(size_t frame_count, size_t frame_width, size_t frame_height);
	~Rewinder();
	void CopyFrame(uint8_t* source_frame, size_t frame_number);
	uint8_t* GetFrame(size_t frame_number);
	
private:
	size_t frame_count_, frame_width_, frame_height_;
	uint8_t** frames_;
};

#endif // REWINDER_H