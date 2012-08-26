/*
 *  Rewinder.cpp
 *  pixelbridge
 *
 *  Created by Dave Estes on 3/20/11.
 *  Copyright 2011 Dave Estes. All rights reserved.
 *
 */
#include "Rewinder.h"

Rewinder::Rewinder(size_t frame_count, size_t frame_width, size_t frame_height)
:	frame_count_(frame_count),
	frame_width_(frame_width),
	frame_height_(frame_height)
{
	frames_ = (uint8_t**)malloc(frame_count_ * sizeof(int*));
	for (int i = 0; i < frame_count_; i++) {
		frames_[i] = (uint8_t*)malloc(frame_width_ * frame_height_ * sizeof(uint8_t) * VIDEO_PIXEL_SIZE);
	}
}


Rewinder::~Rewinder() {
	for (int i = 0; i < frame_count_; i++) {
		free(frames_[i]);
	}
	free(frames_);
}


void Rewinder::CopyFrame(uint8_t* source_frame, size_t frame_number) {
	memcpy(frames_[frame_number], source_frame, frame_width_ * frame_height_ * sizeof(uint8_t) * VIDEO_PIXEL_SIZE);
}


uint8_t* Rewinder::GetFrame(size_t frame_number) {
	return frames_[frame_number];
}

