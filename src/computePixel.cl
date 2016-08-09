inline uchar CLAMP_SIGNED_BYTE(int i) {
    uint ret;

    if (i < 0) {
        ret = 0;
    } else if (i > 0xff) {
        ret = 0xff;
    } else {
        ret = i;
    }

    return (uchar)ret;
}

inline uchar CLAMP_UNSIGNED_BYTE(uint i) {
    uint ret;

    if (i > 0xff) {
        ret = 0xff;
    } else {
        ret = i;
    }

    return ret;
}

__kernel void computePixel(__global int* inputVector,
#ifdef NARROW_DATA_STORES
                           __global short* coefficientPlane,
#else
                           __global int* coefficientPlane,
#endif
#ifdef NARROW_DATA_STORES
                           __global short* scalers,
#else
                           __global int* scalers,
#endif
                           __global uint* frameVolume,
                           __global uint* frameVolumeDims,
                           __write_only image2d_t frameBuffer,
                           const uint cm_width,
                           const uint cm_height,
                           const uint display_width,
                           const uint display_height,
                           const uint display_planes,
                           const uint shifter,
                           const int do_signed)
{
    uint x = get_global_id(0);
    uint y = get_global_id(1);
    int2 coord = { x, y };
    uint cm_size = cm_width * cm_height;
    uint mult = 1;
    uint fvOffset = 0;
    int4 accumulator = {0, 0, 0, 255};

    if ((x < display_width) && (y < display_height)) {

        // Calculate the pixel for each plane and add to the accumulators
        for (uint p = 0; p < display_planes; p++) {
            // Calculate the scaler and coefficient plane offsets
            uint scOffset = (p * display_height * display_width + y * display_width + x);
            uint cpOffset = (p * display_height * display_width + y * display_width + x) * cm_size;

            // Calculate the frame volume offset
            for (uint cmj = 0; cmj < cm_height; cmj++) {
                uint s = 0;
                s += coefficientPlane[cpOffset] * x; cpOffset++;
                s += coefficientPlane[cpOffset] * y; cpOffset++;
                for (uint cmi = 2; cmi < cm_width; cmi++) {
                    s += coefficientPlane[cpOffset] * inputVector[cmi]; cpOffset++;
                }
                fvOffset += s * mult;
                mult *= frameVolumeDims[cmj];
            }

            // Multiply the pixel by the scaler and add to the accumulator
#ifdef USE_ALPHA_CHANNEL
            if (!do_signed) {
                uchar4 pixel = as_uchar4(frameVolume[fvOffset]);

                accumulator.s0 += pixel.s0 * pixel.s3 * scalers[scOffset+0];
                accumulator.s1 += pixel.s1 * pixel.s3 * scalers[scOffset+1];
                accumulator.s2 += pixel.s2 * pixel.s3 * scalers[scOffset+2];
            } else {
                char4 pixel = as_char4(frameVolume[fvOffset]);

                accumulator.s0 += pixel.s0 * pixel.s3 * scalers[scOffset+0];
                accumulator.s1 += pixel.s1 * pixel.s3 * scalers[scOffset+1];
                accumulator.s2 += pixel.s2 * pixel.s3 * scalers[scOffset+2];
            }
#else
            if (!do_signed) {
                uchar4 pixel = as_uchar4(frameVolume[fvOffset]);

                accumulator.s0 += pixel.s0 * scalers[scOffset+0];
                accumulator.s1 += pixel.s1 * scalers[scOffset+1];
                accumulator.s2 += pixel.s2 * scalers[scOffset+2];
            } else {
                char4 pixel = as_char4(frameVolume[fvOffset]);

                accumulator.s0 += pixel.s0 * scalers[scOffset+0];
                accumulator.s1 += pixel.s1 * scalers[scOffset+1];
                accumulator.s2 += pixel.s2 * scalers[scOffset+2];
            }
#endif
        }

        // Shift the accumulators
#ifdef USE_ALPHA_CHANNEL
        if (!do_signed) {
            accumulator.s0 = CLAMP_UNSIGNED_BYTE(accumulator.s0 >> (8 + shifter));
            accumulator.s1 = CLAMP_UNSIGNED_BYTE(accumulator.s1 >> (8 + shifter));
            accumulator.s2 = CLAMP_UNSIGNED_BYTE(accumulator.s2 >> (8 + shifter));
        } else {
            accumulator.s0 = CLAMP_SIGNED_BYTE(accumulator.s0 >> (8 + shifter));
            accumulator.s1 = CLAMP_SIGNED_BYTE(accumulator.s1 >> (8 + shifter));
            accumulator.s2 = CLAMP_SIGNED_BYTE(accumulator.s2 >> (8 + shifter));
        }
#else
        if (!do_signed) {
            accumulator.s0 = CLAMP_UNSIGNED_BYTE(accumulator.s0 >> shifter);
            accumulator.s1 = CLAMP_UNSIGNED_BYTE(accumulator.s1 >> shifter);
            accumulator.s2 = CLAMP_UNSIGNED_BYTE(accumulator.s2 >> shifter);
        } else {
            accumulator.s0 = CLAMP_SIGNED_BYTE(accumulator.s0 >> shifter);
            accumulator.s1 = CLAMP_SIGNED_BYTE(accumulator.s1 >> shifter);
            accumulator.s2 = CLAMP_SIGNED_BYTE(accumulator.s2 >> shifter);
        }
#endif

        uchar4 colori = {accumulator.s0, accumulator.s1, accumulator.s2, accumulator.s3};
        float4 colorf = { (float)colori.s0/255.0f,
                          (float)colori.s1/255.0f,
                          (float)colori.s2/255.0f,
                          1.0f};
        write_imagef( frameBuffer, coord, colorf );
    }
}
