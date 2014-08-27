__kernel void computePixel(__global int* inputVector,
#ifdef NARROW_DATA_STORES
                           __global short* coefficientPlane,
#else
                           __global int* coefficientPlane,
#endif
                           __global uint* frameVolume,
                           __global uint* frameVolumeDims,
                           __write_only image2d_t frameBuffer,
                           const uint cm_width,
                           const uint cm_height,
                           const uint display_width,
                           const uint display_height)
{
    uint x = get_global_id(0);
    uint y = get_global_id(1);
    int2 coord = { x, y };
    uint cm_size = cm_width * cm_height;
    uint cpOffset = (y * display_width + x) * cm_size;
    uint mult = 1;
    uint fvOffset = 0;

    if ((x < display_width) && (y < display_height)) {
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
        uchar4 colori = as_uchar4(frameVolume[fvOffset]);
        float4 colorf = { (float)colori.s0/255.0f,
                          (float)colori.s1/255.0f,
                          (float)colori.s2/255.0f,
                          1.0f};
        write_imagef( frameBuffer, coord, colorf );
    }
}
