// Must match struct in ClNddiDisplay.h
typedef struct {
#ifdef NARROW_DATA_STORES
    short  coefficient;
#else
    int  coefficient;
#endif
    uint posCol;
    uint posRow;
    uint startX;
    uint startY;
    uint sizeW;
    uint sizeH;
} coefficient_update_t;

__kernel void fillCoefficient(__global coefficient_update_t* packet,
#ifdef NARROW_DATA_STORES
                              __global short* coefficientPlane,
#else
                              __global int* coefficientPlane,
#endif
                              const uint cm_width,
                              const uint cm_height,
                              const uint display_width,
                              const uint display_height,
                              uint num)
{
    uint i = get_global_id(0);
    uint cpOffset;
    coefficient_update_t upd;

    if (i < num) {
        upd = packet[i];
        for (uint y = upd.startY; y < upd.startY + upd.sizeH; y++) {
            for (uint x = upd.startX; x < upd.startX + upd.sizeW; x++) {
                cpOffset = (cm_width * cm_height) * (y * display_width + x) + upd.posRow * cm_width + upd.posCol;
                coefficientPlane[cpOffset] = upd.coefficient;
            }
        }
    }
}
