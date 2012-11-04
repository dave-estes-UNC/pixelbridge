__kernel void fillCoefficient(__global uint* packet,
                              __write_only int* coefficientPlane,
                              const uint cm_width,
                              const uint cm_height,
                              const uint display_width,
                              const uint display_height,
                              const uint num)
{
    uint i = get_global_id(0);
    if (i < num) {
    }
}
