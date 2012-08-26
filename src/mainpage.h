/**
 * @mainpage n-Dimensional Display Interface (NDDI): Pixel Bridge Use Case
 *
 * @authors Dave Estes <cdestes@email.unc.edu>
 *
 * @section overview Overview
 *
 * The n-Dimensional Display Interface (NDDI) is a new approach to getting pixels to a display. Typically
 * displays use a frame buffer for receiving pixels to be displayed. That frame buffer must be updated at
 * a frequency equal to the desired refresh rate.
 *
 * Modern displays implement features such as double buffering to prevent to prevent partial updates or even higher
 * order functions to fill or blit. In each case, the entire framebuffer would not need to be updated on every cycle.
 *
 * NDDI builds on these concepts in a simple-to-implement hardware solution.
 *
 * @section frameVolume Pixel Storage: Frame Volume
 *
 * The n-Dimensional aspect of NDDI refers to the Frame Volume. Instead of organizing a display device's
 * pixels into a two-dimensional frame buffer and updating all of those pixels for every frame, a multi-dimensional
 * frame volume is created. The frame volume also contains pixels, but there is a highly-configurable mapping
 * that takes place, instead of the fixed mapping used with a traditional frame buffer.
 *
 * @section coefficientPlane Mapping: Coefficient Plane
 *
 * The mapping is configured through the use of a coefficient plane. The coefficient plane, is a 
 * is a two dimensional array of coefficient matrices. The dimensions of the coefficient plane
 * match the dimensions of the viewable area, which is typically the size of the actual display.
 * The size of the coefficient matrices are dependant on the frame volume and the input vector size.
 * the number of columns in each coefficient matrix matches the size of the input vector. The number
 * of rows matches the number of dimensions in the frame volume.
 *
 * @section inputVector Driving the Display: Input Vector
 *
 * The display output is driven by the input vector. The first two values in the input correspond to
 * to a pixel location on the display as well as the corresponding coefficient matrix in the coefficient
 * plane. These first two values are not driven by the nddi client, but are used by the display to refresh itself.
 * If the input vector is larger than two, then the nddi client can use them to drive particular mappings.
 * As an example, a third value can be used as a clock tick that will help animated the frames of a sprite.
 *
 * @section pixelBridge Pixel Brige Use Case
 *
 * Pixel Brige is the first use case, used to study the viability of NDDI. The goal of Pixel Bridge is is to take a
 * recording computing session and display it on an NDDI display. FFMPEG is leveraged to playback VMNC recordings
 * created with VMWare. VMNC is an RFB stream encapsulated in an AVI containter. RFB is the same protocol used by
 * VNC. Using FFMPEG also allows the playback of another other video format supported by FFMPEG.
 *
 * The NDDI display was implemented in such a way that it returns a framebuffer that is rendered into an OpenGL
 * window. The parallel calculations were sped up with OpenMP. OpenCL was used as an alternative, but the latest
 * status is that it's suffering a lot of rendering errors.
 *
 * The NDDI display is configured as a simple frame buffer, as a flat tiled display, and as a cached tiled display.
 * The FlatTiler and CachedTiler classes are responsible for taking a decoded frame and tiling it. They then
 * update the attached NDDI.
 *
 * Additionally the Pixel Brige can calculate the estimated cost of the same number duration of video renderer
 * at 60 Hz. This serves as a lower bounds, which all three configurations will outperform. For the upper boundary,
 * Pixel Bridge can return the number of bytes (four per pixel) changed. This would be the cost associated
 * if "Perfect Pixel Latching" without any addressing information was possible.
 */
