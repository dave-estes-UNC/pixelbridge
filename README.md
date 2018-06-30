pixelbridge
===========

DEPRECATED: NDDI is moved into a seperate project and the pixelbridge application
is now a sample nddiwall application, additionally under a seperate project.

    https://github.com/dave-estes-UNC/nddi
    https://github.com/dave-estes-UNC/nddiwall

This is the original implementation of a simulated NDDI display embodied in an
application called pixelbridge. Pixelbridge represents the backward comaptible
use case of simply "bridging" pixels from a traditional framebuffer application
to an NDDI display using unique configurations of the display to achieve a cost
savings for pixel transmission.

Building
========

The project does still build, though OpenCL isn't working lately. To build

    $ make NO_CL=1
    $ make NO_CL=1 debug