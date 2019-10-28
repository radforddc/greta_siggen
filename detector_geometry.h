/* detector_geometry.h -- based on m3d2s.f by I-Yang Lee
 * Karin Lagergren
 *
 * This module keeps track of the detector geometry
 * See also: assym_detector.h
 */
#ifndef _DETECTOR_GEOMETRY_H
#define _DETECTOR_GEOMETRY_H

#include "greta_siggen.h"
#include "point.h"

/* geometry_init
   reads information about detector geometry from file given by geometry_fname
   returns total number of segments, -1 for failure
*/
int geometry_init(GRETA_Siggen_Setup *setup);

/* in_crystal
   returns 1 if inside crystal, 0 otherwise
*/
int in_crystal(GRETA_Siggen_Setup *setup, point pt);

/* zmax_detector
   returns the maximum z value for detector
*/
float zmax_detector(GRETA_Siggen_Setup *setup);

/* rmax_detector
   returns the maximum radius for detector
*/
float rmax_detector(GRETA_Siggen_Setup *setup);

/* hole_r_detector
   returns the radius of the central hole
*/
float hole_r_detector(GRETA_Siggen_Setup *setup);

/* hole_z_detector
   returns the z value for the start ("bottom") of the central hole
*/
float hole_z_detector(GRETA_Siggen_Setup *setup);

#endif /*#ifndef _DETECTOR_GEOMETRY_H*/
