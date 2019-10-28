/* fields.h -- based on m3d2s.f by I-Yang Lee
 * Karin Lagergren
 *
 * This module handles the electric field and weighting potential and 
 * calculates drift velocities
 */

#ifndef _FIELDS_H
#define _FIELDS_H

/* calculate anisotropic drift velocities? (vel. depends on angle between
   el. field and crystal axis; otherwise the velocity will always be 
   in the direction of the el. field 
*/
#define DRIFT_VEL_ANISOTROPY 1
#include "greta_siggen.h"
#include "point.h"

/* field_setup
   given a field directory file, read electic field and weighting
   potential tables from files listed in directory
   check that number of segments matches nsegs
   returns 0 for success
*/
int field_setup(GRETA_Siggen_Setup *setup);

/* free malloc()'ed memory and do other cleanup*/
int fields_finalize(GRETA_Siggen_Setup *setup);

/* wpotentials
   gives (interpolated or extrapolated ) weighting potential for each segment
   at point pt. These values are stored in wp, which is assumed to 
   have one element per segment.
   returns 0 for success, 1 on failure.
*/
int wpotentials(GRETA_Siggen_Setup *setup, point pt, float *wp);

/* drift_velocity
   calculates drift velocity for charge q at point pt
   returns 0 on success, 1 if successful but extrapolation was needed,
   and -1 for failure
*/
int drift_velocity(GRETA_Siggen_Setup *setup, point pt, float q, vector *velocity);

/*set detector temperature. 77F (no correction) is the default
   MIN_TEMP & MAX_TEMP defines allowed range*/
void set_temp(GRETA_Siggen_Setup *setup, float temp);


#endif /*#ifndef _FIELDS_H*/
