/* NB: This code is not as flexible as it could be. Non-integer dimensions
   may not be handled properly ("implemented" but not tested)
   Also: Add setup_done flag?
   Modified November 2007: Asymmetric hexagon, two different detector shapes/KL
   Also: setup mode, to compensate for deficiensies in field tables
   setup mode affects in_crystal and its helper function project_to_edges
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "detector_geometry.h"
#include "point.h"
#include "calc_signal.h"
#include "signal_calc_util.h"

/* The positions of the (6) physical corners of the crystal at the
   front of the detector, and to where the corresponding corners would
   be at the back of the crystal if the tapering was constant to
   z=zlen virtual_crystal_corners[CTYPE][IZ][NCORNERS/2] CTYPE = A/B, IZ=
   front or back corner, 
*/
static struct point virtual_crystal_corners[N_CRYSTAL_TYPES][2][NCORNERS/2];

/* project pt onto the six "edges" between front and back of hex. crystal*/
static int project_to_edges(GRETA_Siggen_Setup *setup, struct point pt, 
			    struct point projections[NCORNERS/2]);


float zmax_detector(GRETA_Siggen_Setup *setup){
  return setup->xtal_length;
}

float rmax_detector(GRETA_Siggen_Setup *setup){
  return setup->xtal_radius;
}

float hole_r_detector(GRETA_Siggen_Setup *setup){
  return setup->core_radius;
}

float hole_z_detector(GRETA_Siggen_Setup *setup){
  return setup->core_gap;
}

/* geometry_init
   reads information about detector geometry from file given by geometry_fname
   returns total number of segments, 0 for failure

   GRETINA type detectors are assumed
   Format of init file for asymmetric crystals:
   int x 2    -- number of segments, z and phi directions
   float x nsegz -- thickness for each segment, z direction
   float      -- size in z direction, mm 
   float      -- maximum radius, mm
   float      -- radius of central hole
   float      -- distance detector front => central hole, mm
   float      -- hole end curvature at radius r=rmax
   float      -- hole opening "funnel" width/depth
   float x 3  -- "virtual corners" for asymmetric crystal
    ...       -- 6 x 2 x N_CRYSTAL_TYPES lines of corner data 

   Format of init file for symmetric crystals:
   int x 2    -- number of segments, z and phi directions
   float x nsegz -- thickness for each segment, z direction
   float      -- size in z direction, mm 
   float      -- maximum radius, mm
   float      -- radius of central hole
   float      -- distance detector front => central hole, mm
   float      -- taper vertex, rel front of detector, mm
   float      -- taper angle, degrees

   # means remainder of line is comment, which will be ignored
*/


int geometry_init(GRETA_Siggen_Setup *setup){
  FILE *fp;
  char line[MAX_LINE], *cp, *next;
  int i, j, k;
  struct point pt;

  if ((fp = fopen(setup->geometry_name, "r")) == NULL){
    error("Failed to open geometry configuration file: %s\n", setup->geometry_name);
    return -1;
  }
  if (read_setup_line(fp, line) != 0 ||
      sscanf(line, "%d %d", &setup->nseg_z, &setup->nseg_phi) != 2){
    error("Failed to read number of segments from %s\n", setup->geometry_name);
    return -1;
  }
  tell(NORMAL, "Segments in z-direction: %d, phi: %d\n", setup->nseg_z, setup->nseg_phi);
  if (read_setup_line(fp, line) != 0){
    error("Failed to read segment thickness info from file: %s\n",
          setup->geometry_name);
    return -1;
  }
  if ((setup->zmax_segment = malloc(setup->nseg_z*sizeof(*setup->zmax_segment))) == NULL){
    error("Malloc failed in geometry_init\n");
    return -1; /*FIXME -- or exit?*/
  }
  cp = line;
  for (i = 0; i < setup->nseg_z; i++){
    setup->zmax_segment[i] = strtod(cp, &next);
    if (next == cp){
      error("Failed to parse segment thickness info from file: %s\n",
	    setup->geometry_name);
      return -1;
    }
    cp = next;
  }
  for (i = 1; i < setup->nseg_z; i++){
    setup->zmax_segment[i] += setup->zmax_segment[i-1];
  } 
  tell(NORMAL, "Segment max z, as a fn of segment in z direction: ");
  for (i = 0; i < setup->nseg_z; i++){
    tell(NORMAL, "%.1f ", setup->zmax_segment[i]);
  }
  tell(NORMAL,"\n");

  for (i = 0; i < N_CRYSTAL_TYPES; i++){
    for (j = 0; j < 2; j++){
      for (k = 0; k < NCORNERS/2; k++){
	if (read_setup_line(fp,line) != 0 ||
	    sscanf(line, "%f %f %f", &pt.x, &pt.y, &pt.z) != 3){
	  error("Failed to read crystal corner (%d,%d,%d) from file: %s\n",
		i, j, k, setup->geometry_name);
	  return -1;
	}
	virtual_crystal_corners[i][j][k] = pt;
      }
    }
  }
  tell(NORMAL, "Succesfully read %d corner positions\n", 
       N_CRYSTAL_TYPES*2*NCORNERS/2);
  fclose(fp);

  return setup->nseg_phi * setup->nseg_z + 1;
}

#define SQ(x) ((x)*(x))

/*returns 0 (false) or 1 (true) depending on whether pt is inside the crystal */
int in_crystal(GRETA_Siggen_Setup *setup, point pt){
  float r, r2, rmin;
  struct point ep[NCORNERS/2], vij, vpj, cp, cpold;
  int   i, j;
  float hz, dotp;
  float rmax, zmin, zmax, hr, dx, dy;
  static int first = 1, core_offset = 0;

  if (first) {
    first = 0;
    setup->core_offset_max = 0;
    if (setup->core_offset_x_top !=0 ||
        setup->core_offset_y_top !=0 ||
        setup->core_offset_x_bottom !=0 ||
        setup->core_offset_y_bottom !=0) {
      core_offset = 1;
      setup->core_offset_max = sqrt(SQ(setup->core_offset_x_top) + SQ(setup->core_offset_y_top));
      r2 = sqrt(SQ(setup->core_offset_x_bottom) + SQ(setup->core_offset_y_bottom));
      if (setup->core_offset_max < r2) setup->core_offset_max = r2;
    }
  }

  zmin = 0;
  zmax = setup->xtal_length;
  rmax = setup->xtal_radius;
  hr   = setup->core_radius;

  /* outside crystal in z direction */
  if (pt.z < zmin || pt.z >= zmax) return 0;

  /* outside maximum radius */
  r = sqrt(SQ(pt.x) + SQ(pt.y));
  if (r > rmax) return 0;

  /* inside central hole */
  if (core_offset) {
    dx = setup->core_offset_x_top +
      (setup->core_offset_x_bottom - setup->core_offset_x_top) *
      (setup->xtal_length - pt.z) / setup->core_length;
    dy = setup->core_offset_y_top +
      (setup->core_offset_y_bottom - setup->core_offset_y_top) *
      (setup->xtal_length - pt.z) / setup->core_length;
    r2 = sqrt(SQ(pt.x - dx) + SQ(pt.y - dy));
  } else {
    r2 = r;
  }
  hz = setup->core_gap + setup->core_bullet_radius;
  if (pt.z >= hz && r2 < hr) return 0;

  /* rounded end of hole */
  if (pt.z < hz && pt.z >= setup->core_gap) {
    rmin = sqrt(SQ(setup->core_bullet_radius) - SQ(hz - pt.z)) + hr - setup->core_bullet_radius;
    if (r2 < rmin) return 0;
  }

  /* widening at hole opening -- only for nonsymm crystal */
  if (pt.z > setup->xtal_length - setup->bottom_taper_length) {
    rmin = setup->core_radius + setup->bottom_taper_width *
      (pt.z - (setup->xtal_length - setup->bottom_taper_length)) / setup->bottom_taper_length;
    if (r2 < rmin) return 0;
  }

  /* find points corresponding to same "z" between pairs 
     of virtual corners @ front & back of detector => points at 
     "edges" of crystal */
  project_to_edges(setup, pt, ep);

  /* is pt inside the hexagon formed by ep? */
  for (i = 0; i < NCORNERS/2; i++) {
    j = (i + 1)%6;
    vij = vector_sub(ep[j], ep[i]);
    vpj = vector_sub(pt, ep[j]);
    cp = cross_prod(vij,vpj);
    if (i == 0) {
      cpold = cp;
    } else {
      dotp = dot_prod(cp,cpold);
      if (dotp < 0) return 0; /*outside crystal*/
    }
  }

  return 1;
}
#undef SQ

static int project_to_edges(GRETA_Siggen_Setup *setup, struct point pt, struct point projections[NCORNERS/2]) {
  int i;
  struct point vc1, vc2, dvc;
  float f;

  for (i = 0; i < NCORNERS/2; i++) {
    vc1 = setup->virtual_crystal_corners[setup->xtal_type][0][i];
    vc2 = setup->virtual_crystal_corners[setup->xtal_type][1][i];
    f = (pt.z - vc1.z)/(vc2.z - vc1.z);
    dvc = vector_sub(vc2, vc1);
    projections[i] = vector_add(vc1, vector_scale(dvc,f));
  }
  return 0;
}
/*
static point move_towards_axis(point pt, float distance){
  point new_pt;
  vector v,dv;

  v = pt;
  v.z = 0;
  dv = vector_scale(v, -distance/vector_length(v));
  new_pt = vector_add(pt, dv); 

  return new_pt;
}
*/
