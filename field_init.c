#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

float fminf(float x, float y);

#include "greta_siggen.h"

static int in_crystal(GRETA_Siggen_Setup *setup, point pt);
static int project_to_edges(GRETA_Siggen_Setup *setup, struct point pt, struct point projections[NCORNERS/2]);
static int segment_number(GRETA_Siggen_Setup *setup, point pt);
static int segment_number_new(GRETA_Siggen_Setup *setup, struct point pt,
                              float *distance_xy, float *distance_z, int *seg2_xy, int *seg2_z);


/* -------------------------------------- dist_from_contact ------------------- */
float dist_from_contact(point pt, point delta, GRETA_Siggen_Setup *setup) {
  float  factor = 1, d = 0.5;
  point test;
  int    n;

  for (n=0; n<7; n++) {  // 7 steps => 1/128 precision
    test.x = pt.x + factor * delta.x;
    test.y = pt.y + factor * delta.y;
    test.z = pt.z + factor * delta.z;
    if (!in_crystal(setup, test)) {
      factor -= d;
    } else {
      if (n == 0) return -1.0;
      factor += d;
    } 
    d /= 2.0;
  }
  return factor;
} /* dist_from_contact */

#define SQ(x) ((x)*(x))

int geometry_init(GRETA_Siggen_Setup *setup) {
  FILE  *fp;
  char  line[MAX_LINE], *cp, *next;
  int   i, j, k, n;
  float r, d, grid = setup->xtal_grid, g = grid * 0.05;
  point pt, tp[6];

  if ((fp = fopen(setup->geometry_name, "r")) == NULL) {
    error("Failed to open geometry configuration file: %s\n", setup->geometry_name);
    return -1;
  }
  if (read_setup_line(fp, line) != 0 ||
      sscanf(line, "%d %d", &setup->nseg_z, &setup->nseg_phi) != 2) {
    error("Failed to read number of segments from %s\n", setup->geometry_name);
    return -1;
  }
  printf( "Segments in z-direction: %d, phi: %d\n", setup->nseg_z, setup->nseg_phi);
  if (read_setup_line(fp, line) != 0) {
    error("Failed to read segment thickness info from file: %s\n", 
	  setup->geometry_name);
    return -1;
  }
  if ((setup->zmax_segment = malloc(setup->nseg_z*sizeof(*setup->zmax_segment))) == NULL) {
    error("Malloc failed in geometry_init\n");
    return -1; /*FIXME -- or exit?*/
  }
  cp = line;
  for (i = 0; i < setup->nseg_z; i++) {
    setup->zmax_segment[i] = strtod(cp, &next);
    if (next == cp) {
      error("Failed to parse segment thickness info from file: %s\n",
	    setup->geometry_name);
      return -1;
    }
    cp = next;
  }
  for (i = 1; i < setup->nseg_z; i++) {
    setup->zmax_segment[i] += setup->zmax_segment[i-1];
  } 
  printf( "Segment max z, as a fn of segment in z direction: ");
  for (i = 0; i < setup->nseg_z; i++) {
    printf( "%.1f ", setup->zmax_segment[i]);
  }
  printf("\n");

  for (i = 0; i < N_CRYSTAL_TYPES; i++) {
    for (j = 0; j < 2; j++) {
      for (k = 0; k < NCORNERS/2; k++) {
	if (read_setup_line(fp,line) != 0 ||
	    sscanf(line, "%f %f %f", &pt.x, &pt.y, &pt.z) != 3) {
	  error("Failed to read crystal corner (%d,%d,%d) from file: %s\n",
		i, j, k, setup->geometry_name);
	  return -1;
	}
	setup->virtual_crystal_corners[i][j][k] = pt;
      }
    }
  }
  printf( "Succesfully read %d corner positions\n", N_CRYSTAL_TYPES*2*NCORNERS/2);
  fclose(fp);

  /* find voxel types : OUTSIDE, CONTACT_0, CONTACT_VB, FIXED, INSIDE, CONTACT_EDGE */

  printf("\n  Initializing points in crystal for x = %.1f to %.1f, grid %.2f...\n",
         setup->xmin, setup->xmax, grid);
  for (i = 0; i < setup->numx; i++) {
    printf("\r %d/%d", i, setup->numx-1);
    pt.x = setup->x0 + i*grid;
    for (n=0; n<6; n++) tp[n].x = pt.x;
    tp[0].x -= g; tp[1].x += g;
    for (j = 0; j < setup->numy; j++) {
      pt.y = setup->y0 + j*grid;
      for (n=0; n<6; n++) tp[n].y = pt.y;
      tp[2].y -= g; tp[3].y += g;
      r = sqrt(pt.x*pt.x + pt.y*pt.y);
      for (k = 0; k < setup->numz; k++) {
	pt.z = setup->zmin + k*grid;
        for (n=0; n<6; n++) tp[n].z = pt.z;
        tp[4].z -= g; tp[5].z += g;
        setup->point_type[i][j][k] = INSIDE;      // start by assuming voxel is all inside bulk

        /* see if pixel is ouside (or very nearly outside) the detector bulk */
        if (i == 0 ||
            !in_crystal(setup, tp[0]) || !in_crystal(setup, tp[1]) ||
            !in_crystal(setup, tp[2]) || !in_crystal(setup, tp[3]) ||
            !in_crystal(setup, tp[4]) || !in_crystal(setup, tp[5])) {

	  setup->point_type[i][j][k] = OUTSIDE;   // voxel is all outside
	  if (pt.z < grid) {            /* front face */
	    setup->point_type[i][j][k] = CONTACT_0;
	  } else if (r <= (setup->core_radius + setup->core_offset_max +
                           setup->bottom_taper_width + grid) &&
                     k < setup->numz-1) {          /* core contact; see also below for k == setup->numz-1 */
	    setup->point_type[i][j][k] = CONTACT_VB;
	  } else if (pt.z >= setup->xtal_length - grid &&
                     r <= setup->xtal_radius - grid) {     /* passivated surface */
            setup->point_type[i][j][k] = INSIDE;
	  } else {                                                     /* outside contact */
	    setup->point_type[i][j][k] = CONTACT_0;
	  }
	}
      }
    }
  }

  // top of central hole should be CONTACT_VB
  k = setup->numz-1;
  for (i = 1; i < setup->numx-1; i++) {
    pt.x = setup->x0 + i*grid;
    for (j = 1; j < setup->numy-1; j++) {
      pt.y = setup->y0 + j*grid;
      r = sqrt(SQ(pt.x - setup->core_offset_x_top) + SQ(pt.y - setup->core_offset_y_top));
      if (r <= setup->core_radius + setup->bottom_taper_width) setup->point_type[i][j][k] = CONTACT_VB;
    }
  }

  /* find the pixels next to the contact surfaces */
  for (n=0; n<6; n++) tp[n].x = tp[n].y = tp[n].y = 0;
  tp[0].x = grid; tp[1].x = -grid;
  tp[2].y = grid; tp[3].y = -grid;
  tp[4].z = grid; tp[5].z = -grid;
  for (i = 1; i < setup->numx-1; i++) {
    printf("\r %d/%d", i, setup->numx-1);fflush(stdout);
    pt.x = setup->x0 + i*grid;
    for (j = 1; j < setup->numy-1; j++) {
      pt.y = setup->y0 + j*grid;
      for (k = 1; k < setup->numz; k++) {
	pt.z = setup->zmin + k*grid;
        setup->dx[0][i][j][k] = setup->dx[1][i][j][k] = 1.0;
        setup->dy[0][i][j][k] = setup->dy[1][i][j][k] = 1.0;
        setup->dz[0][i][j][k] = setup->dz[1][i][j][k] = 1.0;

        if (setup->point_type[i][j][k] == INSIDE &&
            (setup->point_type[i+1][j][k] < INSIDE || setup->point_type[i-1][j][k] < INSIDE ||
             setup->point_type[i][j+1][k] < INSIDE || setup->point_type[i][j-1][k] < INSIDE ||
             (k < setup->numz-1 && setup->point_type[i][j][k+1] < INSIDE) ||
             setup->point_type[i][j][k-1] < INSIDE )) {
          setup->point_type[i][j][k] = CONTACT_EDGE;
          /* find distance to contact surface */
          if (setup->point_type[i+1][j][k] < INSIDE && (d = dist_from_contact(pt, tp[0], setup)) > 0)
            setup->dx[1][i][j][k] = 1.0/d;
          if (setup->point_type[i-1][j][k] < INSIDE && (d = dist_from_contact(pt, tp[1], setup)) > 0)
            setup->dx[0][i][j][k] = 1.0/d;
          if (setup->point_type[i][j+1][k] < INSIDE && (d = dist_from_contact(pt, tp[2], setup)) > 0)
            setup->dy[1][i][j][k] = 1.0/d;
          if (setup->point_type[i][j-1][k] < INSIDE && (d = dist_from_contact(pt, tp[3], setup)) > 0)
            setup->dy[0][i][j][k] = 1.0/d;
          if (k < setup->numz-1 &&
              setup->point_type[i][j][k+1] < INSIDE && (d = dist_from_contact(pt, tp[4], setup)) > 0)
            setup->dz[1][i][j][k] = 1.0/d;
          if (setup->point_type[i][j][k-1] < INSIDE && (d = dist_from_contact(pt, tp[5], setup)) > 0)
            setup->dz[0][i][j][k] = 1.0/d;
        }
      }
    }
  }

  printf("\n");
  return setup->nseg_phi * setup->nseg_z + 1;
}

int init_ev_calc(GRETA_Siggen_Setup *setup) {
  int   i, j, k;
  float x, y, r;

  for (i = 0; i < setup->numx; i++) {
    x = setup->x0 + i*setup->xtal_grid;
    for (j = 0; j < setup->numy; j++) {
      y = setup->y0 + j*setup->xtal_grid;
      r = sqrt(x*x + y*y);
      for (k = 0; k < setup->numz; k++) {
        setup->v[0][i][j][k] = setup->v[1][i][j][k] = setup->xtal_HV/5.0;
	if (setup->point_type[i][j][k] == CONTACT_0) {
	  setup->v[0][i][j][k] = setup->v[1][i][j][k] = 0.0;
	} else if (setup->point_type[i][j][k] == CONTACT_VB) {
	  setup->v[0][i][j][k] = setup->v[1][i][j][k] = setup->xtal_HV;
	} else if (setup->point_type[i][j][k] == OUTSIDE) {
	  setup->v[0][i][j][k] = setup->v[1][i][j][k] = 0.0;
	} else {
	  //setup->v[0][i][j][k] = setup->v[1][i][j][k] = setup->xtal_HV*(setup->xtal_radius - r) /
          //  (Setup->xtal_radius - setup->core_radius);
	  if (setup->v[0][i][j][k] >= setup->xtal_HV*0.95) setup->v[0][i][j][k] = setup->v[1][i][j][k] = setup->xtal_HV*0.95;
	  if (setup->v[0][i][j][k] <= setup->xtal_HV*0.05) setup->v[0][i][j][k] = setup->v[1][i][j][k] = setup->xtal_HV*0.05;
	}
      }
    }
  }
  printf("ev ... init done\n");
  return 0;

}

#define HSG 0.25  // half of segmentation boundary gap (mm)

 int init_wp_calc(GRETA_Siggen_Setup *setup, int cnum) {
  int i, j, k;
  point pt;
  int   seg, seg2xy, seg2z;
  float dxy, dz;

  for (i = 0; i < setup->numx; i++) {
    for (j = 0; j < setup->numy; j++) {
      for (k = 0; k < setup->numz; k++) {
	setup->v[0][i][j][k] = setup->v[1][i][j][k] = 0.0;
      }
    }
  }

  for (i = 0; i < setup->numx; i++) {
    pt.x = setup->x0 + i*setup->xtal_grid;
    for (j = 0; j < setup->numy; j++) {
      pt.y = setup->y0 + j*setup->xtal_grid;
      for (k = 0; k < setup->numz; k++) {
	pt.z = setup->zmin + k*setup->xtal_grid;
	if (setup->point_type[i][j][k] == CONTACT_VB) {
	  if (cnum == setup->nseg_z*setup->nseg_phi) setup->v[1][i][j][k] = 1.0;
	  else setup->v[1][i][j][k] = 0.0;
	} else if (setup->point_type[i][j][k] == CONTACT_0) {
	  if (cnum < setup->nseg_z*setup->nseg_phi) {
	    // if (segment_number(setup, pt) == cnum) setup->v[1][i][j][k] = 1.0;
	    if ((seg = segment_number_new(setup, pt, &dxy, &dz, &seg2xy, &seg2z)) == cnum) {
	      setup->v[1][i][j][k] = 1.0;
              if (dxy < HSG) setup->v[1][i][j][k]  = 0.5 * (1.0 + dxy/HSG);
              if (dz  < HSG) setup->v[1][i][j][k] *= 0.5 * (1.0 + dz/HSG);
            } else if (seg2xy == cnum && dxy < HSG) {
              setup->v[1][i][j][k] = 0.5 * (1.0 - dxy/HSG);
            } else if (seg2z == cnum && dz < HSG) {
              setup->v[1][i][j][k] = 0.5 * (1.0 - dz/HSG);
            }
            if (cnum < 6 && fabs(pt.x) < HSG/2.0 && fabs(pt.y) < HSG/2.0)
              setup->v[1][i][j][k] = 0.167;
            else if (cnum < 6 && fabs(pt.x) < HSG*2.0 && fabs(pt.y) < HSG*2.0)
              setup->v[1][i][j][k] *= 0.33*(1.0 + sqrt(pt.x*pt.x + pt.y*pt.y)/HSG);
            // if (cnum < 6 && fabs(pt.x) < 4.1 && fabs(pt.y) < 4.1 && setup->v[1][i][j][k] > 0)
            //   printf("x, y = %5.2f %5.2f ; WP, dxy =  %5.2f %5.2f\n", pt.x, pt.y, setup->v[1][i][j][k], dxy);
	  } else {
	    setup->v[1][i][j][k] = 0.0;
	  }
	} else if (setup->point_type[i][j][k]  == OUTSIDE) {
	  setup->v[1][i][j][k] = 0.0;
	} else {
	  setup->v[1][i][j][k] = 0.5;
	}
        setup->v[0][i][j][k] = setup->v[1][i][j][k];
      }
    }
  }
  printf("wp ... init done\n");
  return 0;
}


/* segment_number_new()
   returns segment number
   seg2* = nearest neighbor segment in xy or z direction
   distance* = distance of projected point from center of nearest segment boundary in xy or z direction */
static int segment_number_new(GRETA_Siggen_Setup *setup, struct point pt,
                              float *distance_xy, float *distance_z, int *seg2_xy, int *seg2_z) {
  int i, seg = -1, segz, segz2;
  struct point c1, c2, corners[NCORNERS/2];
  float a, b, d;

  if (fabsf(pt.x) < 0.001 && fabsf(pt.y) < 0.001) return -1;

  /* find segment number in z direction */
  for (segz = 0; segz < setup->nseg_z; segz++)
    if (pt.z <= setup->zmax_segment[segz]) break;
  if (segz >= setup->nseg_z) return -1;

  /* find nearest-neighbor segz and distance to middle of z-segment line */
  if (segz == 0) {
    segz2 = segz + 1;
    *distance_z = setup->zmax_segment[segz] - pt.z;
  } else if (segz == setup->nseg_z-1) {
    segz2 = segz - 1;
    *distance_z = pt.z - setup->zmax_segment[segz - 1];
  } else {
    if (setup->zmax_segment[segz] - pt.z < pt.z - setup->zmax_segment[segz - 1]) {
      segz2 = segz + 1;
      *distance_z = setup->zmax_segment[segz] - pt.z;
    } else {
      segz2 = segz - 1;
      *distance_z = pt.z - setup->zmax_segment[segz - 1];
    }
  }

  project_to_edges(setup, pt, corners);
  for (i = 0; i < NCORNERS/2; i++) {
    // identify two corners c1, c2
    c1 = corners[i];
    if (i<5) c2 = corners[i+1];
    else c2 = corners[0];

    if (fabsf(pt.x) < 0.001) {
      if (fabsf(c1.x - c2.x) < 0.1) continue;
      b = c1.x / (c1.x - c2.x);
      a = (c1.y + b*(c2.y - c1.y)) / pt.y;
    } else {
      b = (c1.x*pt.y - c1.y*pt.x) / (pt.x*(c2.y - c1.y) - pt.y*(c2.x - c1.x));
      a = (c1.x + b*(c2.x - c1.x)) / pt.x;
    }
    if (b < 0 || b > 1 || a < 0) continue; // projected point does not intersect line between corners
    d = sqrt((c2.y - c1.y)*(c2.y - c1.y) + (c2.x - c1.x)*(c2.x - c1.x));
    if (b < 0.5) {
      seg = i;
      *seg2_xy = i+1;
      if (i == 5) *seg2_xy = 0;
      *distance_xy = (0.5 - b) * d;
    } else {
      *seg2_xy = i;
      seg = i+1;
      if (i == 5) seg = 0;
      *distance_xy = (b - 0.5) * d;
    }
    // projection.x = a*pt.x; projection.y = a*pt.x; projection.z = pt.z;
    // if (pt.z < 0.01 && a > 1.0)
    *distance_xy /= a; // scale distance from segment boundary
    /* if a < 1 then pt is outside detector */

    *seg2_xy = (*seg2_xy + 4) % 6 + 6 * segz;
    *seg2_z  = (seg + 4) % 6 + 6 * segz2;
    seg      = (seg + 4) % 6 + 6 * segz;
    return seg;
  }
  printf("segment_number_new error: x y z %5.2f %5.2f %5.2f -> no solution\n", pt.x, pt.y, pt.z);
  return -1;
}

/* segment_number
   returns the (geometrical) segment number at point pt, or -1 
   if outside crystal
*/
static int segment_number(GRETA_Siggen_Setup *setup, point pt) {
  int seg_phi, seg_z;
  float thp, thn, th, costh;
  struct point xyvect, edge_p[NCORNERS/2], side_pxy, cp, cpold;
  int i, j, seg_phi_n, seg_phi_p;

  //if (!in_crystal(setup, pt)) return -1;

  xyvect = pt;
  xyvect.z = 0.0;
  project_to_edges(setup, pt, edge_p);
  seg_z = 0;
  /* find segment number in z direction */
  for (i = 0; i < setup->nseg_z; i++) {
    if (pt.z <= setup->zmax_segment[i]) {
      seg_z = i;
      break;
    }
  }
  /* find segment number in phi direction */
  thp = thn = M_PI;
  side_pxy.z = 0;
  seg_phi = seg_phi_p = seg_phi_n = 0;
  for (i = 0; i < NCORNERS/2; i++) {
    j = (i+1)%6;
    side_pxy.x = (edge_p[i].x + edge_p[j].x)/2;
    side_pxy.y = (edge_p[i].y + edge_p[j].y)/2;
    costh = dot_prod(side_pxy, xyvect)
      /sqrt(dot_prod(side_pxy, side_pxy)*dot_prod(xyvect, xyvect));
    th = acos(costh);
    cp = cross_prod(side_pxy, xyvect);
    if (i == 0) {
      cpold = cp;
    }      
    if (cpold.z*cp.z >= 0) {
      if (th <= thp) {
	thp = th;
	seg_phi_p = i;
      }
    } else {
      if (th <= thn) {
	thn  = th;
	seg_phi_n = i;
      }
    }
  }
  if (seg_phi_p - seg_phi_n == 1) seg_phi = seg_phi_p;
  if (seg_phi_n - seg_phi_p == 1) seg_phi = seg_phi_n;
  if (abs(seg_phi_n - seg_phi_p) == 5) seg_phi = 0;
  //if (sign == 0.0) seg_phi = 0; //what?
  if (setup->xtal_type == CRYSTAL_A)
    seg_phi = (seg_phi + 4) % 6;
  if (seg_phi < 0) seg_phi += 6;
  
  return setup->nseg_phi*seg_z + seg_phi;
}

/* returns 0 (false) or 1 (true) depending on whether pt is inside the crystal */
static int in_crystal(GRETA_Siggen_Setup *setup, point pt) {
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
