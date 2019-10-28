#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

float fminf(float x, float y);

#include "greta_siggen.h"

static int   in_crystal(GRETA_Siggen_Setup *setup, point pt);
static int   init_point_fractions(GRETA_Siggen_Setup *setup);
static float get_fraction(GRETA_Siggen_Setup *setup, int i, int j, int k, int nsteps);
static int   project_to_edges(GRETA_Siggen_Setup *setup, struct point pt, struct point projections[NCORNERS/2]);
static int   segment_number(GRETA_Siggen_Setup *setup, point pt);


int geometry_init(GRETA_Siggen_Setup *setup) {
  FILE *fp;
  char line[MAX_LINE], *cp, *next;
  int i, j, k;
  struct point pt;

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

  init_point_fractions(setup);
  return setup->nseg_phi * setup->nseg_z + 1;
}

int init_ev_calc(GRETA_Siggen_Setup *setup) {
  int i, j, k;
  point pt;
  float r;

  for (i = 0; i < setup->numx; i++) {
    pt.x = setup->xmin + i*setup->xtal_grid;
    for (j = 0; j < setup->numy; j++) {
      pt.y = setup->ymin + j*setup->xtal_grid;
      r = sqrt(pt.x*pt.x + pt.y*pt.y);
      for (k = 0; k < setup->numz; k++) {
	pt.z = setup->zmin + k*setup->xtal_grid;
	if (setup->point_type[i][j][k] == CONTACT_0) {
	  setup->v[0][i][j][k] = setup->v[1][i][j][k] = 0.0;
	} else if (setup->point_type[i][j][k] == CONTACT_VB) {
	  setup->v[0][i][j][k] = setup->v[1][i][j][k] = setup->xtal_HV;
	} else if (setup->point_type[i][j][k]  == OUTSIDE) {
	  setup->v[0][i][j][k] = setup->v[1][i][j][k] = 0.0;
	} else {
	  setup->v[0][i][j][k] = setup->v[1][i][j][k] = setup->xtal_HV*(setup->xtal_radius - r) /
            (setup->xtal_radius - setup->core_radius);
	  if (setup->v[0][i][j][k] >= setup->xtal_HV*0.95) setup->v[0][i][j][k] = setup->xtal_HV*0.95;
	  if (setup->v[0][i][j][k] <= setup->xtal_HV*0.05) setup->v[0][i][j][k] = setup->xtal_HV*0.05;
	}
      }
    }
  }
  printf("ev ... init done\n"); fflush(stdout);
  return 0;

}

int init_wp_calc(GRETA_Siggen_Setup *setup, int cnum) {
  int i, j, k;
  point pt;
  float r;

  for (i = 0; i < setup->numx; i++) {
    for (j = 0; j < setup->numy; j++) {
      for (k = 0; k < setup->numz; k++) {
	setup->v[0][i][j][k] = setup->v[1][i][j][k] = 0.0;
      }
    }
  }

  for (i = 0; i < setup->numx; i++) {
    pt.x = setup->xmin + i*setup->xtal_grid;
    for (j = 0; j < setup->numy; j++) {
      pt.y = setup->ymin + j*setup->xtal_grid;
      r = sqrt(pt.x*pt.x + pt.y*pt.y);
      for (k = 0; k < setup->numz; k++) {
	pt.z = setup->zmin + k*setup->xtal_grid;
	if (setup->point_type[i][j][k] == CONTACT_VB) {
	  if (cnum == setup->nseg_z*setup->nseg_phi)
	    setup->v[0][i][j][k] = setup->v[1][i][j][k] = 1.0;
	  else
	    setup->v[0][i][j][k] = setup->v[1][i][j][k] = 0.0;
	} else if (setup->point_type[i][j][k] == CONTACT_0) {
	  if (cnum < setup->nseg_z*setup->nseg_phi) {
	    if (segment_number(setup, pt) == cnum)
	      setup->v[0][i][j][k] = setup->v[1][i][j][k] = 1.0;
	  } else {
	    setup->v[0][i][j][k] = setup->v[1][i][j][k] = 0.0;
	  }
	} else if (setup->point_type[i][j][k]  == OUTSIDE) {
	  setup->v[0][i][j][k] = setup->v[1][i][j][k] = 0.0;
	} else {
	  setup->v[0][i][j][k] = setup->v[1][i][j][k] = 0.5;
	}
      }
    }
  }
  printf("wp ... init done\n"); fflush(stdout);
  return 0;


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
  /*find segment number in z direction*/
  for (i = 0; i < setup->nseg_z; i++) {
    if (pt.z <= setup->zmax_segment[i]) {
      seg_z = i;
      break;
    }
  }
  /*find segment number in phi direction*/
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

#define SQ(x) ((x)*(x))

static int init_point_fractions(GRETA_Siggen_Setup *setup) {
  int i, j, k;
  point pt;
  float r, f, r2;

  printf("\n  Initializing points in crystal for x = %.1f to %.1f, grid %.2f...\n",
         setup->xmin, setup->xmax, setup->xtal_grid);
  for (i = 0; i < setup->numx; i++) {
    printf("\r %d/%d", i, setup->numx-1);fflush(stdout);
    pt.x = setup->xmin + i*setup->xtal_grid;
    for (j = 0; j < setup->numy; j++) {
      pt.y = setup->ymin + j*setup->xtal_grid;
      for (k = 0; k < setup->numz; k++) {
	pt.z = setup->zmin + k*setup->xtal_grid;
	f = get_fraction(setup, i, j, k, 3);
	if (f == 1.0) {                           // voxel is all inside bulk
	  setup->point_type[i][j][k] = INSIDE;
	} else if (f == 0.0) {                    // voxel is all outside
	  setup->point_type[i][j][k] = OUTSIDE;
	} else {                                  // voxel is on a boundary
	  r = sqrt(pt.x*pt.x + pt.y*pt.y);
	  if (pt.z < setup->xtal_grid) {            /* front face */
	    setup->point_type[i][j][k] = CONTACT_0;
	  } else if (r <= (setup->core_radius + setup->core_offset_max +
                           setup->bottom_taper_width + setup->xtal_grid) &&
                     k < setup->numz-1) {          /* core contact; see also below for k == setup->numz-1 */
	    setup->point_type[i][j][k] = CONTACT_VB;
	  } else if (pt.z >= setup->xtal_length - setup->xtal_grid &&
                     r <= setup->xtal_radius - setup->xtal_grid) {     /* passivated surface */
            setup->point_type[i][j][k] = INSIDE;
	  } else {                                                     /* outside contact */
	    setup->point_type[i][j][k] = CONTACT_0;
	  }
	}
      }
    }
  }
  for (i = 1; i < setup->numx-1; i++) {
    pt.x = setup->xmin + i*setup->xtal_grid;
    for (j = 1; j < setup->numy-1; j++) {
      pt.y = setup->ymin + j*setup->xtal_grid;
      k = setup->numz-1;
      /* top of central hole should be CONTACT_VB */
      r = sqrt(SQ(pt.x - setup->core_offset_x_top) + SQ(pt.y - setup->core_offset_y_top));
      r2 = setup->core_radius + setup->bottom_taper_width;
      if (r  < r2) setup->point_type[i][j][k] = OUTSIDE;
      if (r <= r2 && r + setup->xtal_grid > r2) setup->point_type[i][j][k] = CONTACT_VB;
      if (0 && setup->point_type[i][j][k] == CONTACT_VB)
        printf(" >> top of core:  %3d %3d %3d, %5.2f %5.2f %5.2f\n", i, j, k,
               setup->xmin + i*setup->xtal_grid,
               setup->ymin + j*setup->xtal_grid,
               setup->zmin + k*setup->xtal_grid);

      for (k = 1; k < setup->numz-1; k++) {
	if (setup->point_type[i][j][k] != INSIDE) continue;
	if ((i < setup->numx/2 && setup->point_type[i-1][j][k] == OUTSIDE) ||
            (i > setup->numx/2 && setup->point_type[i+1][j][k] == OUTSIDE) ||
            (j < setup->numy/2 && setup->point_type[i][j-1][k] == OUTSIDE) ||
            (j > setup->numy/2 && setup->point_type[i][j+1][k] == OUTSIDE)) {
	  setup->point_type[i][j+1][k] = CONTACT_0;
	}
      }
    }
  }

  printf("\n");
  return 0;
}

static float get_fraction(GRETA_Siggen_Setup *setup, int i, int j, int k, int nsteps) {
  int ii, jj, kk;
  point pt;
  int n, m;

  // first try all eight corners of voxel to see if they are inside the detector
  n = m = 0;
  for (ii = 0; ii < 2; ii++) {
    pt.x = setup->xmin + (i + ii - 0.5)*setup->xtal_grid;
    if (pt.x < setup->xmin || pt.x > setup->xmax) {
      n += 4;
      continue;
    }
    for (jj = 0; jj < 2; jj++) {
      pt.y = setup->ymin + (j + jj - 0.5)*setup->xtal_grid;
      if (pt.y < setup->ymin || pt.y > setup->ymax) {
	n += 2;
	continue;
      }
      for (kk = 0; kk < 2; kk++) {
        n++;
	pt.z = setup->zmin + (k + kk - 0.5)*setup->xtal_grid;
	if (pt.z < setup->zmin || pt.z > setup->zmax) continue;
	if (in_crystal(setup, pt)) m++;
      }
    }
  }
  if (m == 0) return 0.0;  // all are outside
  if (m == n) return 1.0;  // all are inside

  /*
  // no definitive answer, so use (nsteps+1)^3 smaller fractions of the voxel
  n = m = 0;
  for (ii = -nsteps; ii <= nsteps; ii++) {
    pt.x = setup->xmin + (i + ii/2.0/nsteps)*setup->xtal_grid;
    if (pt.x < setup->xmin || pt.x > setup->xmax) {
      n += (2*nsteps+1) * (2*nsteps+1);
      continue;
    }
    for (jj = -nsteps; jj <=nsteps; jj++) {
      pt.y = setup->ymin + (j + jj/2.0/nsteps)*setup->xtal_grid;
      if (pt.y < setup->ymin || pt.y > setup->ymax) {
	n += 2*nsteps+1;
	continue;
      }
      for (kk = -nsteps; kk <= nsteps; kk++) {
        n++;
	pt.z = setup->zmin + (k + kk/2.0/nsteps)*setup->xtal_grid;
	if (pt.z < setup->zmin || pt.z > setup->zmax) continue;
	if (in_crystal(setup, pt)) m++;
      }
    }
  }
  */
  return ((float) m) / ((float) n);
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
