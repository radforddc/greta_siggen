/* program to calculate weighting potentials of coaxial Ge detectors
   by relaxation

   to run: ./fieldgen3d <config_file_name> -b bias_volts -w {0,1} -d {0,1} -p {0,1})
                        -z rho_spectrum_file_name

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
float fminf(float x, float y);
#include <unistd.h>
#include <ctype.h>

#include "greta_siggen.h"
#include "point.h"

#define e_over_E (1.413*8/3)

static int grid_init(GRETA_Siggen_Setup *setup);
static int ev_calc(GRETA_Siggen_Setup *setup, GRETA_Siggen_Setup *old_setup);
static int wp_calc(GRETA_Siggen_Setup *setup, GRETA_Siggen_Setup *old_setup, int cnum);
static int malloc_arrays(GRETA_Siggen_Setup *setup);
static int write_ev(GRETA_Siggen_Setup *setup, char*fname);
static int write_wp(GRETA_Siggen_Setup *setup, char*fname);
static int fname_insert(char fname[256], char dest[256], char *str);
int do_relaxation(GRETA_Siggen_Setup *setup, int efld_calc);

int main(int argc, char *argv[]) {
  GRETA_Siggen_Setup setup1, setup2, setup3, *setup = &setup1;
  float rho_z[256] = {0};
  int   i;
  FILE  *fp;


  if (argc%2 != 0 || !strncmp(argv[1], "-h", 2)) {
    printf("usage: %s <config_file_name>\n"
           "        -b bias_volts\n"
           "        -w {0,1}   (do_not/do write the field file)\n"
           "        -p {0,1}   (do_not/do write the WP file)\n"
           "        -f {0,1}   (do_not/do try to fix potential values from larger grid sizes)\n"
           "                    - note that this reduces accuracy of min field and depl voltage\n"
           "        -z rho_spectrum_file_name\n", argv[0]);
    return 0;
  }
  if (read_config(argv[1], &setup1)) return 1;
  strncpy(setup1.config_file_name, argv[1], 256);
  setup1.fix_adaptive = 0;

  for (i=2; i<argc-1; i+=2) {
    if (strstr(argv[i], "-b")) {
      setup1.xtal_HV = atof(argv[i+1]);     // bias volts
    } else if (strstr(argv[i], "-w")) {
      setup1.write_field = atoi(argv[i+1]); // write-out options
    } else if (strstr(argv[i], "-p")) {
      setup1.write_WP = atoi(argv[i+1]);    // weighting-potential options
    } else if (strstr(argv[i], "-f")) {
      setup1.fix_adaptive = atoi(argv[i+1]);      // adaptive-grid-potential fixing option
    } else if (strstr(argv[i], "-r")) {
      if (!(fp = fopen(argv[i+1], "r"))) { // impurity-profile-spectrum file name
        printf("\nERROR: cannot open impurity profile spectrum file %s\n\n", argv[i+1]);
        return 1;
      }
      fread(rho_z, 36, 1, fp);
      fread(rho_z, sizeof(rho_z), 1, fp);
      fclose(fp);
      printf(" z(mm)   rho\n");
      for (i=0; i < 200 && rho_z[i] != 0.0f; i++)  printf(" %3d  %7.3f\n", i, rho_z[i]);
    } else {
      printf("Possible options:\n"
	     "      -b bias_volts\n"
	     "      -w {0,1,2} (for WV options)\n"
	     "      -p {0,1}   (for WP options)\n"
             "      -f {0,1}   (do_not/do try to fix potential values from larger grid sizes)\n"
             "                  - note that this reduces accuracy of min field and depl voltage\n"
             "      -r rho_spectrum_file_name\n");
      return 1;
    }
    // FIXME - implement other option: -z
  }
  printf("write_WP = %d\n", setup1.write_WP);

  /* modify core dimensions by Li thickness */
  setup1.core_radius += setup1.Li_thickness;
  setup1.core_length += setup1.Li_thickness;
  setup1.core_gap    -= setup1.Li_thickness;
  setup1.core_bullet_radius  += setup1.Li_thickness;
  setup1.core_bullet_radius  += setup1.Li_thickness;
  setup1.bottom_taper_width  += setup1.Li_thickness * 0.7;
  setup1.bottom_taper_length += setup1.Li_thickness * 0.7;
  if (setup1.Li_thickness > 0.1)
    printf("*** Note: *** Increased core contact dimensions by Li_thickness = %.2f mm\n",
           setup1.Li_thickness);

  printf("number of iterations: %d\n", setup1.max_iterations);

  memcpy(&setup2, &setup1, sizeof(setup1));
  memcpy(&setup3, &setup1, sizeof(setup1));

  setup1.xtal_grid *= 9.0;
  setup2.xtal_grid *= 3.0;
  if (grid_init(&setup1) != 0 ||
      grid_init(&setup2) != 0 ||
      grid_init(&setup3) != 0) {
    error("failed to init field calculations\n");
    return 1;
  }
  if (setup->write_field) {
    setup1.write_field = 0; // no need to save intermediate calculations
    setup2.write_field = 0;
    if (setup3.xtal_grid > 0.4) {
      ev_calc(&setup2, NULL);
    } else {
      ev_calc(&setup1, NULL);
      ev_calc(&setup2, &setup1);
    }
    ev_calc(&setup3, &setup2);
  }
  if (setup->write_WP) {
    setup1.write_WP = 0; // no need to save intermediate calculations
    setup2.write_WP = 0;
    int cnum;
    for (cnum = 0; cnum < setup->ncontacts; cnum++) {
      if (setup3.xtal_grid > 0.4) {
        wp_calc(&setup2, NULL, cnum);
      } else {
        wp_calc(&setup1, NULL, cnum);
        wp_calc(&setup2, &setup1, cnum);
      }
      wp_calc(&setup3, &setup2, cnum);
    }
  }
 
  return 0;
}


static int grid_init(GRETA_Siggen_Setup *setup) {

  setup->xmin = setup->ymin = - setup->xtal_radius;
  setup->xmax = setup->ymax = setup->xtal_radius;
  setup->zmin = 0; setup->zmax = setup->xtal_length;

  setup->numx = 2*(int) (setup->xmax/setup->xtal_grid + 1.999) - 1;
  setup->numy = 2*(int) (setup->ymax/setup->xtal_grid + 1.999) - 1;
  setup->numz =   (int) (setup->zmax/setup->xtal_grid + 1.999);
  printf("numx, numy, numz: %d %d %d\n", setup->numx, setup->numy, setup->numz);
  setup->x0 = -(setup->numx/2) * setup->xtal_grid;
  setup->y0 = -(setup->numy/2) * setup->xtal_grid;

  if (malloc_arrays(setup) < 0) return -1;

  if ((setup->ncontacts = geometry_init(setup)) < 0) {
    error("failed to init geometry\n");
    return 1;
  }

  if (setup->impurity_lamda != 0 || setup->rho_b != 0) {
    setup->rho_c = -setup->impurity_gradient/10 -
                    setup->rho_b / setup->xtal_length * (1 - exp(-setup->xtal_length / setup->impurity_lamda));
    printf("rho_c = %f\n", setup->rho_c);
  }

  printf("\n\n"
	 "          grid size: %5.2f\n"
	 "         xmin, xmax: %5.2f %5.2f\n"
	 "         ymin, ymax: %5.2f %5.2f\n"
	 "         zmin, zmax: %5.2f %5.2f\n"
	 "     HV, rho0, drho: %.0f %.4f %.4f\n"
	 "lamda, rho_b, rho_c: %.4f %.4f %.4f\n", 
	 setup->xtal_grid, setup->xmin, setup->xmax, setup->ymin, setup->ymax,
         setup->zmin, setup->zmax, setup->xtal_HV, setup->impurity_z0, setup->impurity_gradient,
         setup->impurity_lamda, setup->rho_b, setup->rho_c);
  printf("\n");

  return 0;
}


int do_relaxation(GRETA_Siggen_Setup *setup, int efld_calc) {
  int    i, j, k, *ylo, *yhi, **zlo;
  int    old, new, iter, undep;
  double sum_dif, max_dif, dif, save_dif, mean;
  double ****v = setup->v;
  float  z, ****dx = setup->dx, ****dy = setup->dy, ****dz = setup->dz;
  double OR_fact, d1, d2;
  int    L  = lrint(setup->xtal_length/setup->xtal_grid);
  int    R  = lrint(setup->xtal_radius/setup->xtal_grid);


  if ((ylo = malloc(setup->numx*sizeof(*ylo))) == NULL ||
      (yhi = malloc(setup->numx*sizeof(*yhi))) == NULL ||
      (zlo = malloc(setup->numx*sizeof(*zlo))) == NULL) {
    error("malloc failed\n");
    return -1;
  }
  for (i = 0; i < setup->numx; i++) {
    if ((zlo[i] = malloc(setup->numy*sizeof(*zlo[i]))) == NULL) {
      error("malloc failed\n");
      return -1;
    }
  }
  printf("grid, L, R = %f %d %d\n", setup->xtal_grid, L, R);

  old = 1; new = 0;
  for (i = 0; i < setup->numx; i++) {
    ylo[i] = setup->numy-1;
    yhi[i] = 1;
    for (j = 0; j < setup->numy; j++) {
      zlo[i][j] = setup->numz;
      for (k = 0; k < setup->numz; k++) {
	if (setup->point_type[i][j][k] == OUTSIDE ||
            setup->point_type[i][j][k] == CONTACT_0) continue;
        if (ylo[i] > j) ylo[i] = j;
        if (yhi[i] < j) yhi[i] = j;
        if (zlo[i][j] > k) {
          zlo[i][j] = k;
        }
	if (setup->point_type[i][j][k] != INSIDE &&
            setup->point_type[i][j][k] != FIXED &&
            setup->point_type[i][j][k] != CONTACT_VB &&
            setup->point_type[i][j][k] != CONTACT_EDGE &&
	    (i > 0 && i + 1 < setup->numx) &&
            (j > 0 && j + 1 < setup->numy) &&
            (k > 0 && k < setup->numz)) {  // k = numz is passivated surface
          printf("error for point %d %d %d : %d\n", i, j, k, setup->point_type[i][j][k]);
	}
      }
    }
    //printf("i, ylo[i], yhi[i] = %3d %3d %3d\n", i, ylo[i], yhi[i]); fflush(stdout);
  }
  //printf("grid, L, R = %f %d %d\n", setup->xtal_grid, L, R); fflush(stdout);

  /* calculate voxel impurity value as a function of z */
  for (k = 0; k < setup->numz; k++) {
    z = setup->zmin + k*setup->xtal_grid;
    if (setup->impurity_lamda == 0.0 || setup->rho_b == 0.0) {
      setup->voxel_impurity_z[k] = (setup->impurity_z0 + setup->impurity_gradient * z/10) *
        3.0 * setup->xtal_grid*setup->xtal_grid  * e_over_E/6.0;
    } else {
      setup->voxel_impurity_z[k] = (setup->impurity_z0 -
                                    setup->rho_b*(1-exp(-z/setup->impurity_lamda)) -
                                    setup->rho_c*z) *
        3.0 * setup->xtal_grid*setup->xtal_grid * e_over_E/6.0;
    }
  }
  /* add surface charge to passivated surface */
  setup->voxel_impurity_z[setup->numz-1] +=
    setup->impurity_surface * 3.0 * setup->xtal_grid * e_over_E/6.0;

  printf("starting relaxation...\n");
  for (iter = 0; iter < setup->max_iterations; iter++) {

#define OVER_RELAX_F2   (efld_calc ? (0.88 + 0.02/setup->xtal_grid) : (0.88 + 0.02/setup->xtal_grid))

    if (efld_calc) {
      OR_fact = (0.964  - 236.0/(L*R)) * (1.0 - 3.0/(5.0+iter));
      // OR_fact = OVER_RELAX_F2;
      d1 = 0.6;
    } else {
      OR_fact = (0.983 - 200.0/(L*R)) * (1.0 - 3.0/(5.0+iter));
      d1 = 0.1;
    }
    if (240.0/(L*R) > 0.5) OR_fact = (0.5 * (1.0 - 0.5/(double)(1+iter/6)));
    if (iter < 3)  OR_fact = 0;

    d2 = 1.0 - d1;
    old = 1 - old;
    new = 1 - new;
    sum_dif = 0.0;
    max_dif = 0.0;

    /* do relaxation iteration */
    i = j = setup->numx/2;
    for (i = 1; i < setup->numx-1; i++) {
      for (j = 1; j < setup->numy-1; j++) {
      //for (j = ylo[i]; j <= yhi[i]; j++) {
        /* set up reflection symmetry around z=numz-1 (passivated surface) */
        v[old][i][j][setup->numz] = v[old][i][j][setup->numz-2];
        for (k = 1; k < setup->numz; k++) {
        //for (k = zlo[i][j]; k < setup->numz; k++) {
	  if (setup->point_type[i][j][k] < INSIDE) continue;
          // save step difference from previous iteration
          save_dif = v[old][i][j][k] - v[new][i][j][k];
          setup->e[i][j][k] = d1*setup->e[i][j][k] + d2*save_dif;

          if (setup->point_type[i][j][k] == INSIDE) {
            mean = (v[old][i-1][j][k] + v[old][i+1][j][k] +
                    v[old][i][j-1][k] + v[old][i][j+1][k] +
                    v[old][i][j][k-1] + v[old][i][j][k+1])/6.0;
          } else {
            mean = (v[old][i-1][j][k] * dx[0][i][j][k] +
                    v[old][i+1][j][k] * dx[1][i][j][k] +
                    v[old][i][j-1][k] * dy[0][i][j][k] +
                    v[old][i][j+1][k] * dy[1][i][j][k] +
                    v[old][i][j][k-1] * dz[0][i][j][k] +
                    v[old][i][j][k+1] * dz[1][i][j][k]) /
              (double) (dx[0][i][j][k] + dx[1][i][j][k] +
                        dy[0][i][j][k] + dy[1][i][j][k] +
                        dz[0][i][j][k] + dz[1][i][j][k]);
          }

          if (efld_calc) {
            v[new][i][j][k] = mean;// + setup->voxel_impurity_z[k];
          } else {
            v[new][i][j][k] = mean;
          }
	  dif = v[old][i][j][k] - v[new][i][j][k];
 
          //v[new][i][j][k] += OR_fact * save_dif;
          v[new][i][j][k] += OR_fact * setup->e[i][j][k];   // do over-relaxation

	  if (dif < 0.0) dif = -dif;
	  sum_dif += dif;
	  if (max_dif < dif) max_dif = dif;
	}
      }
    }
    sum_dif /= setup->numx * setup->numy * setup->numz;
    if (iter %100 == 0 ||
        (setup->xtal_grid < 0.4 && iter %50 == 0)) {
      printf("%5d %.3e %.3e\n", iter, max_dif, sum_dif);
      fflush(stdout);
    }
    if (max_dif < 1.0e-8) break;

    if (efld_calc) {
      //if (max_dif/setup->xtal_HV < 1.0e-3/5000.0 &&
      //    sum_dif/setup->xtal_HV < 6.0e-5/5000.0) break;
      if (max_dif/setup->xtal_HV < 2.0e-3/5000.0 &&
          sum_dif/setup->xtal_HV < 3.0e-4/5000.0) break;
   } else {
      if (max_dif < 5.0e-7 &&
          sum_dif < 1.0e-7) break;
    }
  }
  printf("%5d %.3e %.3e\n", 
	 iter, max_dif, sum_dif);

  undep = 0;
  for (i = 0; i < setup->numx && !undep; i++) {
    for (j = 0; j < setup->numy && !undep; j++) {
      for (k = 0; k < setup->numz && !undep; k++) {
	if (v[new][i][j][k] > setup->xtal_HV) {
	  printf("detector is undepleted!\n");
	  undep = 1;
	}
      }
    }
  }
  printf("...  relaxation done\n");
  for (i = 0; i < setup->numx; i++) free(zlo[i]);
  free(ylo); free(yhi); free(zlo);

  return 0;
}

static int v_interpolate(GRETA_Siggen_Setup *setup, GRETA_Siggen_Setup *old_setup) {
  int   i, j, k, n, i2, j2, k2;
  int   xmin, xmax, ymin, ymax, zmin, zmax;
  float f, fx, fy, fz;
  double ***ov = old_setup->v[1];

  /* the previous calculation was on a coarser grid...
     now copy/expand the potential to the new finer grid
  */
  n = (int) (old_setup->xtal_grid / setup->xtal_grid + 0.5);
  f = 1.0 / (float) n;
  printf("\ngrid %.4f -> %.4f; ratio = %d %.3f\n\n",
         old_setup->xtal_grid, setup->xtal_grid, n, f);

  for (i = i2 = 1; i < old_setup->numx-1; i++) {
    xmin = i2;
    xmax = i2 + n;
    if (xmax > setup->numx-1) xmax = setup->numx-1;
    for (j = j2 = 1; j < old_setup->numy-1; j++) {
      ymin = j2;
      ymax = j2 + n;
      if (ymax > setup->numy-1) ymax = setup->numy-1;
      for (k = k2 = 1; k < old_setup->numz-1; k++) {
        zmin = k2;
        zmax = k2 + n;
        if (zmax > setup->numz-1) zmax = setup->numz-1;
        fx = 1.0;
        for (i2 = xmin; i2 < xmax; i2++) {
          fy = 1.0;
          for (j2 = ymin; j2 < ymax; j2++) {
            fz = 1.0;
            for (k2 = zmin; k2 < zmax; k2++) {
              if (setup->point_type[i2][j2][k2] >= INSIDE)
                setup->v[0][i2][j2][k2] = setup->v[1][i2][j2][k2] =
                  fx       * fy       * fz       * ov[i  ][j  ][k  ] +
                  (1.0-fx) * fy       * fz       * ov[i+1][j  ][k  ] +
                  fx       * (1.0-fy) * fz       * ov[i  ][j+1][k  ] +
                  (1.0-fx) * (1.0-fy) * fz       * ov[i+1][j+1][k  ] +
                  fx       * fy       * (1.0-fz) * ov[i  ][j  ][k+1] +
                  (1.0-fx) * fy       * (1.0-fz) * ov[i+1][j  ][k+1] +
                  fx       * (1.0-fy) * (1.0-fz) * ov[i  ][j+1][k+1] +
                  (1.0-fx) * (1.0-fy) * (1.0-fz) * ov[i+1][j+1][k+1];
              fz -= f;
            }
            fy -= f;
          }
          fx -= f;
        }
        k2 = zmax;
      }
      j2 = ymax;
    }
    i2 = xmax;
  }

  return 0;
}

static int ev_calc(GRETA_Siggen_Setup *setup, GRETA_Siggen_Setup *old_setup) {

  if (!old_setup) printf("\n\n ---- starting EV calculation --- \n");
  init_ev_calc(setup);
  if (old_setup) v_interpolate(setup, old_setup);
  do_relaxation(setup, 1);
  if (setup->write_field) {
    printf("saving E, V in file %s\n", setup->field_name);
    write_ev(setup, setup->field_name);
  }
  return 0;
}

static int wp_calc(GRETA_Siggen_Setup *setup, GRETA_Siggen_Setup *old_setup, int cnum) {
  char fname[256], str[256];

  if (!old_setup)
    printf("\n\n ---- starting WP calculation for contact number %d --- \n", cnum);
  init_wp_calc(setup, cnum);
  if (old_setup) v_interpolate(setup, old_setup);
  do_relaxation(setup, 0);
  if (setup->write_WP) {
    snprintf(str, 256, "-seg%2.2d", cnum);
    fname_insert(setup->wp_name, fname, str);
    printf("saving WP in file %s\n", fname);
    write_wp(setup, fname);
  }
  return 0;
}

int report_config(FILE *fp_out, char *config_file_name) {

  char  *c, line[256];
  FILE  *file;

  fprintf(fp_out, "# Config file: %s\n", config_file_name);
  if (!(file = fopen(config_file_name, "r"))) return 1;

  while (fgets(line, sizeof(line), file)) {
    if (strlen(line) < 3 || *line == ' ' || *line == '\t' || *line == '#') continue;
    if ((c = strchr(line, '#')) || (c = strchr(line, '\n'))) *c = '\0';
    fprintf(fp_out, "# %s\n", line);
  }
  fclose(file);
  return 0;
}

int write_ev(GRETA_Siggen_Setup *setup, char *fname) {
  int    i,j,k;
  FILE   *fp = NULL, *fp2 = NULL;
  struct efld{float x; float y; float z;} ***e;
  float  grid = setup->xtal_grid;
  double ***v = setup->v[0];


  if ((fp = fopen(fname, "w")) == NULL) {
    error("Failed to open output file: %s\n", fname);
    return -1;
  }
  if ((fp2 = fopen("fields/ev2d.dat", "w")) == NULL) {
    error("Failed to open output file fields/ev2d.dat\n");
    return -1;
  }

  /* copy configuration parameters to output file */
  report_config(fp, setup->config_file_name);
  fprintf(fp, "#\n# HV bias in fieldgen: %.1f V\n", setup->xtal_HV);

  if ((e = malloc(setup->numx*sizeof(*e))) == NULL) {
    error("malloc failed\n");
    fclose(fp);
    return -1;
  }
  for (i = 0; i < setup->numx; i++) {
    if ((e[i] = malloc(setup->numy*sizeof(*e[i]))) == NULL) {
      error("malloc failed\n");
      fclose(fp);
      return -1;
    }
    for (j = 0; j < setup->numy; j++) {
      if ((e[i][j] = calloc(setup->numz, sizeof(*e[i][j]))) == NULL) {
	error("malloc failed\n");
	fclose(fp);
	return -1;
      }
    }
  }
  for (i = 1; i < setup->numx-1; i++) {
    for (j = 1; j < setup->numy-1; j++) {
      // E_z for k = 0
      e[i][j][0].z = 10.0*(v[i][j][0] - v[i][j][1])/grid;
      
      for (k = 1; k < setup->numz; k++) {
        if (setup->point_type[i][j][k] < INSIDE) continue;

        // calculate E_x
        if (setup->point_type[i][j][k] == CONTACT_EDGE) {
          e[i][j][k].x = 5.0*((v[i][j][k] - v[i+1][j][k])*setup->dx[1][i][j][k] +
                              (v[i-1][j][k] - v[i][j][k])*setup->dx[0][i][j][k])/grid;
        } else if (setup->point_type[i][j][k] == INSIDE) {
          e[i][j][k].x = 5.0*(v[i-1][j][k] - v[i+1][j][k])/grid;
        }

        // calculate E_y
        if (setup->point_type[i][j][k] == CONTACT_EDGE) {
          e[i][j][k].y = 5.0*((v[i][j][k] - v[i][j+1][k])*setup->dy[1][i][j][k] +
                              (v[i][j-1][k] - v[i][j][k])*setup->dy[0][i][j][k])/grid;
        } else if (setup->point_type[i][j][k] == INSIDE) {
          e[i][j][k].y = 5.0*(v[i][j-1][k] - v[i][j+1][k])/grid;
        }

        if (k == setup->numz-1) {
          // E_z for k = numz-1 (passivated surface)
          if (setup->point_type[i][j][k] == CONTACT_EDGE) {
            e[i][j][k].z = 10.0*(v[i][j][k-1] - v[i][j][k])*setup->dz[0][i][j][k]/grid;
          } else {
            e[i][j][k].z = 10.0*(v[i][j][k-1] - v[i][j][k])/grid;
          }
        } else {
          if (setup->point_type[i][j][k] == CONTACT_EDGE) {
            e[i][j][k].z = 5.0*((v[i][j][k] - v[i][j][k+1])*setup->dz[1][i][j][k] +
                                (v[i][j][k-1] - v[i][j][k])*setup->dz[0][i][j][k])/grid;
          } else if (setup->point_type[i][j][k] == INSIDE) {
            e[i][j][k].z = 5.0*(v[i][j][k-1] - v[i][j][k+1])/grid;
          }
        }
      }
    }
  }
  /*
  // calculate field AT the contact surfaces
  for (i = 1; i < setup->numx-1; i++) {
    for (j = 1; j < setup->numy-1; j++) {
      for (k = 1; k < setup->numz; k++) {
        if (setup->point_type[i][j][k] != CONTACT_0 &&
            setup->point_type[i][j][k] != CONTACT_VB) continue;

        if (setup->point_type[i+1][j][k] == CONTACT_EDGE)
          e[i][j][k].x = 10.0*(v[i][j][k] - v[i+1][j][k])*setup->dx[0][i+1][j][k]/grid;
        if (setup->point_type[i-1][j][k] == CONTACT_EDGE)
          e[i][j][k].x = 10.0*(v[i-1][j][k] - v[i][j][k])*setup->dx[1][i-1][j][k]/grid;
 
        if (setup->point_type[i][j+1][k] == CONTACT_EDGE)
          e[i][j][k].y = 10.0*(v[i][j][k] - v[i][j+1][k])*setup->dy[0][i][j+1][k]/grid;
        if (setup->point_type[i][j-1][k] == CONTACT_EDGE)
          e[i][j][k].y = 10.0*(v[i][j-1][k] - v[i][j][k])*setup->dy[1][i][j-1][k]/grid;

        if (setup->point_type[i][j][k+1] == CONTACT_EDGE)
          e[i][j][k].z = 10.0*(v[i][j][k] - v[i][j][k+1])*setup->dz[0][i][j][k+1]/grid;
        if (setup->point_type[i][j][k-1] == CONTACT_EDGE)
          e[i][j][k].z = 10.0*(v[i][j][k-1] - v[i][j][k])*setup->dz[1][i][j][k-1]/grid;
      }
    }
  }
  */
  /*
  // calculate field AT the contact surfaces by extrapolation
  for (i = 1; i < setup->numx-1; i++) {
    for (j = 1; j < setup->numy-1; j++) {
      for (k = 1; k < setup->numz; k++) {
        if (setup->point_type[i][j][k] != CONTACT_0 &&
            setup->point_type[i][j][k] != CONTACT_VB) continue;
        float f = 0;
        
        if (setup->point_type[i+1][j][k] == CONTACT_EDGE) {
          e[i][j][k].x += 2.0*e[i+1][j][k].x - e[i+2][j][k].x;
          e[i][j][k].y += 2.0*e[i+1][j][k].y - e[i+2][j][k].y;
          e[i][j][k].z += 2.0*e[i+1][j][k].z - e[i+2][j][k].z;
          f++;
        } else if (setup->point_type[i-1][j][k] == CONTACT_EDGE) {
          e[i][j][k].x += 2.0*e[i-1][j][k].x - e[i-2][j][k].x;
          e[i][j][k].y += 2.0*e[i-1][j][k].y - e[i-2][j][k].y;
          e[i][j][k].z += 2.0*e[i-1][j][k].z - e[i-2][j][k].z;
          f++;
        }

        if (setup->point_type[i][j+1][k] == CONTACT_EDGE) {
          e[i][j][k].x += 2.0*e[i][j+1][k].x - e[i][j+2][k].x;
          e[i][j][k].y += 2.0*e[i][j+1][k].y - e[i][j+2][k].y;
          e[i][j][k].z += 2.0*e[i][j+1][k].z - e[i][j+2][k].z;
          f++;
        } else if (setup->point_type[i][j-1][k] == CONTACT_EDGE) {
          e[i][j][k].x += 2.0*e[i][j-1][k].x - e[i][j-2][k].x;
          e[i][j][k].y += 2.0*e[i][j-1][k].y - e[i][j-2][k].y;
          e[i][j][k].z += 2.0*e[i][j-1][k].z - e[i][j-2][k].z;
          f++;
        }

        if (setup->point_type[i][j][k+1] == CONTACT_EDGE) {
          e[i][j][k].x += 2.0*e[i][j][k+1].x - e[i][j][k+2].x;
          e[i][j][k].y += 2.0*e[i][j][k+1].y - e[i][j][k+2].y;
          e[i][j][k].z += 2.0*e[i][j][k+1].z - e[i][j][k+2].z;
          f++;
        } else if (setup->point_type[i][j][k-1] == CONTACT_EDGE) {
          e[i][j][k].x += 2.0*e[i][j][k-1].x - e[i][j][k-2].x;
          e[i][j][k].y += 2.0*e[i][j][k-1].y - e[i][j][k-2].y;
          e[i][j][k].z += 2.0*e[i][j][k-1].z - e[i][j][k-2].z;
          f++;
        }
        if (f > 1.0) {
          e[i][j][k].x /= f;
          e[i][j][k].y /= f;
          e[i][j][k].z /= f;
        }
      }
    }
  }
  */
  // calculate field AT the contact surfaces: take mean of neighbors
  for (i = 1; i < setup->numx-1; i++) {
    for (j = 1; j < setup->numy-1; j++) {
      for (k = 1; k < setup->numz; k++) {
        if (setup->point_type[i][j][k] != CONTACT_0 &&
            setup->point_type[i][j][k] != CONTACT_VB) continue;
        float f = 0;
        
        if (setup->point_type[i+1][j][k] == CONTACT_EDGE) {
          e[i][j][k].x += e[i+1][j][k].x;
          e[i][j][k].y += e[i+1][j][k].y;
          e[i][j][k].z += e[i+1][j][k].z;
          f++;
        } else if (setup->point_type[i-1][j][k] == CONTACT_EDGE) {
          e[i][j][k].x += e[i-1][j][k].x;
          e[i][j][k].y += e[i-1][j][k].y;
          e[i][j][k].z += e[i-1][j][k].z;
          f++;
        }

        if (setup->point_type[i][j+1][k] == CONTACT_EDGE) {
          e[i][j][k].x += e[i][j+1][k].x;
          e[i][j][k].y += e[i][j+1][k].y;
          e[i][j][k].z += e[i][j+1][k].z;
          f++;
        } else if (setup->point_type[i][j-1][k] == CONTACT_EDGE) {
          e[i][j][k].x += e[i][j-1][k].x;
          e[i][j][k].y += e[i][j-1][k].y;
          e[i][j][k].z += e[i][j-1][k].z;
          f++;
        }

        if (setup->point_type[i][j][k+1] == CONTACT_EDGE) {
          e[i][j][k].x += e[i][j][k+1].x;
          e[i][j][k].y += e[i][j][k+1].y;
          e[i][j][k].z += e[i][j][k+1].z;
          f++;
        } else if (setup->point_type[i][j][k-1] == CONTACT_EDGE) {
          e[i][j][k].x += e[i][j][k-1].x;
          e[i][j][k].y += e[i][j][k-1].y;
          e[i][j][k].z += e[i][j][k-1].z;
          f++;
        }
        if (f > 1.0) {
          e[i][j][k].x /= f;
          e[i][j][k].y /= f;
          e[i][j][k].z /= f;
        }
      }
    }
  }
  /*
  // calculate field AT the contact surfaces: take field from nearest bulk point
  for (i = 1; i < setup->numx-1; i++) {
    for (j = 1; j < setup->numy-1; j++) {
      for (k = 1; k < setup->numz; k++) {
        if (setup->point_type[i][j][k] != CONTACT_0 &&
            setup->point_type[i][j][k] != CONTACT_VB) continue;
        float d = 0.5/grid;
        
        if (setup->point_type[i+1][j][k] == CONTACT_EDGE) {
          e[i][j][k].x = e[i+1][j][k].x;
          e[i][j][k].y = e[i+1][j][k].y;
          e[i][j][k].z = e[i+1][j][k].z;
          d = setup->dx[0][i+1][j][k];
        } else if (setup->point_type[i-1][j][k] == CONTACT_EDGE) {
          e[i][j][k].x = e[i-1][j][k].x;
          e[i][j][k].y = e[i-1][j][k].y;
          e[i][j][k].z = e[i-1][j][k].z;
          d = setup->dx[1][i-1][j][k];
        }

        if (setup->point_type[i][j+1][k] == CONTACT_EDGE &&
            d < setup->dy[0][i][j+1][k]) {
          e[i][j][k].x = e[i][j+1][k].x;
          e[i][j][k].y = e[i][j+1][k].y;
          e[i][j][k].z = e[i][j+1][k].z;
          d = setup->dy[0][i][j+1][k];
        } else if (setup->point_type[i][j-1][k] == CONTACT_EDGE &&
                   d < setup->dy[1][i][j-1][k]) {
          e[i][j][k].x = e[i][j-1][k].x;
          e[i][j][k].y = e[i][j-1][k].y;
          e[i][j][k].z = e[i][j-1][k].z;
          d = setup->dy[1][i][j-1][k];
        }

        if (setup->point_type[i][j][k+1] == CONTACT_EDGE &&
            d < setup->dz[0][i][j][k+1]) {
          e[i][j][k].x = e[i][j][k+1].x;
          e[i][j][k].y = e[i][j][k+1].y;
          e[i][j][k].z = e[i][j][k+1].z;
          d = setup->dz[0][i][j][k+1];
        } else if (setup->point_type[i][j][k-1] == CONTACT_EDGE &&
                   d < setup->dz[1][i][j][k-1]) {
          e[i][j][k].x = e[i][j][k-1].x;
          e[i][j][k].y = e[i][j][k-1].y;
          e[i][j][k].z = e[i][j][k-1].z;
          d = setup->dz[1][i][j][k-1];
        }
      }
    }
  }
  */
  if (strstr(fname, "unf")) {
    fprintf(fp, "#\n## start of unformatted data\n");
    for (i = 0; i < setup->numx; i++) {
      for (j = 0; j < setup->numy; j++) {
        if (fwrite(e[i][j], sizeof(vector), setup->numz, fp) != setup->numz) {
          error("Error while writing %s\n", fname);
        }
      }
    }
  } else {
  fprintf(fp, "#\n## x[mm]  y[mm]  z[mm]  Ex[V/cm]  Ey[V/cm]  Ez[V/cm]  |E|[V/cm] V[V]\n");
#define SQ(x) ((x)*(x))
    for (i = 0; i < setup->numx; i++) {
      for (j = 0; j < setup->numy; j++) {
        for (k = 0; k < setup->numz; k++) {
          fprintf(fp, "%6.2f %6.2f %6.2f\t %10e %10e %10e %10e  %10e\n",
                  setup->x0 + i*grid, setup->y0 + j*grid, setup->zmin + k*grid, 		
                  e[i][j][k].x, e[i][j][k].y, e[i][j][k].z,
                  sqrt(SQ(e[i][j][k].x) + SQ(e[i][j][k].y) + SQ(e[i][j][k].z)),
                  v[i][j][k]);
        }
      }
    }
  }

  fprintf(fp2, "#\n## x[mm]  y[mm]  z[mm]  Ex[V/cm]  Ey[V/cm]  Ez[V/cm]  |E|[V/cm] V[V]\n");
  i = j = setup->numx/2;
  for (j = 0; j < setup->numy; j++) {
    // i = j;
    for (k = 0; k < setup->numz; k++) {
      fprintf(fp2, "%6.2f %6.2f %6.2f\t %10e %10e %10e %10e  %10e\n",
              setup->x0 + i*grid, setup->y0 + j*grid, setup->zmin + k*grid, 		
              e[i][j][k].x, e[i][j][k].y, e[i][j][k].z,
              sqrt(SQ(e[i][j][k].x) + SQ(e[i][j][k].y) + SQ(e[i][j][k].z)),
              v[i][j][k]);
    }
  }
  for (i = 0; i < setup->numx; i++) {
    for (j = 0; j < setup->numy; j++) { 
      free(e[i][j]);
    }
    free(e[i]);
  }
  free(e);
  
  fclose(fp);
  fclose(fp2);

  if (1) { /* write point_type to output file */
    i = j = setup->numx/2;
    fp = fopen("fields/point_type.dat", "w");
    for (j = 0; j < setup->numy; j++) {
      // i = j;
      for (k = 0; k < setup->numz; k++) {
        fprintf(fp, "%6.2f %6.2f %6.2f 0 0 0 %2d\n",
                setup->x0+i*grid, setup->y0+j*grid, setup->zmin+k*grid, 		
                setup->point_type[i][j][k]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }

 return 0;
 }

 int write_wp(GRETA_Siggen_Setup *setup, char *fname) {
  int i,j,k;
  FILE *fp = NULL;
  float grid = setup->xtal_grid;

  if ((fp = fopen(fname, "w")) == NULL) {
    error("Failed to open output file: %s\n", fname);
    return -1;
  }
  /* copy configuration parameters to output file */
  report_config(fp, setup->config_file_name);

  if (strstr(fname, "unf")) {
    fprintf(fp, "#\n## start of unformatted data\n");
    for (i = 0; i < setup->numx; i++) {
      for (j = 0; j < setup->numy; j++) {
        for (k = 0; k < setup->numz; k++) {
          setup->dz[0][0][0][k] = setup->v[0][i][j][k]; // re-using float array to convert double to float
        }
        if (fwrite(setup->dz[0][0][0], sizeof(float), setup->numz, fp) != setup->numz) {
          error("Error while writing %s\n", fname);
        }
      }
    }
  } else {
    fprintf(fp, "#\n## x[mm]  y[mm]  z[mm]   WP\n");
    for (i = 0; i < setup->numx; i++) {
      for (j = 0; j < setup->numy; j++) {
        for (k = 0; k < setup->numz; k++) {
          fprintf(fp, "%6.2f %6.2f %6.2f\t %10e\n",
                  setup->x0 + i*grid, setup->y0 + j*grid, setup->zmin + k*grid, 		
                  setup->v[0][i][j][k]);
        }
      }
    }
  }
  
  fclose(fp);
  return 0;
}


static int malloc_arrays(GRETA_Siggen_Setup *setup) {
  int i, j, k;

  if ((setup->v[0]       = malloc(setup->numx * sizeof(setup->v[0][0]))) == NULL ||
      (setup->v[1]       = malloc(setup->numx * sizeof(setup->v[1][0]))) == NULL) {
    error("malloc failed\n");
    return -1;
  }
  for (i = 0; i < setup->numx; i++) {
    if ((setup->v[0][i]       = malloc(setup->numy * sizeof(*setup->v[0][i]))) == NULL ||
	(setup->v[1][i]       = malloc(setup->numy * sizeof(*setup->v[1][i]))) == NULL) {
      error("malloc failed\n");
      return -1;
    }
    for (j = 0; j < setup->numy; j++) {
      if ((setup->v[0][i][j]       = malloc((setup->numz+1) * sizeof(*setup->v[0][i][j]))) == NULL ||
          (setup->v[1][i][j]       = malloc((setup->numz+1) * sizeof(*setup->v[1][i][j]))) == NULL) {
	error("malloc failed\n");
	return -1;
      }
    }
  }

  if ((setup->voxel_impurity_z = malloc(setup->numz * sizeof(double)))    == NULL ||
      (setup->v[0]       = malloc(setup->numx * sizeof(setup->v[0][0])))  == NULL ||
      (setup->v[1]       = malloc(setup->numx * sizeof(setup->v[1][0])))  == NULL ||
      (setup->e          = malloc(setup->numx * sizeof(setup->e[0])))     == NULL ||
      (setup->dx[0]      = malloc(setup->numx * sizeof(setup->dx[0][0]))) == NULL ||
      (setup->dx[1]      = malloc(setup->numx * sizeof(setup->dx[1][0]))) == NULL ||
      (setup->dy[0]      = malloc(setup->numx * sizeof(setup->dy[0][0]))) == NULL ||
      (setup->dy[1]      = malloc(setup->numx * sizeof(setup->dy[1][0]))) == NULL ||
      (setup->dz[0]      = malloc(setup->numx * sizeof(setup->dz[0][0]))) == NULL ||
      (setup->dz[1]      = malloc(setup->numx * sizeof(setup->dz[1][0]))) == NULL ||
      (setup->point_type = malloc(setup->numx * sizeof(setup->point_type[0]))) == NULL) {
    error("malloc failed\n");
    return -1;
  }
  for (i = 0; i < setup->numx; i++) {
    if ((setup->v[0][i]       = malloc(setup->numy * sizeof(*setup->v[0][i])))  == NULL ||
	(setup->v[1][i]       = malloc(setup->numy * sizeof(*setup->v[1][i])))  == NULL ||
        (setup->e[i]          = malloc(setup->numy * sizeof(*setup->e[i])))     == NULL ||
        (setup->dx[0][i]      = malloc(setup->numy * sizeof(*setup->dx[0][i]))) == NULL ||
        (setup->dx[1][i]      = malloc(setup->numy * sizeof(*setup->dx[1][i]))) == NULL ||
        (setup->dy[0][i]      = malloc(setup->numy * sizeof(*setup->dy[0][i]))) == NULL ||
        (setup->dy[1][i]      = malloc(setup->numy * sizeof(*setup->dy[1][i]))) == NULL ||
        (setup->dz[0][i]      = malloc(setup->numy * sizeof(*setup->dz[0][i]))) == NULL ||
        (setup->dz[1][i]      = malloc(setup->numy * sizeof(*setup->dz[1][i]))) == NULL ||
        (setup->point_type[i] = malloc(setup->numy * sizeof(*setup->point_type[i]))) == NULL) {
      error("malloc failed\n");
      return -1;
    }
    for (j = 0; j < setup->numy; j++) {
      if ((setup->v[0][i][j]       = malloc((setup->numz+1) * sizeof(*setup->v[0][i][j]))) == NULL ||
          (setup->v[1][i][j]       = malloc((setup->numz+1) * sizeof(*setup->v[1][i][j]))) == NULL ||
          (setup->e[i][j]          = malloc(setup->numz * sizeof(*setup->e[i][j])))     == NULL ||
          (setup->dx[0][i][j]      = malloc(setup->numz * sizeof(*setup->dx[0][i][j]))) == NULL ||
          (setup->dx[1][i][j]      = malloc(setup->numz * sizeof(*setup->dx[1][i][j]))) == NULL ||
          (setup->dy[0][i][j]      = malloc(setup->numz * sizeof(*setup->dy[0][i][j]))) == NULL ||
          (setup->dy[1][i][j]      = malloc(setup->numz * sizeof(*setup->dy[1][i][j]))) == NULL ||
          (setup->dz[0][i][j]      = malloc(setup->numz * sizeof(*setup->dz[0][i][j]))) == NULL ||
          (setup->dz[1][i][j]      = malloc(setup->numz * sizeof(*setup->dz[1][i][j]))) == NULL ||
          (setup->point_type[i][j] = malloc(setup->numz * sizeof(*setup->point_type[i][j]))) == NULL) {
	error("malloc failed\n");
	return -1;
      }
      for (k = 0; k < setup->numz; k++) {
        setup->dx[0][i][j][k] = setup->dy[0][i][j][k] = setup->dz[0][i][j][k] = 1.0;
        setup->dx[1][i][j][k] = setup->dy[1][i][j][k] = setup->dz[1][i][j][k] = 1.0;
        setup->e[i][j][k] = 0;
        setup->point_type[i][j][k] = 99;
      }
    }
  }

  return 0;
}

static int fname_insert(char fname[256], char dest[256], char *str) {
  char *cp1, *cp2;

  strncpy(dest, fname, 256);
  cp1 = rindex(dest,'.');
  cp2 = rindex(fname,'.');
  if (cp1 == NULL) cp1 = dest + strlen(dest);
  strncpy(cp1, str, 256 - (cp1 - dest));
  cp1 = dest + strlen(dest);
  if (cp2 != NULL)
    strncpy(cp1, cp2, 256 - (cp1 - dest));

  return 0;
}
