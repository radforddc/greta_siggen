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
//#define OVER_RELAX_F2 0.90
//#define OVER_RELAX_F2 (0.95 * (1.0 - 1.0/(double)(1+iter/4)))
//#define OVER_RELAX_F2   (0.8  * (1.0 - 1.0/(double)(1+iter/4)))
//#define OVER_RELAX_F2   (0.96 * (5000.0/(5000.0+iter)))
static float overra[] = {0.96,0.97,0.94,0.40,0.35,0.30,0.25,0.25};  // grid = 0.5
static float overrb[] = {0.985,0.98,0.975,0.96,0.95,0.92,0.30,0.80};  // grid = 0.25
#define OVER_RELAX_FA   (iter<1800 ? overra[(iter)/251] : 0.2)
#define OVER_RELAX_FB   (iter<1950 ? overrb[(iter)/251] : 0.2)

#define OVER_RELAX_F2   (setup->xtal_grid > 0.4 ? OVER_RELAX_FA : OVER_RELAX_FB)



static int grid_init(GRETA_Siggen_Setup *setup);
static int ev_calc(GRETA_Siggen_Setup *setup);
static int wp_calc(GRETA_Siggen_Setup *setup);
static int malloc_arrays(GRETA_Siggen_Setup *setup);
static int free_arrays(GRETA_Siggen_Setup *setup);
static int write_ev(GRETA_Siggen_Setup *setup, char*fname);
static int write_wp(GRETA_Siggen_Setup *setup, char*fname);
static int fname_insert(char fname[256], char dest[256], char *str);
int do_relaxation(GRETA_Siggen_Setup *setup, int efld_calc);

int main(int argc, char *argv[]) {
  GRETA_Siggen_Setup setup;

  char opt;
  char optstring[] = "n:r:g:f:v:w:ht:";
  extern char *optarg;
  extern int optind, opterr, optopt;

  if (argc%2 != 0 || !strncmp(argv[1], "-h", 2)) {
    printf("usage: %s <config_file_name>\n"
           "        -b bias_volts\n"
           "        -w {0,1}    (do_not/do write the field file)\n"
           "        -p {0,1}    (do_not/do write the WP file)\n"
           "        -z rho_spectrum_file_name\n", argv[0]);
    return 0;
  }
  if (read_config(argv[1], &setup)) return 1;
  strncpy(setup.config_file_name, argv[1], 256);

  while ((opt = getopt(argc, argv, optstring)) >= 0) {
    if (opt == 'h' || opt == '?') {
      printf("usage: %s <config_file_name>\n"
             "        -b bias_volts\n"
             "        -w {0,1}    (do_not/do write the field file)\n"
             "        -p {0,1}    (do_not/do write the WP file)\n"
             "        -z rho_spectrum_file_name\n", argv[0]);
      return 0;
    }
    if (opt == 'b') {
      setup.xtal_HV = atoi(optarg);
    } else if (opt == 'w') {
      setup.write_field = atoi(optarg);
    } else if (opt == 'p') {
      setup.write_WP= atoi(optarg);
    } else {
      printf("unknown option: %c\n", opt);
    } 
    // FIXME - implement other option: -z
  }

  /* modify core dimensions by Li thickness */
  setup.core_radius += setup.Li_thickness;
  setup.core_length += setup.Li_thickness;
  setup.core_gap    -= setup.Li_thickness;
  setup.core_bullet_radius  += setup.Li_thickness;
  setup.core_bullet_radius  += setup.Li_thickness;
  setup.bottom_taper_width  += setup.Li_thickness * 0.7;
  setup.bottom_taper_length += setup.Li_thickness * 0.7;
  if (setup.Li_thickness > 0.1)
    printf("*** Note: *** Increased core contact dimensions by Li_thickness = %.2f mm\n",
           setup.Li_thickness);

  printf("number of iterations: %d\n", setup.max_iterations);

  if (grid_init(&setup) != 0) {
    error("failed to init field calculations\n");
    return 1;
  }
  if (setup.write_field) ev_calc(&setup);
  if (setup.write_WP)    wp_calc(&setup);
     
  return 0;
}


static int grid_init(GRETA_Siggen_Setup *setup) {

  setup->xmin = setup->ymin = - setup->xtal_radius;
  setup->xmax = setup->ymax = setup->xtal_radius;
  setup->zmin = 0; setup->zmax = setup->xtal_length;

  setup->numx = (int)rint((setup->xmax - setup->xmin)/setup->xtal_grid) + 1;
  setup->numy = (int)rint((setup->ymax - setup->ymin)/setup->xtal_grid) + 1;
  setup->numz = (int)rint((setup->zmax - setup->zmin)/setup->xtal_grid) + 1;
  printf("numx, numy, numz: %d %d %d\n", setup->numx, setup->numy, setup->numz);
  
  if (malloc_arrays(setup) < 0) return -1;

  if ((setup->ncontacts = geometry_init(setup)) <0) {
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
  int   i, j, k;
  int   old, new, iter, undep;
  float sum_dif, max_dif, dif, save_dif;
  float mean, z, x, y, r;


  old = 1; new = 0;
  for (i = 0; i < setup->numx; i++) {
    for (j = 0; j < setup->numy; j++) {
      for (k = 0; k < setup->numz; k++) {
	if (setup->point_type[i][j][k] == OUTSIDE) continue;
	if (setup->point_type[i][j][k] != INSIDE &&
            setup->point_type[i][j][k] != CONTACT_0 &&
            setup->point_type[i][j][k] != CONTACT_VB &&
	    (i > 0 && i + 1 < setup->numx) &&
            (j > 0 && j + 1 < setup->numy) &&
            (k > 0 && k < setup->numz)) {  // k = numz is passivated surface
          printf("error for point %d %d %d\n", i, j, k);
	}

	if ((i > 0             && setup->point_type[i-1][j][k] == OUTSIDE &&
	    i+1 < setup->numx  && setup->point_type[i+1][j][k] == OUTSIDE) ||
	    (j > 0             && setup->point_type[i][j-1][k] == OUTSIDE &&
             j+1 < setup->numy && setup->point_type[i][j+1][k] == OUTSIDE)) {
	  printf("error2 for point %d %d %d\n", i, j, k);
	}
      }
    }
  }

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
    old = old == 0;
    new = new == 0;
    sum_dif = 0.0;
    max_dif = 0.0;

    /* do relaxation iteration */
    for (i = 1; i < setup->numx-1; i++) {
      x = setup->xmin + i*setup->xtal_grid;
      for (j = 1; j < setup->numy-1; j++) {
	y = setup->ymin + j*setup->xtal_grid;
	r = sqrt(x*x + y*y);
	for (k = 1; k < setup->numz; k++) {
	  z = setup->zmin + k*setup->xtal_grid;
	  if (setup->point_type[i][j][k] != INSIDE) continue;
          // save step difference from previous iteration
          save_dif = setup->v[old][i][j][k] - setup->v[new][i][j][k];
          // if (iter < 2) save_dif = 0; 
          mean = setup->v[old][i-1][j][k] +
                 setup->v[old][i+1][j][k] +
                 setup->v[old][i][j-1][k] +
                 setup->v[old][i][j+1][k] +
                 setup->v[old][i][j][k-1];
          if (k < setup->numz-1) {  // bulk
	    mean += setup->v[old][i][j][k+1];
	  } else {                // passivated surface
	    mean += setup->v[old][i][j][k-1];
	  }

          if (efld_calc) {
            setup->v[new][i][j][k] = mean/6.0 + setup->voxel_impurity_z[k];
          } else {
            setup->v[new][i][j][k] = mean/6.0;
          }
	  dif = setup->v[old][i][j][k] - setup->v[new][i][j][k];
	  if (dif < 0.0) dif = -dif;
	  sum_dif += dif;
	  if (max_dif < dif) max_dif = dif;

          setup->v[new][i][j][k] += OVER_RELAX_F2*save_dif;   // do over-relaxation
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
      if (max_dif/setup->xtal_HV < 1.0e-3/5000.0 && max_dif/sum_dif > 200) break;
    } else {
      if (max_dif < 2.0e-7 && max_dif/sum_dif > 15) break;
    }
  }
  printf("%5d %.3e %.3e\n", 
	 iter, max_dif, sum_dif/(setup->numx*setup->numy*setup->numz));

  undep = 0;
  for (i = 0; i < setup->numx && !undep; i++) {
    for (j = 0; j < setup->numy && !undep; j++) {
      for (k = 0; k < setup->numz && !undep; k++) {
	if (setup->v[new][i][j][k] > setup->xtal_HV) {
	  printf("detector is undepleted!\n");
	  undep = 1;
	}
      }
    }
  }
  printf("...  relaxation done\n");

  return 0;
}


static int ev_calc(GRETA_Siggen_Setup *setup) {
  
  printf("\n\n ---- starting EV calculation --- \n");
  fflush(stdout);
  init_ev_calc(setup);

  printf("Starting potential/field calculation...\n");
  do_relaxation(setup, 1);
  printf("saving E, V in file %s\n", setup->field_name);
  write_ev(setup, setup->field_name);

  return 0;
}

static int wp_calc(GRETA_Siggen_Setup *setup) {
  int  cnum;
  char fname[256], str[256];

  printf("\n\n ---- starting WP calculation --- \n");
  for (cnum = 0; cnum < setup->ncontacts; cnum++) {
    //for (cnum = 33; cnum < setup->ncontacts; cnum++) {
    init_wp_calc(setup, cnum);
    snprintf(str, 256, "-seg%2.2d_%d", cnum, 0);
    fname_insert(setup->wp_name, fname, str);
    
    printf("Starting WP calculation for contact number %d...\n", cnum);
    do_relaxation(setup, 0);
    snprintf(str, 256, "-seg%2.2d", cnum);
    fname_insert(setup->wp_name, fname, str);
    printf("saving WP in file %s\n", fname);
    write_wp(setup, fname);
  }

  free_arrays(setup);

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
  int i,j,k;
  FILE *fp = NULL;
  struct efld{float x; float y; float z;} ***e;


  if ((fp = fopen(fname, "w")) == NULL) {
    error("Failed to open output file: %s\n", fname);
    return -1;
  }
  /* copy configuration parameters to output file */
  report_config(fp, setup->config_file_name);
  fprintf(fp, "#\n# HV bias in fieldgen: %.1f V\n", setup->xtal_HV);
  fprintf(fp, "#\n## x[mm]  y[mm]  z[mm]  Ex[V/cm]  Ey[V/cm]  Ez[V/cm]  |E|[V/cm] V[V]\n");

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
  for (i = 0; i < setup->numx; i++) {
    for (j = 0; j < setup->numy; j++) {
      for (k = 0; k < setup->numz; k++) {
	if (setup->point_type[i][j][k] == OUTSIDE) continue;
	if (setup->point_type[i][j][k] == CONTACT_0) continue;
	if (setup->point_type[i][j][k] == CONTACT_VB) continue;
	if (i == 0 || (setup->point_type[i-1][j][k] == OUTSIDE))
	  e[i][j][k].x = -(setup->v[0][i+1][j][k] - setup->v[0][i][j][k])/setup->xtal_grid;
	else if (i == setup->numx -1 || (setup->point_type[i+1][j][k] == OUTSIDE))
	  e[i][j][k].x = -(setup->v[0][i][j][k] - setup->v[0][i-1][j][k])/setup->xtal_grid;
	else
	  e[i][j][k].x = -(setup->v[0][i+1][j][k] - setup->v[0][i-1][j][k])/setup->xtal_grid/2;
	if (j == 0 || (setup->point_type[i][j-1][k] == OUTSIDE))
	  e[i][j][k].y = -(setup->v[0][i][j+1][k] - setup->v[0][i][j][k])/setup->xtal_grid;
	else if (j == setup->numy -1 || setup->point_type[i][j+1][k] == OUTSIDE)
	  e[i][j][k].y = -(setup->v[0][i][j][k] - setup->v[0][i][j-1][k])/setup->xtal_grid;
	else
	e[i][j][k].y = -(setup->v[0][i][j+1][k] - setup->v[0][i][j-1][k])/setup->xtal_grid/2;
	if (k == 0 || (setup->point_type[i][j][k-1] == OUTSIDE))
	  e[i][j][k].z = -(setup->v[0][i][j][k+1] - setup->v[0][i][j][k])/setup->xtal_grid;
	else if (k == setup->numz -1 || (setup->point_type[i][j][k+1] == OUTSIDE))
	  e[i][j][k].z = -(setup->v[0][i][j][k] - setup->v[0][i][j][k-1])/setup->xtal_grid;
	else
	  e[i][j][k].z = -(setup->v[0][i][j][k+1] - setup->v[0][i][j][k-1])/setup->xtal_grid/2;
      }
    }
  }
#define SQ(x) ((x)*(x))
  for (i = 0; i < setup->numx; i++) {
    for (j = 0; j < setup->numy; j++) {
      for (k = 0; k < setup->numz; k++) {
	fprintf(fp, "%6.2f %6.2f %6.2f\t %10e %10e %10e %10e  %10e\n",
		(setup->xmin+i*setup->xtal_grid), 
		(setup->ymin+j*setup->xtal_grid),
                (setup->zmin+k*setup->xtal_grid), 		
		e[i][j][k].x*10, e[i][j][k].y*10, e[i][j][k].z*10,
		sqrt(SQ(e[i][j][k].x) + SQ(e[i][j][k].y) 
		     + SQ(e[i][j][k].z))*10,
		setup->v[0][i][j][k]);
      }
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
  return 0;
 }

 int write_wp(GRETA_Siggen_Setup *setup, char *fname) {
  int i,j,k;
  FILE *fp = NULL;

  if ((fp = fopen(fname, "w")) == NULL) {
    error("Failed to open output file: %s\n", fname);
    return -1;
  }
  /* copy configuration parameters to output file */
  report_config(fp, setup->config_file_name);
  fprintf(fp, "#\n## x[mm]  y[mm]  z[mm]   WP\n");

  for (i = 0; i < setup->numx; i++) {
    for (j = 0; j < setup->numy; j++) {
      for (k = 0; k < setup->numz; k++) {
	fprintf(fp, "%6.2f %6.2f %6.2f\t %10e\n",
		(setup->xmin + i * setup->xtal_grid),
                (setup->ymin + j * setup->xtal_grid),
                (setup->zmin + k * setup->xtal_grid), 
		setup->v[0][i][j][k]);
      }
    }
  }
  
  fclose(fp);
  return 0;
}


static int malloc_arrays(GRETA_Siggen_Setup *setup) {
  int i,j;

  if ((setup->voxel_impurity_z = malloc(setup->numz * sizeof(double))) == NULL ||
      (setup->v[0] = malloc(setup->numx * sizeof(setup->v[0][0]))) == NULL ||
      (setup->v[1] = malloc(setup->numx * sizeof(setup->v[0][0]))) == NULL) {
    error("malloc failed\n");
    return -1;
  }
  if ((setup->point_type = malloc(setup->numx*sizeof(setup->point_type))) == NULL) {
    error("malloc failed\n");
    return -1;
  }
  for (i = 0; i < setup->numx; i++) {
    if ((setup->v[0][i] = malloc(setup->numy*sizeof(*setup->v[0][i]))) == NULL ||
	(setup->v[1][i] = malloc(setup->numy*sizeof(*setup->v[1][i]))) == NULL ||
	(setup->point_type[i] = malloc(setup->numy*sizeof(setup->point_type[i]))) == NULL) {
      error("malloc failed\n");
      return -1;
    }
    for (j = 0; j < setup->numy; j++) {
      if ((setup->v[0][i][j] = malloc(setup->numz*sizeof(*setup->v[0][i][j]))) == NULL ||
	  (setup->v[1][i][j] = malloc(setup->numz*sizeof(*setup->v[1][i][j]))) == NULL ||
	  (setup->point_type[i][j] = calloc(setup->numz,sizeof(setup->point_type[i][j]))) == NULL) {
	error("malloc failed\n");
	return -1;
      }     
    }
  }

  return 0;
}

static int free_arrays(GRETA_Siggen_Setup *setup) {
  int i, j;
  for (i = 0; i < setup->numx; i++) {
    for (j = 0; j < setup->numy; j++) {
      free(setup->v[0][i][j]);
      free(setup->v[1][i][j]);
      free(setup->point_type[i][j]);
    }
    free(setup->v[0][i]);
    free(setup->v[1][i]);
    free(setup->point_type[i]);
  }
  free(setup->v[0]);
  free(setup->v[1]);
  free(setup->point_type);

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
