#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "greta_siggen.h"

int read_config(char *config_file_name, GRETA_Siggen_Setup *setup) {

  /* reads and parses configuration file of name config_file_name
     fills in values of GRETA_Siggen_Setup setup, defined in greta_siggen.h
     returns 0 on success, 1 otherwise
  */

  char key_word[][32] = {
    "xtal_type",
    "xtal_length",
    "xtal_radius",
    "core_length_gap",     // can use gap (i.e. xtal_length - core_length) instead.
                           // This keyword must be before "core_length".
    "core_length",
    "core_radius",
    "core_bullet_radius",
    "bottom_taper_length",
    "bottom_taper_angle",
    "bottom_taper_width",
    "top_bullet_radius",
    "bottom_bullet_radius",
    "Li_thickness",
    "core_offset_x_top",
    "core_offset_y_top",
    "core_offset_x_bottom",
    "core_offset_y_bottom",
    "xtal_axis_rotation",

    "xtal_grid",
    "impurity_z0",
    "impurity_gradient",
    "impurity_quadratic",
    "impurity_surface",
    "impurity_radial_add",
    "impurity_radial_mult",
    "impurity_rpower",
    "xtal_HV",
    "surface_drift_vel_factor",

    "geometry_name",
    "drift_name",
    "field_name",
    "wp_name",

    "xtal_temp",
    "preamp_tau",
    "time_steps_calc",
    "step_time_calc",
    "step_time_out",
    "charge_cloud_size",
    "use_diffusion",
    "energy",
    "verbosity_level",
    "max_iterations",
    "write_field",
    "write_WP",
    ""
  };

  int   ii, i, l, n=0, ok, iint = 0, diameter = 0;
  float fi;
  char  *c, line[256], line2[256], name[256];
  FILE  *file;


  /* initialize everything to zero... */
  memset(setup, 0, sizeof(*setup));
  /* ...except for impurity_radial_mult */
  setup->impurity_radial_mult = 1.0f;      // 1.0 is neutral (no radial gradient)
  setup->surface_drift_vel_factor = 1.0f;  // normal fast drift, same as bulk

  setup->xtal_type = DEFAULT_CRYSTAL_TYPE;


  if (!(file = fopen(config_file_name, "r"))) {
    printf("\nERROR: config file %s does not exist?\n", config_file_name);
    return 1;
  }
  /* read config file */
  while (fgets(line, sizeof(line), file)) {
    n++;
    /* ignore comments and blank lines */
    if (strlen(line) < 3 || *line == ' ' || *line == '\t' || *line == '#') continue;
    /* if line contains "_diam" (for diameter) replace with "_radius"
       this allows the user to specify diametwrs instead of radii */
    diameter = 0;
    if ((c = strstr(line, "_diam"))) {
      diameter = 1;
      if (setup->verbosity >= CHATTY) printf("Line: %s", line);
      strcpy(line2, c+5);
      strcpy(c, "_radius");
      strcpy(c+7, line2);
      if (setup->verbosity >= CHATTY) printf("   -->  %s", line);      
    }

    for (i=0; (l=strlen(key_word[i])) > 0; i++) {
      if (!strncmp(line, key_word[i], l)) {
	/* line starts with key_word[i] */
        if (setup->verbosity >= CHATTY)printf("Line: %s", line);
	if (line[l] != ' ' && line[l] != '\t') {
	  ok = 0;
	} else {
	  // for (c = line + l; *c != ' ' && *c != '\t'; c++) ;
	  /* find next non-white-space char */
	  for (c = line + l; *c == ' ' || *c == '\t'; c++) ;
	  name[0] = 0;
	  ii = iint = 0;
	  fi = 0;
	  if (strstr(key_word[i], "_name")) {
	    /* extract character string for file name */
	    for (ok=0; ok<256 && *c != ' ' && *c != '\t' && *c != '\n' &&  *c != '\r'; ok++) {
	      name[ok] = *c;
	      c++;
	    }
	    name[ok] = '\0';	    
	  } else if (!strncmp("xtal_type", key_word[i], l) ||
                     !strncmp("time_steps_calc", key_word[i], l) ||
		     !strncmp("use_diffusion", key_word[i], l) ||
		     !strncmp("verbosity_level", key_word[i], l) ||
		     !strncmp("max_iterations", key_word[i], l) ||
		     !strncmp("write_field", key_word[i], l) ||
		     !strncmp("write_WP", key_word[i], l)) {
	    /* extract integer value */
	    ok = sscanf(c, "%d", &ii);
	    iint = 1;
	  } else {
	    /* extract float value */
	    ok = sscanf(c, "%f", &fi);
            if (diameter) fi /= 2.0;
	  }
	}
	if (ok < 1) {
	  printf("ERROR reading %s from config file %s\n"
		 "   ...line number %d is: %s",
		 key_word[i], config_file_name, n, line);
	  return 1;
	}

	if (strstr(key_word[i], "verbosity_level")) {
	  setup->verbosity = ii;
	} else if (strstr(key_word[i], "xtal_type")) {
	  setup->xtal_type = ii;
	} else if (strstr(key_word[i], "xtal_length")) {
	  setup->xtal_length = fi;
	} else if (strstr(key_word[i], "xtal_radius")) {
	  setup->xtal_radius = fi;
	} else if (strstr(key_word[i], "top_bullet_radius")) {
	  setup->top_bullet_radius = fi;
	} else if (strstr(key_word[i], "bottom_bullet_radius")) {
	  setup->bottom_bullet_radius = fi;
	} else if (strstr(key_word[i], "core_length")) {
	  setup->core_length = fi;
          /* the user can specify "core_length_gap" = xtal_length - core_length, instead of "core_length" */
          if (strstr(line, "core_length_gap")  && fi < setup->xtal_length) {
            setup->core_length = setup->xtal_length - fi;
            if (setup->verbosity >= CHATTY) printf("   -->  core_length: %f\n", setup->core_length);
          }
	} else if (strstr(key_word[i], "core_radius")) {
	  setup->core_radius = fi;
	} else if (strstr(key_word[i], "core_bullet_radius")) {
	  setup->core_bullet_radius = fi;
	} else if (strstr(key_word[i], "bottom_taper_length")) {
	  setup->bottom_taper_length = fi;
	} else if (strstr(key_word[i], "bottom_taper_angle")) {
	  setup->bottom_taper_angle = fi;
	} else if (strstr(key_word[i], "bottom_taper_width")) {
	  setup->bottom_taper_width = fi;
	} else if (strstr(key_word[i], "Li_thickness")) {
	  setup->Li_thickness = fi;
	} else if (strstr(key_word[i], "core_offset_x_top")) {
	  setup->core_offset_x_top = fi;
	} else if (strstr(key_word[i], "core_offset_y_top")) {
	  setup->core_offset_y_top = fi;
	} else if (strstr(key_word[i], "core_offset_x_bottom")) {
	  setup->core_offset_x_bottom = fi;
	} else if (strstr(key_word[i], "core_offset_y_bottom")) {
	  setup->core_offset_y_bottom = fi;
	} else if (strstr(key_word[i], "xtal_axis_rotation")) {
	  setup->xtal_axis_rotation = fi;

	} else if (strstr(key_word[i], "xtal_grid")) {
	  setup->xtal_grid = fi;
	} else if (strstr(key_word[i], "impurity_z0")) {
	  setup->impurity_z0 = fi;
	} else if (strstr(key_word[i], "impurity_gradient")) {
	  setup->impurity_gradient = fi;
	} else if (strstr(key_word[i], "impurity_quadratic")) {
	  setup->impurity_quadratic = fi;
	} else if (strstr(key_word[i], "impurity_surface")) {
	  setup->impurity_surface = fi;
	} else if (strstr(key_word[i], "impurity_radial_add")) {
	  setup->impurity_radial_add = fi;
	} else if (strstr(key_word[i], "impurity_radial_mult")) {
	  setup->impurity_radial_mult = fi;
	} else if (strstr(key_word[i], "impurity_rpower")) {
	  setup->impurity_rpower = fi;
	} else if (strstr(key_word[i], "xtal_HV")) {
	  setup->xtal_HV = fi;
	} else if (strstr(key_word[i], "surface_drift_vel_factor")) {
	  setup->surface_drift_vel_factor = fi;
	} else if (strstr(key_word[i], "geometry_name")) {
	  strncpy(setup->geometry_name, name, 256);
	} else if (strstr(key_word[i], "drift_name")) {
	  strncpy(setup->drift_name, name, 256);
	} else if (strstr(key_word[i], "field_name")) {
	  strncpy(setup->field_name, name, 256);
	} else if (strstr(key_word[i], "wp_name")) {
	  strncpy(setup->wp_name, name, 256);
	} else if (strstr(key_word[i], "xtal_temp")) {
	  setup->xtal_temp = fi;
	} else if (strstr(key_word[i], "preamp_tau")) {
	  setup->preamp_tau = fi;
	} else if (strstr(key_word[i], "time_steps_calc")) {
	  setup->time_steps_calc = ii;
	} else if (strstr(key_word[i], "step_time_calc")) {
	  setup->step_time_calc = fi;
	} else if (strstr(key_word[i], "step_time_out")) {
	  setup->step_time_out = fi;
	} else if (strstr(key_word[i], "charge_cloud_size")) {
	  setup->charge_cloud_size = fi;
	} else if (strstr(key_word[i], "use_diffusion")) {
	  setup->use_diffusion = ii;
	} else if (strstr(key_word[i], "energy")) {
	  setup->energy = fi;
	} else if (strstr(key_word[i], "max_iterations")) {
	  setup->max_iterations = ii;
	} else if (strstr(key_word[i], "write_field")) {
	  setup->write_field = ii;
	} else if (strstr(key_word[i], "write_WP")) {
	  setup->write_WP = ii;
	} else {
	  printf("ERROR; unrecognized keyword %s\n", key_word[i]);
	  return 1;
	}

	if (setup->verbosity >= CHATTY) {
	  // printf("%s", line);
	  if (iint) {
	    printf("%s: %d\n", key_word[i], ii);
	  } else if (strlen(name) > 0) {
	    printf("%s: %s\n", key_word[i], name);
	  } else {
	    printf("%s: %f\n", key_word[i], fi);
	  }
	}
	break;
      }
    }
  }
  fclose(file);

  if (setup->bottom_taper_angle > 0) {
    /* convert taper angle to taper width */
    if (setup->bottom_taper_length > 0) {
      setup->bottom_taper_width =
        setup->bottom_taper_length * tan(setup->bottom_taper_angle * 3.14159/180.0);
      if (setup->verbosity >= CHATTY)
        printf("  -->  bottom taper width: %f\n", setup->bottom_taper_width);
    }
  } else {
    /* convert taper width to taper angle */
    if (setup->bottom_taper_length > 0 &&
        setup->bottom_taper_width > 0)
      setup->bottom_taper_angle =
        atan(setup->bottom_taper_width/setup->bottom_taper_length) * 180.0/3.14159;
    if (setup->bottom_taper_angle > 0)
      if (setup->verbosity >= CHATTY)
        printf("  -->  bottomtaper angle: %f\n", setup->bottom_taper_angle);
  }
  setup->core_gap = setup->xtal_length - setup->core_length;

  /* some consistency checks */
  if (setup->core_bullet_radius > setup->core_radius)
    setup->core_bullet_radius = setup->core_radius;
  if (setup->bottom_taper_length > setup->core_length ||
      setup->core_length > setup->xtal_length ||
      (setup->core_radius + setup->bottom_taper_width) > setup->xtal_radius) {
    printf("\nERROR: Inconsistent detector dimensions:\n"
           "    crystal length and radius: %5.2f %5.2f\n"
           "       core length and radius: %5.2f %5.2f\n"
           "bottom taper length and width: %5.2f %5.2f\n\n",
           setup->xtal_length, setup->xtal_radius,
           setup->core_length, setup->core_radius,
           setup->bottom_taper_length, setup->bottom_taper_width);
    return 1;
  }

  return 0;
}
