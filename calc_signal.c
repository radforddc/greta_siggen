 /* calc_signal.c -- based on m3d2s.f by I-Yang Lee
 * Karin Lagergren
 *
 * This module contains the main interface to the signal calculation
 * code. 
 *
 * To use: 
 * -- call signal_calc_init. This will initialize geometry, fields,
 *       drift velocities etc.
 * -- call hit_segment/get_signal
 *
 *  David Radford   Oct 2019
 *  Updated code by Karin Lagergren to match new config file arrangement and new fieldgen
 *
 */
/* TODO: see FIXME's below
   charge_trapping is just a placeholder ATM. Should it be defined
   in the fields module?
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "greta_siggen.h"
#include "calc_signal.h"
#include "point.h"
#include "detector_geometry.h"
#include "fields.h"

#define HOLE_CHARGE 1.0
#define ELECTRON_CHARGE -1.0
#define ONE_HIT_SEG 1 // try to allow only one outside segment to have net charge

/* prototypes for module-private functions */
static int   segment_max_wp(GRETA_Siggen_Setup *setup, float *wp, float thresh);
static int   make_signal(GRETA_Siggen_Setup *setup, point pt, float **signal, float q);
static int   sum_signal(GRETA_Siggen_Setup *setup, float **s);         //modifies s
static float charge_trapping(point pt, float distance, float q);
static int   zero_signal(GRETA_Siggen_Setup *setup, float **signal);
// static int   rc_integrate(float *s_in, float *s_out, float tau, int time_steps);


/* signal_calc_init
   read setup from configuration file,
   then read the electric field and weighting potentials,
   and initialize the signal calculation variables and nsegments
   returns 0 for success
*/
int signal_calc_init(GRETA_Siggen_Setup *setup, int *nsegs)
{
  char fname[256];

  strncpy(fname, setup->config_file_name, 256);
  tell(TERSE, "Reading configuration file: %s\n", fname);
  if (read_config(fname, setup) ||
      (setup->nsegments = geometry_init(setup)) <= 0) {
    error("Setup of detector geometry failed\n");
    return -1;
  }
  strncpy(setup->config_file_name, fname, 256);
  setup->ntsteps_out = setup->time_steps_calc * setup->step_time_calc / setup->step_time_out;

  /* modify core dimensions by Li thickness */
  setup->core_radius += setup->Li_thickness;
  setup->core_length += setup->Li_thickness;
  setup->core_gap    -= setup->Li_thickness;
  setup->core_bullet_radius  += setup->Li_thickness;
  setup->core_bullet_radius  += setup->Li_thickness;
  setup->bottom_taper_width  += setup->Li_thickness * 0.7;
  setup->bottom_taper_length += setup->Li_thickness * 0.7;
  if (setup->Li_thickness > 0.1)
    tell(NORMAL,
         "*** Note: *** Increased core contact dimensions by Li_thickness = %.2f mm\n",
         setup->Li_thickness);

  tell(TERSE, "Reading fields\n");
  if (field_setup(setup) != 0) {
    error("Field setup failed\n");
    return -1;
  }
  *nsegs = setup->nsegments;
  
  if ((setup->dpath_e = malloc(setup->time_steps_calc * sizeof(*setup->dpath_e))) == NULL ||
      (setup->dpath_h = malloc(setup->time_steps_calc * sizeof(*setup->dpath_h))) == NULL) {
    error("Malloc failed\n");
    return -1;
  }

  tell(TERSE, "Setup of signal calculation done\n");
  return 0;
}

/* hit_segment
   return the segment number for interaction at point pt
   returns -1 if point is outside crystal
*/
int hit_segment(GRETA_Siggen_Setup *setup, point pt)
{
  return get_signal(setup, pt, NULL);
}

/* get_signal
   calculate signal for point pt. Result is placed in signal_out array
   which is assumed to have the appropriate size (nsegments * setup->ntsteps_out)
   returns segment number or -1 if outside crystal
   if signal_out == NULL => no signal is stored
*/
int get_signal(GRETA_Siggen_Setup *setup, point pt, float **signal_out)
{
  static float **signal;
  static int tsteps = 0, nsegs = 0;
  float **tmp = NULL;
  char tmpstr[MAX_LINE];
  int  segment, i, j, comp_f;

  /*first time -- allocate memory for signal array */
  if (nsegs != setup->nsegments || tsteps != setup->time_steps_calc) {
    if (nsegs != 0) tmp = signal;
    if ((signal = malloc(setup->nsegments*sizeof(*signal))) == NULL) {
      error("malloc failed in hit_segment\n");
      return -1;
    }
    for (j = 0; j < setup->nsegments; j++) {
      if (tsteps != 0) free(tmp[j]);
      if ((signal[j] = malloc(setup->time_steps_calc*sizeof(*signal[j]))) == NULL)
	return -1;
    }
    if (tmp) free(tmp);
    nsegs = setup->nsegments;
    tsteps = setup->time_steps_calc;
  }

  pt_to_str(tmpstr, MAX_LINE, pt);
  if (!in_crystal(setup, pt)) return -1;
  tell(CHATTY, "point %s is in crystal\n", tmpstr);

  zero_signal(setup, signal);
  segment = -1;
  memset(setup->dpath_e, 0, setup->time_steps_calc*sizeof(*setup->dpath_e));
  memset(setup->dpath_h, 0, setup->time_steps_calc*sizeof(*setup->dpath_h));

  tell(CHATTY, " @@@@@ Signal for %s\n", tmpstr);
  if (make_signal(setup, pt, signal, ELECTRON_CHARGE) ||
      make_signal(setup, pt, signal, HOLE_CHARGE)) return -1;
  /* make_signal returns 0 for success */
  sum_signal(setup, signal);

  /*figure out which segment has net charge*/
  for (i = 0; i < setup->nsegments-1; i++) {
    if (signal[i][setup->time_steps_calc-1] > NET_SIGNAL_THRESH) {
      if (segment >= 0) {
	error("found more than one segment with net charge\n");
	return -1;
      }
      segment = i;
    }
  }

  if (signal_out != NULL) {
    /* now, compress the signal and place it in the signal_out array */
    comp_f = setup->time_steps_calc/setup->ntsteps_out;
    for (i = 0; i < setup->nsegments; i++) {
      for (j = 0; j < setup->ntsteps_out; j++) signal_out[i][j] = 0;
      /* truncate the signal if setup->time_steps_calc % setup->ntsteps_out != 0 */
      for (j = 0; j < setup->ntsteps_out*comp_f; j++)
	signal_out[i][j/comp_f] += signal[i][j]/comp_f;
    }
  }
  return segment;
}


static int zero_signal(GRETA_Siggen_Setup *setup, float **signal)
{
  int i, j;
  
  for (i = 0; i < setup->nsegments; i++) {
    for (j = 0; j < setup->time_steps_calc; j++) signal[i][j] = 0.0;
  }
  return 0;
}


/* make_signal
   Generates the signal originating at point pt, for charge q
   returns 0 for success
*/
static int make_signal(GRETA_Siggen_Setup *setup, point pt, float **signal, float q)
{
  static float *wpot, *wpot2, *wpot_old, *dwpot;
  static int wp_size;

  char tmpstr[MAX_LINE];
  point new_pt, prev_pt;
  vector v, dx;
  float dist;
  int t, n, largest_wp_seg, i, j, vr, keep_drifting = 1;


  if (wp_size != setup->nsegments) { /* first time called */
    if (wpot != NULL) free(wpot);
    if (wpot2 != NULL) free(wpot2);
    if (wpot_old != NULL) free(wpot_old);
    if (dwpot != NULL) free(dwpot);
    if ((wpot = malloc(setup->nsegments*sizeof(*wpot))) == NULL ||
        (wpot2 = malloc(setup->nsegments*sizeof(*wpot2))) == NULL ||
	(wpot_old = malloc(setup->nsegments*sizeof(*wpot_old))) == NULL ||
	(dwpot = malloc(setup->nsegments*sizeof(*dwpot))) == NULL) {
      error("malloc failed in make_signal\n");
      exit(1);
    }
    wp_size = setup->nsegments;
  }

  prev_pt = new_pt = pt;
  // int cf = setup->step_time_out / setup->step_time_calc;
  // printf("cf = %d\n", cf); fflush(stdout);
  for (t = 0; ((vr = drift_velocity(setup, new_pt, q, &v)) >= 0 && keep_drifting); t++) { 
    if (q > 0) setup->dpath_h[t] = new_pt;
    else setup->dpath_e[t] = new_pt;
    tell(CHATTY, "t: %d  pt: (%.2f %.2f %.2f)\n", t, new_pt.x,new_pt.y, new_pt.z);
    tell(CHATTY, "v: (%e %e %e)\n", v.x, v.y, v.z);
    // if ((t%cf == 0 || t < 2) && wpotentials(setup, new_pt, wpot) != 0) {
    if (wpotentials(setup, new_pt, wpot) != 0) {
      pt_to_str(tmpstr, MAX_LINE, new_pt);
      tell(NORMAL,
	   "Can calculate velocity but not weighting potentials at %s!\n",
	   tmpstr);
      return -1;
    }
    if (0) tell(CHATTY, " >>>>>>>    t: %d WP: %.4f %.4f %.4f %.4f %.4f\n",
                t, wpot[4], wpot[5], wpot[10], wpot[11], wpot[36]);

    /* ------------- DCR added Oct 2019: if core WP is very small or large, then stop drifting */
    // FIXME: hardcoded seg 36 = core
    if ((q * setup->impurity_z0 < 0 &&    // drifting to core    (e in n-type or h in p-type)
         1.0-wpot[36] < 5.0e-5) ||
        (q * setup->impurity_z0 > 0 &&    // drifting to outside (h in n-type or e in p-type)
         wpot[36] < 5.0e-5)) {
      pt_to_str(tmpstr, MAX_LINE, new_pt);
      tell(CHATTY, "Reached full WP at %s; WP[36] = %9.6f\n", tmpstr, wpot[36]);
      keep_drifting = 0;
    }
    /* relax the limit on this check if we have reached a boundary
       i.e field extrapolation is require) */
    if (vr > 0 && new_pt.z < setup->xtal_length - setup->xtal_grid) {  // CHECKME
      if ((q * setup->impurity_z0 < 0 &&      // drifting to core    (e in n-type or h in p-type)
           1.0-wpot[36] < 0.05 &&
           wpot[36] - wpot_old[36] > 0.0) ||  // make sure local drift is in correct direction
          (q * setup->impurity_z0 > 0 &&      // drifting to outside (h in n-type or e in p-type)
           wpot[36] < 0.05 &&
           wpot[36] - wpot_old[36] < 0.0)) {  // make sure local drift is in correct direction
        pt_to_str(tmpstr, MAX_LINE, new_pt);
        tell(CHATTY, "Reached boundary at %s; WP[36] = %9.6f\n", tmpstr, wpot[36]);
        keep_drifting = 0;
      }
    }
    if (t >= setup->time_steps_calc - 2) keep_drifting = 0; // have run out of time...

    for (i = 0; i < setup->nsegments; i++) {
      if (t > 0) signal[i][t] += q*(wpot[i] - wpot_old[i]);
      wpot_old[i] = wpot[i];
    }
    dx = vector_scale(v, setup->step_time_calc);
    prev_pt = new_pt;
    new_pt = vector_add(new_pt, dx);
    dist = vector_length(dx);
    q = charge_trapping(new_pt, dist, q);   // FIXME
  }

  if (t == 0) {
    pt_to_str(tmpstr, MAX_LINE, pt);
    tell(CHATTY, "The starting point %s is outside the field.\n", tmpstr);
    return -1;
  }
  /*check if we have drifted out of back of detector (not to contact )*/
  if (new_pt.z > setup->xtal_length) {
    tell(CHATTY, "Drifted out of back end of detector.\n");
    // return -1;
    new_pt.z = setup->xtal_length; // FIXME: change horizontal velocity??
  }

  largest_wp_seg = segment_max_wp(setup, wpot, WP_THRESH);
  pt_to_str(tmpstr, MAX_LINE, new_pt);
  tell(CHATTY, "Drifted to edge of field grid, point: %s segment: %d q: %.2f  kd: %d\n", 
       tmpstr, largest_wp_seg, q, keep_drifting);

  /* now we are outside the electric field grid
     decide whether we need to keep drifting to make WPs go to zero */
  if (!ONE_HIT_SEG && largest_wp_seg != 36) return 0;     // FIXME: hardcoded seg 36 = core

  /* figure out how much we must drift to get to the crystal boundary */
  for (n = 0; in_crystal(setup, new_pt) &&  n+t < setup->time_steps_calc; n++) {
    new_pt = vector_add(new_pt, dx);
    if (q > 0) setup->dpath_h[t+n] = new_pt;
    else setup->dpath_e[t+n] = new_pt;
    if (//largest_wp_seg == 36 &&          // drifting to core
        n * setup->step_time_calc > 50 &&  // final drift is longer than 50 ns
        new_pt.z < setup->xtal_length)     // not stuck on the passivated surface
      break;                               // no extra steps beyond 50 ns
    
  }
  if (n == 0) n = 1; /* always drift at least one more step */

  tell(CHATTY, "q: %.1f t: %d n: %d ((%.2f %.2f %.2f)=>(%.2f %.2f %.2f))\n", 
       q,t,n, pt.x, pt.y, pt.z, new_pt.x, new_pt.y, new_pt.z);

  if (n + t >= setup->time_steps_calc) {
    tell(CHATTY, "Exceeded maximum number of time steps (%d)\n", setup->time_steps_calc);
    /* check WP's to see if we have produced most of the signal */
    if (wpot[largest_wp_seg] < 0.95 &&
        wpotentials(setup, new_pt, wpot2) != 0) {
      tell(CHATTY, "Cannot finish drifting to make at least 95\% of signal.\n");
      return -1;  /* FIXME DCR: could this be improved? */
    }
    /* drift to new_pt and wpot2 */
    for (i = 0; i < setup->nsegments; i++) {
      dwpot[i] = (wpot2[i] - wpot[i])/n;
    }
  } else {
    /* weighting pot. is 1 at edge for hit segment, 0 for other segments.
       Make it so, gradually if applicable */
    for (i = 0; i < setup->nsegments; i++) {
      dwpot[i] = ((i == largest_wp_seg) - wpot[i])/n;
    }
  }
  if (dwpot[largest_wp_seg] <  0) {
    tell(CHATTY, "Cannot complete drifting;"
         " WP[%2d], dWP[%2d]: %7.4f %7.4f; WP[36], dWP[36]: %7.4f %7.4f\n",
         largest_wp_seg, largest_wp_seg,
         wpot[largest_wp_seg], dwpot[largest_wp_seg], wpot[36], dwpot[36]);
    return -1;  /* FIXME DCR: could this be improved? */
  }

  /* now drift the final n steps */
  tell(CHATTY, " >>> completing drift: t, n = %4d, %4d (%4d);"
       "  WP[%2d], dWP[%2d]: %7.4f %7.4f; WP[36], dWP[36]: %7.4f %7.4f\n",
       t, n, t+n, largest_wp_seg, largest_wp_seg,
       wpot[largest_wp_seg], dwpot[largest_wp_seg], wpot[36], dwpot[36]);
  dx = vector_scale(v, setup->step_time_calc);
  dist = vector_length(dx);
  for (i = 1; i <= n; i++) {
    for (j = 0; j < setup->nsegments; j++) {
      signal[j][i+t-1] += q*dwpot[j];
    }
    q = charge_trapping(prev_pt,dist, q);   // FIXME    
  }

  pt_to_str(tmpstr, MAX_LINE, pt);
  tell(CHATTY, "q:%.2f pt: %s segment: %d\n", q, tmpstr, largest_wp_seg);

  return 0;
}

/* modifies s. each time step will contain the summed signals of 
   all previous time steps */
static int sum_signal(GRETA_Siggen_Setup *setup, float **s)
{
  int i, j;

  for (i = 0; i < setup->nsegments; i++) {
    for (j = 1; j < setup->time_steps_calc; j++) {
      s[i][j] += s[i][j-1];
    }
  }
  return 0;
}

// FIXME -- placeholder function. Even parameter list is dodgy
static float charge_trapping(point pt, float distance, float q)
{
  return q;
}

/* segment_max_wp 
 *  Return the segment number corresponding to the largest w.p.
 * make sure no more than one segment has wp higher than thresh
 */
static int segment_max_wp(GRETA_Siggen_Setup *setup, float *wp, float thresh)
{
  int n, i;
  int segno;
  float wpmax;

  n = 0;
  for (i = 0; i < setup->nsegments; i++) {
    if (wp[i] > thresh) {
      segno = i;
      n++;
      tell(CHATTY, "Segment %d over threshold\n", i);
    }
  }
  if (n == 1) return segno;
  if (n > 1) {
    error(" %d segments over threshold. Your weigthing potential is broken!\n", n);
    return -1;
  }
  n = 0;
  wpmax = thresh/10; //OK? -- FIXME
  for (i = 0; i < setup->nsegments; i++) {
    if (wp[i] > wpmax) {
      segno = i;
      wpmax = wp[i];
      n++;
    }
  }
  if (n) {
    tell(CHATTY, "Largest wp for segment %d\n", segno);
    return segno;
  }
  tell(CHATTY, "Segment_max_wp: no charge collected!\n");
  return -1;
}

/*
int rc_integrate(float *s_in, float *s_out, float tau, int time_steps)
{
  int   j;
  float s_in_old, s;  // DCR: added so that it's okay to
		      //   call this function with s_out == s_in
  
  if (tau < 1.0f) {
    for (j = time_steps-1; j > 0; j--) s_out[j] = s_in[j-1];
    s_out[0] = 0.0;
  } else {
    s_in_old = s_in[0];
    s_out[0] = 0.0;
    for (j = 1; j < time_steps; j++) {
      s = s_out[j-1] + (s_in_old - s_out[j-1])/tau;
      s_in_old = s_in[j];
      s_out[j] = s;
    }
  }
  return 0;
}
*/

/* signal_calc_finalize
 * Clean up (free arrays, close open files...)
 */
int signal_calc_finalize(GRETA_Siggen_Setup *setup)
{
  fields_finalize(setup);
  free(setup->dpath_h);
  free(setup->dpath_e);
  return 0;
}

int drift_path_e(GRETA_Siggen_Setup *setup,point **pp)
{
  *pp = setup->dpath_e;
  return setup->time_steps_calc;
}
int drift_path_h(GRETA_Siggen_Setup *setup,point **pp)
{
  *pp = setup->dpath_h;
  return setup->time_steps_calc;
}
