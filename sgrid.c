
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "greta_siggen.h"
#include "calc_signal.h"

#define STEP 1.0
#define TRY_BYTES 0


int main(int argc, char **argv)
{

  GRETA_Siggen_Setup setup;
  FILE   *fp;
  struct point cart;
  int    seg, nsegs, time_steps, i, j, k, n, len, use_byte = 0, snum = 0;
  float  **s;
  short  ss[10000] = {0}, a, b;
  signed char *sc;
  
  sc = (signed char *) &ss[0];
  /* read config file and initialize */
  if (argc < 2) {
    printf("Usage: %s <config_file_name>\n", argv[0]);
    return 1;
  }

  strncpy(setup.config_file_name, argv[1], 256);
  if (signal_calc_init(&setup, &nsegs) != 0) return 1;
  nsegs = setup.nsegments;
  time_steps = setup.ntsteps_out;

  if ((s = malloc(nsegs*sizeof(*s))) == NULL) {
    printf("Malloc failed\n");
    return 1;
  }
  for (i = 0; i < nsegs; i++){
    if ((s[i] = malloc(time_steps*sizeof(*s[i]))) == NULL) {
      printf("malloc failed\n");
      return 1;
    }
  }

  /* open output file */
  fp = fopen("sgrid.mat", "w");
  if (fp == NULL){
    fprintf(stderr,"Error! Unable to open sgrid.mat\n");
    return 0;
  }
  ss[0] = nsegs;
  ss[1] = time_steps;
  ss[2] = STEP*100.0 + 0.5;  // grid step size in 0.01 mm
  ss[3] = setup.xtal_radius*100.0 + 0.5;  // xtal radius in 0.01 mm
  ss[4] = setup.xtal_length*100.0 + 0.5;  // xtal length in 0.01 mm
  fwrite(ss, sizeof(short), 5, fp);
  
  for (cart.x = -setup.xtal_radius; cart.x < setup.xtal_radius; cart.x += STEP) {
    printf("x = %7.3f\n", cart.x); 
    ss[0] = (cart.x + setup.xtal_radius) * STEP + 0.5;
   for (cart.y = -setup.xtal_radius; cart.y < setup.xtal_radius; cart.y += STEP) {
     if (cart.x*cart.x + cart.y*cart.y > setup.xtal_radius*setup.xtal_radius) continue;
     ss[1] = (cart.y + setup.xtal_radius) * STEP + 0.5;
     for (cart.z = 0; cart.z < setup.xtal_length; cart.z += STEP) {
       ss[2] = cart.z * STEP + 0.5;
       ss[3] = seg = get_signal(&setup, cart, s);
       if (seg < 0) {
         if (0) printf("point not in crystal: (x = %.1f, y = %.1f, z = %.1f)\n", 
                       cart.x, cart.y, cart.z);
         continue;
       }
       /* copy signal data to output array (shorts) */
       k = 5;
       for (i = 0; i < nsegs; i++) {
         for (j = 0; j < time_steps; j++) {
           ss[k + j + 1] = lrintf(s[i][j] * 10000.0);
         }
         for (n = time_steps; n > 0; n--) {
           if (ss[k + n] != 0 && ss[k + n] != 10000 && ss[k + n] != -10000) break;
         }

         if (TRY_BYTES) {
           use_byte = 1;
           a = 0;
           for (j = 0; j < n; j++) {
             b = ss[k + j] - a;
             a = ss[k + j];
             if (b < -128 || b > 127) {  // can't store result in a single byte
               use_byte = 0;
               break;
             }
           }
         }

         if (use_byte) {
           a = 0;
           for (j = 0; j < n; j ++) {
             b = ss[k + j] - a;
             a = ss[k + j];
             sc[2*k + j] = b;
           }
           ss[k] = -(n+1);
           k += (n+3)/2;
         } else {
           ss[k] = n+1;
           k += n+1;
         }
       }
       ss[4] = len = k;
       fwrite(ss, sizeof(short), len, fp);
       if (snum%1000 == 0)
         printf("snum %5d: %4d shorts of data saved; x, y, x = %7.3f %7.3f %7.3f  ->  seg %2d\n",
                snum, len, cart.x, cart.y, cart.z, seg);
       snum++;
     }
   }
  }

  printf("\n All done. %d signals saved.\n\n", snum);
  return 0;
}
