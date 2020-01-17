
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "greta_siggen.h"
#include "calc_signal.h"

#define STEP 0.2
#define TRY_BYTES 0


int merge_quads(void);

int main(int argc, char **argv)
{

  GRETA_Siggen_Setup setup;
  FILE   *fp;
  struct point cart;
  int    seg, nsegs, time_steps, i, j, k, n, len, use_byte = 0, snum = 0, quadrant=0;
  float  **s;
  short  ss[32000] = {0}, a, b;
  signed char *sc;
  char   fname[256];
  float  xlo, ylo, dxy;
  
  sc = (signed char *) &ss[0];
  /* read config file and initialize */
  if (argc < 2) {
    printf("Usage: %s <config_file_name> [quadrant]\n"
           "       quadrant = 0 or absent to calculate all quadrants\n"
           "       quadrant = 1-4 for (x<0,y<0), (x<0,y>0), (x>0,y<0), (x>0,y>0) respectively\n"
           "       Then run %s merge to merge the four quadrant files.\n",
           argv[0],argv[0]);
    return 1;
  }

 
  if (!strcmp(argv[1], "merge")) return merge_quads(); /* do the merge of the four quadrants */

  if (argc == 3) quadrant = atoi(argv[2]);
  if (quadrant < 0 || quadrant > 4) quadrant = 0;
  if (quadrant) {
    printf("Calculating only quadrant %d\n", quadrant);
    sprintf(fname, "sgrid%d.mat", quadrant);
  } else {
    sprintf(fname, "sgrid.mat");
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
  fp = fopen(fname, "w");
  if (fp == NULL) {
    fprintf(stderr,"Error! Unable to open %s\n", fname);
    return 0;
  }
  ss[0] = nsegs;
  ss[1] = time_steps;
  ss[2] = STEP*100.0 + 0.5;  // grid step size in 0.01 mm
  ss[3] = setup.xtal_radius*100.0 + 0.5;  // xtal radius in 0.01 mm
  ss[4] = setup.xtal_length*100.0 + 0.5;  // xtal length in 0.01 mm
  fwrite(ss, sizeof(short), 5, fp);
  
  xlo = ylo = -setup.xtal_radius + STEP/2.0;
  if (quadrant > 2) xlo = STEP/2.0;
  if (quadrant == 2 || quadrant == 4) ylo = STEP/2.0;
  dxy = 2.0 * setup.xtal_radius - STEP/2.0;
  if (quadrant) dxy = setup.xtal_radius - STEP/2.0;

  for (cart.x = xlo; cart.x < xlo+dxy; cart.x += STEP) {
    printf("x = %7.3f\n", cart.x); 
    ss[0] = (cart.x + setup.xtal_radius) / STEP + 0.1;
    for (cart.y = ylo; cart.y < ylo+dxy; cart.y += STEP) {
      if (cart.x*cart.x + cart.y*cart.y > setup.xtal_radius*setup.xtal_radius) continue;
      ss[1] = (cart.y + setup.xtal_radius) / STEP + 0.1;
      for (cart.z = STEP/2.0; cart.z < setup.xtal_length; cart.z += STEP) {
        ss[2] = cart.z / STEP + 0.1;
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
          printf("snum %5d: %4d shorts of data saved; x, y, z = %6.2f %6.2f %6.2f  ->  seg %2d\n",
                 snum, len, cart.x, cart.y, cart.z, seg);
        snum++;
      }
    }
  }

  printf("\n All done. %d signals saved.\n\n", snum);
  return 0;
}


int merge_quads(void)
{

  FILE   *fp, *fp2;
  int    seg, len, quad, j, ix, iy, iz;
  short  ss[32000];
  char   fname[256], line[256];
  long long int   pos_end[4];

  int    rm_used_files = 0;

  // first check that all four files (sgrid[1-4].mat) exist, and are consistent

  fp = fopen("sgrid1.mat", "r+");
  if (fp == NULL) {
    printf("Error! Unable to open %s\n", fname);
    return 0;
  }
  fread(ss, sizeof(short), 5, fp);
  fseek(fp, 0, SEEK_END);
  pos_end[0] = ftell(fp);
  printf("Quadrant 1 ends at byte %lld\n", pos_end[0]);
  
  for (quad = 2; quad <= 4; quad++) {
    sprintf(fname, "sgrid%d.mat", quad);
    fp2 = fopen(fname, "r");
    if (fp2 == NULL) {
      printf("Error! Unable to open %s\n", fname);
      return 0;
    }
    fread(ss+5, sizeof(short), 5, fp2);
    for (j=0; j<5; j++) {
      if (ss[j] != ss[j+5]) {
        printf("Error! time steps/step size/crystal size does not match! %s\n", fname);
        return 0;
      }
    }
    fclose(fp2);
  }

  /* files all exist and data match */
  printf("files all exist and data match.\n");
  for (quad = 2; quad <= 4; quad++) {
    sprintf(fname, "sgrid%d.mat", quad);

    printf("\n ... copying data from %s\n", fname);
    fp2 = fopen(fname, "r");
    fread(ss+5, sizeof(short), 5, fp2);
    while ((j=fread(ss, 2, sizeof(ss)/2, fp2)))
      fwrite(ss, 2, j, fp);
    fclose(fp2);

    printf(" ... validating copied data\n");
    fseek(fp, pos_end[quad-2], SEEK_SET);
    while (fread(ss, sizeof(short), 5, fp)) {  // small header for signal
      ix   = ss[0];
      iy   = ss[1];
      iz   = ss[2];
      seg  = ss[3];
      len  = ss[4];
      if (ix < 0 || iy < 0 || iz < 0 ||
          seg < -1 || seg > 37 ||
          len < 42 || len > 4000) {
        printf("Error!: ix,iy,iz, seg, len: = %2d %2d %2d, %2d, %4d\n Aborting!\n",
               ix,iy,iz, seg, len);
        return -1;
      }
      fread(ss, sizeof(short), len-5, fp);  // remainder of signal
    }

    fseek(fp, 0, SEEK_END);
    pos_end[quad-1] = ftell(fp);
    printf("Quadrant %d ends at byte %lld\n", quad, pos_end[quad-1]);

    sprintf(line, "rm %s", fname);
    if (rm_used_files) system(line);
  }

  sprintf(line, "mv sgrid1.mat sgrid.mat");
  if (rm_used_files) system(line);

  return 0;
}
