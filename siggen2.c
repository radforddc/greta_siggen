
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <readline/readline.h>
#include <readline/history.h>


/*readline wrapper*/
#define FREE(p) do{free(p); p = NULL; }while(0)

static int rl_gets(char *line, int MAX_LEN){
  static char *line_read = (char *) NULL;

  if (line_read != NULL){
    FREE(line_read);
    line_read = (char *) NULL;
  }

  line_read = readline("> ");
  if (line_read != NULL && *line_read != '\0'){
    add_history(line_read);
  }
  strncpy(line,line_read,MAX_LEN);
  line[MAX_LEN -1] = '\0';
  FREE(line_read);

  return 1;
}


int main(int argc, char **argv)
{

  FILE   *fpi, *fpo;
  int    seg, nsegs, time_steps; //, len;
  int    ix, iy, iz, nx, ny, nz;
  float  radius, length, step;  // x, y, z;
  short  ***seg_dir;
  long long int ***pos_dir, pos;


  if (argc < 2) {
    printf("Usage: %s <input_sg_dir_file_name> <sgrid_file_name>\n", argv[0]);
    return 1;
  }

  /* open input directory file */
  fpi = fopen(argv[1], "r");
  if (fpi == NULL){
    fprintf(stderr,"Error! Unable to open input file argv[1]\n");
    return 0;
  }

  fread(&nsegs     , sizeof(nsegs     ), 1, fpi);
  fread(&time_steps, sizeof(time_steps), 1, fpi);
  fread(&step      , sizeof(step      ), 1, fpi);
  fread(&radius    , sizeof(radius    ), 1, fpi);
  fread(&length    , sizeof(length    ), 1, fpi);

  printf(" %d segs, %d time steps\n"
         " xtal radius x length = %.2f x %.2f; step = %.3f\n",
         nsegs, time_steps, radius, length, step);

  /* malloc space for 3D directory arrays */
  nx = ny = 2.0*radius/step + 0.55;
  nz = length/step + 0.55;
  printf("nx, ny, nz = %d, %d, %d\n\n", nx, ny, nz);

  if ((seg_dir = malloc(nx*sizeof(*seg_dir))) == NULL ||
      (pos_dir = malloc(nx*sizeof(*pos_dir))) == NULL) {
    printf("malloc1 failed\n");
    return 1;
  }
  for (ix = 0; ix < nx; ix++) {
    if ((seg_dir[ix] = malloc(ny*sizeof(*seg_dir[ix]))) == NULL ||
        (pos_dir[ix] = malloc(ny*sizeof(*pos_dir[ix]))) == NULL) {
      printf("malloc2 failed; ix = %d\n", ix);
      return 1;
    }
    for (iy = 0; iy < ny; iy++) {
      if ((seg_dir[ix][iy] = malloc(nz*sizeof(*seg_dir[ix][iy]))) == NULL ||
          (pos_dir[ix][iy] = malloc(nz*sizeof(*pos_dir[ix][iy]))) == NULL) {
        printf("malloc3 failed; ix, iy = %d, %d\n", ix, iy);
        return 1;
      }
      for (iz = 0; iz < nz; iz++) {
        seg_dir[ix][iy][iz]  = -1;
        pos_dir[ix][iy][iz]  = -1;
      }
    }
  }

  /* read segment number directory */
  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      fread(seg_dir[ix][iy], sizeof(***seg_dir), nz, fpi);
    }
  }
  /* read file_position directory */
  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      fread(pos_dir[ix][iy], sizeof(***pos_dir), nz, fpi);
    }
  }
  fclose(fpi);
  printf(" Reading done.\n");
  printf(" %d segs, %d time steps\n"
         " xtal radius x length = %.2f x %.2f; step = %.3f\n\n",
         nsegs, time_steps, radius, length, step);

  /* open input signal file */
  fpi = fopen(argv[2], "r");
  if (fpi == NULL){
    fprintf(stderr,"Error! Unable to open input file argv[2]\n");
    return 0;
  }

  char  cmd[256];
  int   imode = 0;
  short ss[4000] = {0};
  float fs[4096];
  float x, y, z, fwhm=0, tau=0, f;
  int   i, j, k, n, len;

  char  *help = "Available commands:\n i to switch imode\n f <value> to give CC FWHM\n t <value> to give tau\n q to quit\n";
  printf("%s", help);
  while (1) {
    printf("\nEnter %s for signal\n", (imode ? "ix iy iz [grid]" : "x y z [mm]"));
    rl_gets(cmd, sizeof(cmd));
    if (cmd[0] == 'i') {
      imode += 1 - 2*imode;
      printf(" >> %s mode\n", (imode ? "ix iy iz" : "x y z"));
      continue;
    } else if (cmd[0] == 'q') {
      break;
    } else if (cmd[0] == 'h') {
      printf("%s", help);
      continue;
    } else if (cmd[0] == 'f') {
      sscanf(cmd+1, "%f", &f);
      if (f < 0 || f > 10.0*step) printf("Error: FWHM value %.2f is negative or too large.\n", f);
      else fwhm = f;
      printf(" >> FWHM = %.2f\n", fwhm);
      continue;
    } else if (cmd[0] == 't') {
      sscanf(cmd+1, "%f", &tau);
      printf(" >> tau = %.1f\n", tau);
      continue;
    } else {

      if (imode) {
        if (sscanf(cmd, "%d %d %d", &ix, &iy, &iz) < 3) {
          printf("can't decode line %s", cmd);
          continue;
        }
        printf(" >> x y z %3d %3d %3d -> %5.1f %5.1f %5.1f\n",
               ix, iy, iz, (ix-nx/2 + 0.5)*step, (iy-ny/2 + 0.5)*step, (iz + 0.5)*step);
      } else {
        if (sscanf(cmd, "%f %f %f", &x, &y, &z) < 3) {
          printf("can't decode line %s", cmd);
          continue;
        }
        ix = lrintf(x/step - 0.5) + nx/2;
        iy = lrintf(y/step - 0.5) + ny/2;
        iz = lrintf(z/step - 0.5);
        printf(" >> x y z %5.1f %5.1f %5.1f -> %3d %3d %3d\n",
               x, y, z, ix, iy, iz);
      }
      if (ix < 0 || ix >= nx || iy <= 0 || iy >= ny || iz <= 0 || iz >= nz ||
          (seg = seg_dir[ix][iy][iz]) < 0) {
        printf("point is outside crystal\n");
        continue;
      }
      printf("point is in segment %2d\n", seg);

      if (fwhm < 0.03) {  // use only one grid point
        pos = pos_dir[ix][iy][iz];
        fseek(fpi, pos*2, SEEK_SET);
        fread(ss, sizeof(short), 5, fpi);
        len  = ss[4];
        if (ix  != ss[0] || iy  != ss[1] || iz  != ss[2] || seg != ss[3] || len < 42 || len > 4000) {
          printf("ERROR reading signal header!\n");
          printf("ix,iy,iz, seg, len: = %2d %2d %2d, %2d, %4d  pos %lld\n", ix,iy,iz, seg, len, pos);
          continue;
        }
        fread(ss, sizeof(short), len-5, fpi);
        for (i=0; i<4096; i++) fs[i] = 0;
        j = 0;
        for (i=0; i<nsegs; i++) {
          n = ss[j++];
          // printf("iseg j n = %2d %4d %3d\n", i, j, n);
          for (k=0; k<n-1; k++)    fs[i*100 + k] = ss[j++] / 10.0;
          if (abs(ss[j-1]) < 5000) fs[i*100 + k] = 0;
          if (ss[j-1]    <= -5000) fs[i*100 + k] = -1000;
          if (ss[j-1]    >=  5000) fs[i*100 + k] =  1000;
          for (k++; k<100;  k++)    fs[i*100 + k] = fs[i*100 + k-1];
        }

      } else {       // sum several grid point signals
        float fact[40] = {0}, sum1 = 0, sum2 = 0, f3;
        float w  = step * 2.35482 / fwhm;
        int   ix1, iy1, iz1, ns;
        for (ns=0; ns<30 && ns<3.5*fwhm/step; ns++) {
          fact[ns] = exp(-ns*ns * w*w);
          sum1 += fact[ns];
        }
        sum1 = sum1*sum1*sum1;

        for (i=0; i<4096; i++) fs[i] = 0;
        for (ix1 = ix - ns + 1; ix1 < ix + ns; ix1++) {
          for (iy1 = iy - ns + 1; iy1 < iy + ns; iy1++) {
            for (iz1 = iz - ns + 1; iz1 < iz + ns; iz1++) {
              f3 = fact[lrint(abs(ix-ix1))] * fact[lrint(abs(iy-iy1))] * fact[lrint(abs(iz-iz1))];
              if (f3/sum1 < 0.0002) continue;
              if (ix1 < 0 || ix1 >= nx || iy1 <= 0 || iy1 >= ny || iz1 <= 0 || iz1 >= nz ||
                  (seg = seg_dir[ix1][iy1][iz1]) < 0) {
                continue;
              }
              pos = pos_dir[ix1][iy1][iz1];
              fseek(fpi, pos*2, SEEK_SET);
              fread(ss, sizeof(short), 5, fpi);
              len  = ss[4];
              if (ix1  != ss[0] || iy1  != ss[1] || iz1  != ss[2] || seg != ss[3] || len < 42 || len > 4000) {
                printf("ERROR reading signal header!\n");
                printf("ix1,iy1,iz1, seg, len: = %2d %2d %2d, %2d, %4d  pos %lld\n", ix1,iy1,iz1, seg, len, pos);
                continue;
              }
              fread(ss, sizeof(short), len-5, fpi);
              j = 0;
              for (i=0; i<nsegs; i++) {
                n = ss[j++];
                for (k=0; k<n-1; k++) fs[i*100 + k] += ss[j++] * f3 / 10.0;
                for (; k<100;  k++) {
                  if      (ss[j-1] <= -5000) fs[i*100 + k] += -1000.0*f3;
                  else if (ss[j-1] >=  5000) fs[i*100 + k] +=  1000.0*f3;
                }
              }
              sum2 += f3;  // CHECKME: move to *before* check on point in detector?
            }
          }
        }
        for (i=0; i<4096; i++) fs[i] /= sum2;
      }

      fpo = fopen("s.dat", "w");
      fwrite(fs, sizeof(fs), 1, fpo);
      fclose(fpo);
    }

  }

  fclose(fpi);
  
  return 0;
}
