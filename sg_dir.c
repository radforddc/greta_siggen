
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


int main(int argc, char **argv)
{

  FILE   *fpi, *fpo;
  int    seg, nsegs, time_steps, len, snum=0;
  int    ix, iy, iz, nx, ny, nz, pos=5;
  float  x, y, z, radius, length, step;
  short  ss[4000] = {0};
  struct {
    int   seg, snum, pos;
  } ***dir;
    
  
  if (argc < 2) {
    printf("Usage: %s <input_sgrid.mat_file_name>\n", argv[0]);
    return 1;
  }

  /* open input file */
  fpi = fopen(argv[1], "r");
  if (fpi == NULL){
    fprintf(stderr,"Error! Unable to open input file argv[1]\n");
    return 0;
  }
  /* open output file */
  fpo = fopen("sgrid.dir", "w");
  if (fpo == NULL){
    fprintf(stderr,"Error! Unable to open output file sgrid.dir\n");
    return 0;
  }

  fread(ss, sizeof(short), 5, fpi);
  nsegs = ss[0];
  time_steps = ss[1];
  step   = ((float)ss[2])/100.0;  // grid step size in mm
  radius = ((float)ss[3])/100.0;  // xtal radius in mm
  length = ((float)ss[4])/100.0;  // xtal lengthin mm

  printf(" %d segs, %d time steps\n"
         " xtal radius x length = %.2f x %.2f; step = %.3f\n",
         nsegs, time_steps, radius, length, step);

  /* malloc space for 3D dir struct array */
  nx = ny = 2.0*radius/step + 1.5;
  nz = length/step + 1.5;
  printf("nx, ny, nz = %d, %d, %d\n\n", nx, ny, nz);

  if ((dir = malloc(nx*sizeof(*dir))) == NULL) {
    printf("malloc1 failed\n");
    return 1;
  }
  for (ix = 0; ix < nx; ix++) {
    if ((dir[ix] = malloc(ny*sizeof(*dir[ix]))) == NULL) {
      printf("malloc2 failed; ix = %d\n", ix);
      return 1;
    }
    for (iy = 0; iy < ny; iy++) {
      if ((dir[ix][iy] = malloc(nz*sizeof(*dir[ix][iy]))) == NULL) {
        printf("malloc3 failed; ix, iy = %d, %d\n", ix, iy);
        return 1;
      }
      for (iz = 0; iz < nz; iz++) {
        dir[ix][iy][iz].seg  = -1;
        dir[ix][iy][iz].snum = -1;
        dir[ix][iy][iz].pos  = -1;
      }
    }
  }

  while (fread(ss, sizeof(short), 5, fpi)) {
    ix   = ss[0];
    iy   = ss[1];
    iz   = ss[2];
    seg  = ss[3];
    len  = ss[4];
    // printf("snum %8d: ix,iy,iz, seg, len: = %2d %2d %2d, %2d, %4d\n", snum, ix,iy,iz, seg, len);
    fread(ss, sizeof(short), len-5, fpi);

    dir[ix][iy][iz].seg  = seg;
    dir[ix][iy][iz].snum = snum;
    dir[ix][iy][iz].pos  = pos;    
    pos += len;
    snum++;

    x = ix*step - radius;
    y = iy*step - radius;
    z = iz*step;
  
    if (snum%20000 == 0)
      printf("snum %9d: %4d shorts of data read; x, y, x = %7.3f %7.3f %7.3f  ->  seg %2d\n",
             snum, len, x, y, z, seg);
    if (len > 1000) {
      int maxl = 0, l, k = 0;
      for (int i=0; i<37; i++) {
        l = ss[k];
        k += l;
        if (maxl < l) maxl = l;
      }
      if (maxl > 50)
        printf("Long signal: snum %9d: maxl = %3d; x, y, x = %7.3f %7.3f %7.3f  ->  seg %2d\n",
               snum, maxl-1, x, y, z, seg);
    }
  }
  
  printf("\n Reading done. %d signals read.\n\n", snum);

  fwrite(&nsegs     , sizeof(nsegs     ), 1, fpo);
  fwrite(&time_steps, sizeof(time_steps), 1, fpo);
  fwrite(&step      , sizeof(step      ), 1, fpo);
  fwrite(&radius    , sizeof(radius    ), 1, fpo);
  fwrite(&length    , sizeof(length    ), 1, fpo);

  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      fwrite(dir[ix][iy], sizeof(***dir), nz, fpo);
    }
  }

  /*
  int seg2, iz2;
  for (int seg = 0; seg < 36; seg++) {
    for (ix = 0; ix < nx; ix++) {
      for (iy = 0; iy < ny; iy++) {
        for (iz = 0; iz < nz; iz++) {
          if (dir[ix][iy][iz].seg == seg) {
            if ((iz2 = iz-1) >= 0 && (seg2 = dir[ix][iy][iz2].seg) != seg && seg2 >= 0)
              printf("seg %2d   z-bdry %6.2f %6.2f %6.2f  z2, seg2 =  %6.2f, %2d\n",
                     seg, ix*step - radius, iy*step - radius, iz*step, iz2*step, seg2);
            if ((iz2 = iz+1) < nz && (seg2 = dir[ix][iy][iz2].seg) != seg)
              printf("seg %2d   z+bdry %6.2f %6.2f %6.2f  z2, seg2 =  %6.2f, %2d\n",
                     seg, ix*step - radius, iy*step - radius, iz*step, iz2*step, seg2);
  */ /*
 (//(ix > 0  && dir[ix-1][iy][iz].seg != seg) ||
 //(iy > 0  && dir[ix][iy-1][iz].seg != seg) ||
 (iz > 0  && dir[ix][iy][iz-1].seg != seg) ||
 //(ix < nx-1 && dir[ix+1][iy][iz].seg != seg) ||
 //(iy < ny-1 && dir[ix][iy+1][iz].seg != seg) ||
 (iz < nz-1 && dir[ix][iy][iz+1].seg != seg))) {
&& dir[ix-1][iy][iz].seg >= 0
 && dir[ix][iy-1][iz].seg >= 0
 && dir[ix][iy][iz-1].seg >= 0
 && dir[ix+1][iy][iz].seg >= 0
 && dir[ix][iy+1][iz].seg >= 0
                */                       
  /*          }
        }
      }
    }
  }
  */

  for (ix = 1; ix < nx-1; ix++) {
    for (iy = 1; iy < ny-1; iy++) {
      for (iz = 1; iz < nz-1; iz++) {
        if (dir[ix][iy][iz].seg == -1) {
          if (dir[ix-1][iy][iz].seg >= 0 &&
              dir[ix+1][iy][iz].seg >= 0 &&
              dir[ix][iy-1][iz].seg >= 0 &&
              dir[ix][iy+1][iz].seg >= 0 &&
              dir[ix][iy][iz-1].seg >= 0 &&
              dir[ix][iy][iz+1].seg >= 0)
            printf("no signal at %6.2f %6.2f %6.2f\n",
                   ix*step - radius, iy*step - radius, iz*step);
        }
      }
    }
  }
  
  return 0;
}
