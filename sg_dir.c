
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


int main(int argc, char **argv)
{

  FILE   *fpi, *fpo;
  int    seg, nsegs, time_steps, len, snum=0;
  int    ix, iy, iz, nx, ny, nz;
  float  x, y, z, radius, length, step;
  short  ss[4000] = {0};
  short  ***seg_dir;
  long long int ***pos_dir, pos=5;


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
  length = ((float)ss[4])/100.0;  // xtal length in mm

  printf(" %d segs, %d time steps\n"
         " xtal radius x length = %.2f x %.2f; step = %.3f\n",
         nsegs, time_steps, radius, length, step);

  /* malloc space for 3D dir arrays */
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

  while (fread(ss, sizeof(short), 5, fpi)) {
    int i, j = len;
    ix   = ss[0];
    iy   = ss[1];
    iz   = ss[2];
    seg  = ss[3];
    len  = ss[4];
    if (seg < -1 || seg > 37 || len < 42 || len > 4000) {
      printf("Error in seg or len! snum %8d: ix,iy,iz, seg, len: = %2d %2d %2d, %2d, %4d  pos %lld\n",
             snum, ix,iy,iz, seg, len, pos);
      printf("trying to recover backwards\n");
      for (i=0; i < 5; i++) ss[j+i] = ss[i];
      for (i=j-6; i > 5; i--) {
        if (ss[i] < -5000) break;
      }
      printf("i= %d  j= %d, ss[i] = %d\n", i, j, ss[i]);
      if (i < 6) {
        printf("trying to recover forwards\n");
        i = j;
        while (fread(ss, sizeof(short), 1, fpi) && ss[0] > -5000) i++;
        printf("i= %d  j= %d, ss[0] = %d\n", i, j, ss[0]);
        while (fread(ss, sizeof(short), 1, fpi) && ss[0] < -1000) i++;
        printf("i= %d  j= %d, ss[0] = %d\n", i, j, ss[0]);
        fread(ss+1, sizeof(short), 4, fpi);
        ix   = ss[0];
        iy   = ss[1];
        iz   = ss[2];
        seg  = ss[3];
        len  = ss[4];
      } else {
        ix   = ss[i+1];
        iy   = ss[i+2];
        iz   = ss[i+3];
        seg  = ss[i+4];
        len  = ss[1+5];
      }
      printf("NEW snum %8d: ix,iy,iz, seg, len: = %2d %2d %2d, %2d, %4d  pos %lld\n",
             snum, ix,iy,iz, seg, len, pos);
      if (seg < -1 || seg > 37 || len < 42 || len > 4000) {
        printf("Sorry, can't recover! Aborting...\n");
        break;
      }
    }
    fread(ss, sizeof(short), len-5, fpi);

    seg_dir[ix][iy][iz] = seg;
    pos_dir[ix][iy][iz] = pos;    
    pos += len;
    snum++;

    x = ix*step - radius + step/2.0;
    y = iy*step - radius + step/2.0;
    z = iz*step + step/2.0;
  
    if (snum%100000 == 0)
      printf("snum %9d: %4d shorts of data read; x, y, z = %6.2f %6.2f %6.2f  ->  seg %2d   pos %lld\n",
             snum, len, x, y, z, seg_dir[ix][iy][iz], pos_dir[ix][iy][iz]);
    if (len > 1000) {
      int maxl = 0, l, k = 0;
      for (int i=0; i<37; i++) {
        l = ss[k];
        k += l;
        if (maxl < l) maxl = l;
      }
      if (maxl > 50)
        printf("Long signal: snum %9d: maxl = %3d; x, y, z = %6.2f %6.2f %6.2f  ->  seg %2d\n",
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
      fwrite(seg_dir[ix][iy], sizeof(***seg_dir), nz, fpo);
    }
  }
  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      fwrite(pos_dir[ix][iy], sizeof(***pos_dir), nz, fpo);
    }
  }

  /*
  int seg2, iz2;
  for (int seg = 0; seg < 36; seg++) {
    for (ix = 0; ix < nx; ix++) {
      for (iy = 0; iy < ny; iy++) {
        for (iz = 0; iz < nz; iz++) {
          if (seg_dir[ix][iy][iz] == seg) {
            if ((iz2 = iz-1) >= 0 && (seg2 = seg_dir[ix][iy][iz2]) != seg && seg2 >= 0)
              printf("seg %2d   z-bdry %6.2f %6.2f %6.2f  z2, seg2 =  %6.2f, %2d\n",
                     seg, ix*step - radius, iy*step - radius, iz*step, iz2*step, seg2);
            if ((iz2 = iz+1) < nz && (seg2 = seg_dir[ix][iy][iz2]) != seg)
              printf("seg %2d   z+bdry %6.2f %6.2f %6.2f  z2, seg2 =  %6.2f, %2d\n",
                     seg, ix*step - radius, iy*step - radius, iz*step, iz2*step, seg2);
  */ /*
 (//(ix > 0    && seg_dir[ix-1][iy][iz] != seg) ||
  //(iy > 0    && seg_dir[ix][iy-1][iz] != seg) ||
    (iz > 0    && seg_dir[ix][iy][iz-1] != seg) ||
 // (ix < nx-1 && seg_dir[ix+1][iy][iz] != seg) ||
 // (iy < ny-1 && seg_dir[ix][iy+1][iz] != seg) ||
    (iz < nz-1 && seg_dir[ix][iy][iz+1] != seg))) {
 && seg_dir[ix-1][iy][iz]g >= 0
 && seg_dir[ix][iy-1][iz]g >= 0
 && seg_dir[ix][iy][iz-1]g >= 0
 && seg_dir[ix+1][iy][iz]g >= 0
 && seg_dir[ix][iy+1][iz]g >= 0
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
        if (seg_dir[ix][iy][iz] == -1) {
          if (seg_dir[ix-1][iy][iz] >= 0 &&
              seg_dir[ix+1][iy][iz] >= 0 &&
              seg_dir[ix][iy-1][iz] >= 0 &&
              seg_dir[ix][iy+1][iz] >= 0 &&
              seg_dir[ix][iy][iz-1] >= 0 &&
              seg_dir[ix][iy][iz+1] >= 0)
            printf("no signal at %6.2f %6.2f %6.2f\n",
                   ix*step - radius, iy*step - radius, iz*step);
        }
      }
    }
  }
  
  return 0;
}
