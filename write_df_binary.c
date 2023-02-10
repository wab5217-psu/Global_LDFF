#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <math.h>


int main(int argc, char *argv[]){

  int yr,mo,dy,hr,mt,sc;
  int ngrid,dcount;
  int jc;
  float bx,by,bz,v_sw,au,al;
  float center_lat,center_lon,vx,vy,mod_vx,mod_vy;
  int date_vec[6];
  float b_vec[6];
  float dvec[7];
    
  
  char out_file[20];

  
  /* read instruction file
     file should have date, time, lat and lon of area corners, radars to contribute
     filtering instructions
  */

  FILE *infp;
  FILE *outfp;

  if ((infp=fopen("vel_out","r")) == 0) {
    fprintf(stderr,"******FILE %s NOT FOUND******\n",argv[1]);
    exit(-1);
  }


  fscanf(infp,"%d %d %d %d %d %d",&yr,&mo,&dy,&hr,&mt,&sc);
  
  sprintf(out_file,"%d%02d%02d.gdf_vel",yr,mo,dy);
  fprintf(stderr,"%s\n",out_file);
  outfp=fopen(out_file,"wb");
  
  rewind(infp);
  while( fscanf(infp,"%d %d %d %d %d %d",&yr,&mo,&dy,&hr,&mt,&sc) == 6){

    date_vec[0]=yr;
    date_vec[1]=mo;
    date_vec[2]=dy;
    date_vec[3]=hr;
    date_vec[4]=mt;
    date_vec[5]=sc;
    fwrite(&date_vec,1,sizeof(date_vec),outfp);
    
    fscanf(infp,"%f %f %f %f %f %f",&bx,&by,&bz,&v_sw,&au,&al);
    
    b_vec[0]=bx;
    b_vec[1]=by;
    b_vec[2]=bz;
    b_vec[3]=v_sw;
    b_vec[4]=au;
    b_vec[5]=al;
    fwrite(&b_vec,1,sizeof(b_vec),outfp);
    
    fscanf(infp,"%d",&ngrid);
    fwrite(&ngrid,1,sizeof(int),outfp);
    
    for( jc=0; jc<ngrid; jc++){
      fscanf(infp,"%f %f %f %f %f %f %d",&center_lat,&center_lon,&vx,&vy,&mod_vx,&mod_vy,&dcount);
      dvec[0]=center_lat;
      dvec[1]=center_lon;
      dvec[2]=vx;
      dvec[3]=vy;
      dvec[4]=mod_vx;
      dvec[5]=mod_vy;
      dvec[6]=(float)dcount;
      fwrite(&dvec,1,sizeof(dvec),outfp);
      
    }
  }
	  

  fclose(outfp);
  fclose(infp);
}  
