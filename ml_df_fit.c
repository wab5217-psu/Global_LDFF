/* #include <stdlib.h> */
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include <zlib.h>
#include "rtypes.h"
#include "dmap.h"
#include "option.h"
#include "rtime.h"
#include "radar.h"
#include "rprm.h"
#include "rpos.h"
#include "fitdata.h"
#include "cfitdata.h"
#include "scandata.h"
#include "fitread.h"
#include "fitscan.h"
#include "fitindex.h"
#include "fitseek.h"
#include "rtypes.h"
#include "dmap.h"
#include "invmag.h"
#include "griddata.h"
#include "gridread.h"

#include "ml_df.h"

#include "cnvgrid.h"
#include "cnvmap.h"
#include "cnvmapindex.h"
#include "cnvmapseek.h"
#include "cnvmapread.h"
#include "cnvmapsolve.h"
#include "aacgmlib_v2.h"
#include "aacgm.h"
#include "mlt_v2.h"
#include "igrflib.h"
#include "terminator.h"

#include <gsl/gsl_math.h> 
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#define LENGTH(x,y) sqrt(x*x+y*y)

#define INERTIAL 0

#define C 299792458.0 
#define PI 3.14159265359
#define LRE 6371.
#define MIN_RANGE 600
#define MAX_RANGE 3000
#define SECS_P_DAY 86400 /* 24*60*60 */
#define MAX(q,p) (((q)>(p))?(q):(p))
#define sind(x) (sin(fmod((x),360)*PI/180))
#define cosd(x) (cos(fmod((x),360)*PI/180))
#define tand(x) (tan(fmod((x),360)*PI/180))
#define MIN_ERR 200.0
#define ERR_SCALE 1.0
#define MIN_ML_ERR 100.0
#define MAX_ML_COV 1000.0
#define DIV_ERR 0.1
#define MAX_V 3000.0
#define MAX_V_ERR 500.0
#define MIN_V 30.
#define MIN_COUNT 4
#define MAX_BEAMS 24
#define MIN_LAT 55
#define MAX_LAT 90

#define ML_LAT_0 55
#define ML_DLAT 2
#define ML_DLON 15
#define ML_NLON 24

#define BAD_VALUE -99999.9
#define BAD_INT -9999

#define WRITE_COEF 0

#define MAX_DATA_LAT 89

double dlat;
double dlon;
double dlat_ng;
double dlon_ng;
double min_lon_ng;
double max_lon_ng;
long start_time;
long end_time;
int avg_ival;
int smooth;
double MODEL_SCALE;
char *radar_list[30];
int nrad;
int nbp;
B_POINT* bp;
int ngrid;
int nedge;
M_ARRAY* kazm_array;
double *start_lon;
CELL* grid;
NEIGHBOR* neighbors;
COROTATION* corotation_velocity;
int *edge; 
M_ARRAY* edge_array;
double* edge_data;
double* los_data;
double* in_los_data;
double* los_kazm;
double* los_lats;
double* los_lons;
double* los_err;
M_ARRAY* div_array;
M_ARRAY* smooth_array;
int mdays[]={31,59,90,120,151,181,212,243,273,304,334};
FILE *grd_file;
char ml_file[128];
char grid_file[128];
char hemisphere[16]="north";

int dayofweek(int d, int m, int y)
{
    static int t[] = { 0, 3, 2, 5, 0, 3, 5, 1, 4, 6, 2, 4 };
    y -= m < 3;
    return ( y + y/4 - y/100 + y/400 + t[m-1] + d) % 7;
}
 
char *choppy( char *s )
{
  char *n = malloc( strlen( s ? s : "\n" ) );
  if( s )
    strcpy( n, s );
  if( s[strlen(s)-1] =='\n')
    n[strlen(n)-1]='\0';
  return n;
}


void parse_instructions(FILE *fp)
{
  char *line=NULL;
  size_t len=0;
  char *param=NULL;
  char *token;

  nrad=0;
  nbp=0;
  while( getline(&line,&len,fp)!= EOF ){
    if(line[0]=='#' || line[0]==' ')continue;
    param=strtok(line," ");
    if(strcmp(param,"boundary_point")==0)
      {
	nbp++;
	bp=realloc(bp,nbp*sizeof(B_POINT));
	sscanf(strtok(NULL," "),"%lf",&bp[nbp-1].lat);
	sscanf(strtok(NULL," "),"%lf",&bp[nbp-1].lon);
	if( bp[nbp-1].lon<0 )bp[nbp-1].lon+=360;
      }
   else if(strcmp(param,"lat_dl")==0)
     {
	sscanf(strtok(NULL," "),"%lf",&dlat);
     }
   else if(strcmp(param,"lon_dl")==0)
     {
	sscanf(strtok(NULL," "),"%lf",&dlon);
     }
   else if(strcmp(param,"ng_lat_dl")==0)
     {
	sscanf(strtok(NULL," "),"%lf",&dlat_ng);
     }
   else if(strcmp(param,"ng_lon_dl")==0)
     {
	sscanf(strtok(NULL," "),"%lf",&dlon_ng);
     }
   else if(strcmp(param,"start_time")==0)
     {
	sscanf(strtok(NULL," "),"%ld",&start_time);
     }
   else if(strcmp(param,"end_time")==0)
     {
	sscanf(strtok(NULL," "),"%ld",&end_time);
     }
   else if(strcmp(param,"avg_interval")==0)
     {
	sscanf(strtok(NULL," "),"%d",&avg_ival);
     }
   else if(strcmp(param,"ml_file")==0)
     {
	sscanf(strtok(NULL," "),"%s",ml_file);
     }
   else if(strcmp(param,"grid_file")==0)
     {
	sscanf(strtok(NULL," "),"%s",grid_file);
     }
   else if(strcmp(param,"hemisphere")==0)
     {
	sscanf(strtok(NULL," "),"%s",hemisphere);
     }
   else if(strcmp(param,"smooth")==0)
     {
	sscanf(strtok(NULL," "),"%d",&smooth);
     }
   else if(strcmp(param,"model_scale")==0)
     {
	sscanf(strtok(NULL," "),"%lf",&MODEL_SCALE);
     }
   else if(strcmp(param,"radar_list")==0)
     {
       while( (token=strtok(NULL," ")) != NULL)
	 {
	   radar_list[nrad]=choppy(token);
	   nrad++;
	 }
     }
 }
  if(line)
    free(line);
}


int poly_in_out(double lat, double lon )
{
  int j,c=0;
  double lonh,lon1h,lon2h;
  double lat1,lon1;
  double lat2,lon2;
  if( nbp==0 )return 0;
  
  lat1=bp[0].lat;
  lon1=bp[0].lon;
  while( lon1<0 )lon1+=360.0;
  c=0;
  
  for( j=1; j<=nbp; j++)
    {
      lat2=bp[j % nbp].lat;
      lon2=bp[j % nbp].lon;
      while( lon2<0 )lon2+=360.0;

      if( lon1-lon2>180 ){ lon2h=lon2+360; }else{ lon2h=lon2; }
      if( lon2-lon1>180 ){ lon1h=lon1+360; }else{ lon1h=lon1; }
      if( (lon1h-lon>180)||(lon2h-lon>180) ){ lonh=lon+360; }else{lonh=lon;} 
      
      if( ((lonh>lon1h)&&(lonh<lon2h))||((lonh<lon1h)&&(lonh>lon2h)) ){
	if( (lat>lat1)||(lat>lat2) )
	  c=!c;
      }
      lat1=lat2;
      lon1=lon2;
    }
  return c;
}


/* determines if a point is inside a grid cell returns "1" for inside "0" for outside*/

int cell_in_out(double lat, double lon, CELL cell )
{
  double circ=360.0;
  double in_lon;
  
  in_lon=lon+.001;
  if(in_lon == circ) in_lon=359.98;
  if(in_lon == 0.0) in_lon=.02;

  if( (fabs(lat)<fabs(cell.lat[0])) || (fabs(lat)>fabs(cell.lat[2])) )return(0);
    
  if( (in_lon>=cell.lon[0])&&(in_lon<cell.lon[1]) )
    return(1);


  if( ((in_lon-circ)>=cell.lon[0])&&((in_lon-circ)<cell.lon[1]) )
    return(1);

  if( ((in_lon+circ)>=cell.lon[0])&&((in_lon+circ)<cell.lon[1]) )
    return(1);
  
  return(0);
}  

void make_box(double lat, double lon, double ldlat, double ldlon, CELL *box)
{
  box->lat[0]=lat-ldlat/2;
  box->lat[1]=lat-ldlat/2;
  box->lat[2]=lat+ldlat/2;
  box->lat[3]=lat+ldlat/2;
  box->lon[0]=lon-ldlon/2;
  box->lon[1]=lon+ldlon/2;
  box->lon[2]=lon+ldlon/2;
  box->lon[3]=lon-ldlon/2;
}  

double get_random() { return (double)rand() / (double)RAND_MAX; }

void get_grid_size(double min_lat, double max_lat){
  double lat,lon;
  double min_lon=0;
  double max_lon=360;

  double dlon_l,dlon_l_ng;
  double start_ng_lon, end_ng_lon;
  int nlon,nlon_ng,nlon_m,nlon_p;
  int nlat=(max_lat-min_lat)/dlat;
  int jlat,jlon,jlat_ng;

  int lat_ratio=(int)(dlat/dlat_ng);
  
  nlon_ng=0;
  nlon_m=0;
  nlon_p=0;
  
  ngrid=0;

  for( jlat=0; jlat<nlat; jlat++ )
    {
      lat=min_lat+((double)jlat+.5)*dlat;      
      dlon_l=dlat/cosd(lat);
      nlon=(int)((max_lon-min_lon)/dlon_l);
      dlon_l=360/(double)nlon;
      lon=start_lon[jlat];
      start_ng_lon=BAD_VALUE;
      end_ng_lon=BAD_VALUE;
      
      dlon_l_ng=dlat_ng/cosd(lat);

      for( jlon=0; jlon<nlon; jlon++)
	{
	  lon=start_lon[jlat]+((double)jlon+.5)*dlon_l;
	  if( poly_in_out( lat, lon) && (start_ng_lon == BAD_VALUE)){
	    start_ng_lon=lon-dlon_l/2;
	    nlon_m=jlon;
	  }
	  if( poly_in_out( lat, lon)){
	    end_ng_lon=lon+dlon_l/2;
	    nlon_p=jlon;
	  }
	}

      nlon_ng=(int)((end_ng_lon-start_ng_lon)/dlon_l_ng);
      dlon_l_ng=(end_ng_lon-start_ng_lon)/(double)nlon_ng;

      if( start_ng_lon != BAD_VALUE ){
	for( jlon=0; jlon<nlon_m; jlon++)
	  {
	    lon=start_lon[jlat]+((double)jlon+.5)*dlon_l;
	    if( lon > 360 )lon-=360;	    
	    ngrid++;
	  }	

	for( jlat_ng=0; jlat_ng<lat_ratio; jlat_ng++ ){
	  for( jlon=0; jlon<nlon_ng; jlon++)
	    {
	      lon=start_ng_lon+((double)jlon+.5)*dlon_l_ng;
	      ngrid++;
	    }		
	}

	for( jlon=nlon_p+1; jlon<nlon; jlon++)
	  {
	    lon=start_lon[jlat]+((double)jlon+.5)*dlon_l;
	    if( lon > 360 )lon-=360;	    
	    ngrid++;
	  }	
	
      }else{	
	for( jlon=0; jlon<nlon; jlon++)
	  {
	    lon=start_lon[jlat]+((double)jlon+.5)*dlon_l;
	    if( lon > 360 )lon-=360;	    
	    ngrid++;
	  }	
      }
    }  
}


void make_grid(double min_lat, double max_lat){
  double lat,lon,lat_ng;
  double min_lon=0;
  double max_lon=360;

  double dlon_l,dlon_l_ng;
  double start_ng_lon, end_ng_lon;
  int nlon,nlon_ng,nlon_m,nlon_p;
  int nlat=(max_lat-min_lat)/dlat;
  int jlat,jlon,jlat_ng,jgrid;
  int j;
  
  int lat_ratio=(int)(dlat/dlat_ng);

  nlon_ng=0;
  nlon_m=0;
  nlon_p=0;
  
  /* srand(time(NULL)); // randomize seed */
  jgrid=0;

  dlon_l=dlat/cosd(min_lat+dlat/2);
  nedge=(max_lon-min_lon)/dlon_l;
  
  for( jlat=0; jlat<nlat; jlat++ )
    {
      lat=min_lat+((double)jlat+.5)*dlat;      
      dlon_l=dlat/cosd(lat);

      nlon=(int)((max_lon-min_lon)/dlon_l);
      dlon_l=360/(double)nlon;
      lon=start_lon[jlat];
      start_ng_lon=BAD_VALUE;
      end_ng_lon=BAD_VALUE;

      dlon_l_ng=dlat_ng/cosd(lat);

      for( jlon=0; jlon<nlon; jlon++)
	{
	  lon=start_lon[jlat]+((double)jlon+.5)*dlon_l;
	  if( poly_in_out( lat, lon) && (start_ng_lon == BAD_VALUE)){
	    start_ng_lon=lon-dlon_l/2;
	    nlon_m=jlon;
	  }
	  if( poly_in_out( lat, lon)){
	    end_ng_lon=lon+dlon_l/2;
	    nlon_p=jlon;
	  }
	}

      nlon_ng=(int)((end_ng_lon-start_ng_lon)/dlon_l_ng);
      dlon_l_ng=(end_ng_lon-start_ng_lon)/(double)nlon_ng;
      

      if( start_ng_lon != BAD_VALUE ){
	dlon_l=(start_ng_lon-start_lon[jlat])/(double)nlon_m;
	for( jlon=0; jlon<nlon_m; jlon++)
	  {
	    lon=start_lon[jlat]+((double)jlon+.5)*dlon_l;
	    /* if( lon > 360 )lon-=360;	     */
	    make_box(lat,lon,dlat,dlon_l,&grid[jgrid]);
	    grid[jgrid].lat_indx=jlat;
	    grid[jgrid].lon_indx=jlon;
	    grid[jgrid].center_lat=lat;
	    grid[jgrid].center_lon=lon;
	    jgrid++;
	  }	

	for( jlat_ng=0; jlat_ng<lat_ratio; jlat_ng++ ){
	  lat_ng=lat-dlat/2+((double)jlat_ng+.5)*dlat_ng;	  
	  for( jlon=0; jlon<nlon_ng; jlon++)
	    {
	      lon=start_ng_lon+((double)jlon+.5)*dlon_l_ng;
	      make_box(lat_ng,lon,dlat_ng,dlon_l_ng,&grid[jgrid]);
	      grid[jgrid].lat_indx=jlat;
	      grid[jgrid].lon_indx=jlon;
	      grid[jgrid].center_lat=lat_ng;
	      grid[jgrid].center_lon=lon;
	      jgrid++;
	    }		
	}

	dlon_l=(360.0+start_lon[jlat]-end_ng_lon)/(nlon-nlon_p-1);	  
	for( jlon=nlon_p+1; jlon<nlon; jlon++)	
	  {
	    lon=start_lon[jlat]+((double)jlon+.5)*dlon_l;
	    /* if( lon > 360 )lon-=360;	     */
	    make_box(lat,lon,dlat,dlon_l,&grid[jgrid]);
	    grid[jgrid].lat_indx=jlat;
	    grid[jgrid].lon_indx=jlon;
	    grid[jgrid].center_lat=lat;
	    grid[jgrid].center_lon=lon;
	    jgrid++;
	  }	
	
      }else{	
	for( jlon=0; jlon<nlon; jlon++)
	  {
	    lon=start_lon[jlat]+((double)jlon+.5)*dlon_l;
	    /* if( lon > 360 )lon-=360;	     */
	    make_box(lat,lon,dlat,dlon_l,&grid[jgrid]);
	    grid[jgrid].lat_indx=jlat;
	    grid[jgrid].lon_indx=jlon;
	    grid[jgrid].center_lat=lat;
	    grid[jgrid].center_lon=lon;
	    jgrid++;
	  }	
      }
    }

  if( strcmp(hemisphere,"south")==0 ){
    for( j=0; j<jgrid; j++){
      grid[j].center_lat*=-1;
      grid[j].lat[0]*=-1;
      grid[j].lat[1]*=-1;
      grid[j].lat[2]*=-1;
      grid[j].lat[3]*=-1;
    }
  }
  
}

void find_neighbors(){
  int jg,jjg,ig;
  int jc;
  double lat,lon;
  double dlat_l;
  double h_factor=1;
  if( strcmp(hemisphere,"south")==0 ){
    h_factor=-1;
  }

  for( jg=0; jg<ngrid; jg++ ){
    for( jc=0; jc<3; jc++ )neighbors[jg].lower[jc]=BAD_INT;
    for( jc=0; jc<3; jc++ )neighbors[jg].upper[jc]=BAD_INT;
    for( jc=0; jc<2; jc++ )neighbors[jg].left[jc]=BAD_INT;
    for( jc=0; jc<2; jc++ )neighbors[jg].right[jc]=BAD_INT;
    dlat_l=grid[jg].lat[2]-grid[jg].lat[0];
    
    /* find adjacent right cells */
    lat=grid[jg].lat[1]+.1*h_factor;
    lon=grid[jg].lon[1]+.1;
    for(ig=0; ig<ngrid; ig++)if(cell_in_out(lat,lon,grid[ig])){
	neighbors[jg].right[0]=ig;
	break;
      }
    
    if(neighbors[jg].right[0]<0)fprintf(stderr,"bad right neighbor %d %d %6.2f %6.2f\n",jg,neighbors[jg].right[0],grid[jg].center_lat,grid[jg].center_lon);
    if( dlat_l > (grid[neighbors[jg].right[0]].lat[2]-grid[neighbors[jg].right[0]].lat[0]) ){
      lat=grid[jg].lat[2]-.1*h_factor;
      lon=grid[jg].lon[1]+.1;
      for(ig=0; ig<ngrid; ig++)if(cell_in_out(lat,lon,grid[ig])){
	  neighbors[jg].right[1]=ig;
	  break;
	}
    }    
    
    /* find adjacent left cells */ 
    lat=grid[jg].lat[0]+.1*h_factor;
    lon=grid[jg].lon[0]-.1;
    for(ig=0; ig<ngrid; ig++)if(cell_in_out(lat,lon,grid[ig])){
	neighbors[jg].left[0]=ig;
	break;
      }

    if(neighbors[jg].left[0]<0)fprintf(stderr,"bad left neighbor %d %d %6.2f %6.2f\n",jg,neighbors[jg].left[0],grid[jg].center_lat,grid[jg].center_lon);
    if( dlat_l > (grid[neighbors[jg].left[0]].lat[2]-grid[neighbors[jg].left[0]].lat[0]) ){
      lat=grid[jg].lat[3]-.1*h_factor;
      lon=grid[jg].lon[0]-.1;
      for(ig=0; ig<ngrid; ig++)if(cell_in_out(lat,lon,grid[ig])){
	  neighbors[jg].left[1]=ig;
	  break;
	}
    }
  }

  for( jg=nedge; jg<ngrid; jg++ ){

    /* find lower neighbor of lower-left corner: */
    lat=grid[jg].lat[0]-.1*h_factor;
    lon=grid[jg].lon[0]+.001;
    for(ig=0; ig<ngrid; ig++)if(cell_in_out(lat,lon,grid[ig])){
	neighbors[jg].lower[0]=ig;
	break;
      }
    
    /* step in longitude to find all neighbor cells below current cell */
    ig=neighbors[jg].lower[0];
        
    if( (ig<0) || (ig>ngrid-1) )continue;
    
    if( neighbors[ig].right[1] != BAD_INT){
      jjg=neighbors[ig].right[1];
    }else{
      jjg=neighbors[ig].right[0];
    }
    lat=grid[jjg].lat[3]+.1*h_factor;
    lon=grid[jjg].lon[3]+.001;
    if(cell_in_out(lat,lon,grid[jg])){
      neighbors[jg].lower[1]=jjg;
    }

    ig=jjg;
    if( neighbors[ig].right[1] != BAD_INT){
      jjg=neighbors[ig].right[1];
    }else{
      jjg=neighbors[ig].right[0];
    }
    lat=grid[jjg].lat[3]+.1*h_factor;
    lon=grid[jjg].lon[3]+.001;
    if(cell_in_out(lat,lon,grid[jg])){
      neighbors[jg].lower[2]=jjg;
    }
  }
}


void make_div_ar(){
  int igr,inr0,inr1,inr2;
  double fr0,fr1,fr2,chk;
  int idv;
  int ind;
  double dx,dy,dlon_l;
  double dlon0;
  FILE *divfile;

  divfile=fopen("divergence_info.dat","w");
  

  idv=0;
  fprintf(stderr,"make_div_ar\n");
  div_array=calloc(ngrid-nedge,sizeof(struct mod_array));
  /* first nedge points are edge points, so start calculating divergence after those */
  for( igr=nedge; igr<ngrid; igr++){

    /* calculate d/dx term as vx_cell - vx_cell_left
     if edge of nested grid area, use half from lower half from upper */
    
    div_array[idv].coef=calloc(2*ngrid,sizeof(double));

    /* calculate d/dx term as vx_cell_right - vx_cell
     if edge of nested grid area, use half from lower half from upper */    

    inr0=neighbors[igr].right[0];
    inr1=neighbors[igr].right[1];

    dlon_l=grid[igr].lon[1]-grid[igr].lon[0];
      
    dx=fabs(dlon_l*DTOR*LRE*cosd(grid[igr].center_lat));
    div_array[idv].coef[igr]=-1/dx;
    if( inr1==BAD_INT){
      div_array[idv].coef[inr0]=1/dx;
    }else{
      div_array[idv].coef[inr0]=1/(2*dx);
      div_array[idv].coef[inr1]=1/(2*dx);
    }
    fprintf(divfile,"%d %d %d %f ",igr,inr0,inr1,dx);
    
    
    /* calculate d/dy term as vy_cell - vy_cell_lower */

    dy=fabs((grid[igr].lat[2]-grid[igr].lat[1])*DTOR*LRE);    
    div_array[idv].coef[igr+ngrid]=1/dy;

    inr0=neighbors[igr].lower[0];
    inr1=neighbors[igr].lower[1];
    inr2=neighbors[igr].lower[2];
    fr0=0;
    fr1=0;
    fr2=0;
    if( inr0 == BAD_INT ){ 
      div_array[idv].coef[inr1+ngrid]=-1/dy;
    }else if(inr1 == BAD_INT){
      fr0=1;
      div_array[idv].coef[inr0+ngrid]=-1/dy;
    } else if( inr2==BAD_INT) { 
      dlon0=(grid[inr0].lon[2]-grid[igr].lon[0]);
      if( dlon0>360 )dlon0-=360;
      fr0=dlon0/(grid[igr].lon[1]-grid[igr].lon[0]);
      fr1=1-fr0;
      div_array[idv].coef[inr0+ngrid]=-fr0/dy;
      div_array[idv].coef[inr1+ngrid]=-fr1/dy;
    }else{
      dlon0=(grid[inr0].lon[2]-grid[igr].lon[0]);
      if( dlon0>360 )dlon0-=360;
      fr0=dlon0/(grid[igr].lon[1]-grid[igr].lon[0]);
      
      dlon0=(grid[inr1].lon[2]-grid[inr1].lon[3]);
      if( dlon0>360 )dlon0-=360;
      fr1=dlon0/(grid[igr].lon[1]-grid[igr].lon[0]);
      
      fr2=1-(fr0+fr1);
      div_array[idv].coef[inr0+ngrid]=-fr0/dy;
      div_array[idv].coef[inr1+ngrid]=-fr1/dy;
      div_array[idv].coef[inr2+ngrid]=-fr2/dy;      
    }
    fprintf(divfile,"%d %d %d %f %f %f %f\n",inr0,inr1,inr2,dy,fr0,fr1,fr2);
    if( fabs(fr0+fr1+fr2-1)>.001 )fprintf(stderr,"make_div_ar: %d %d %7.3f %7.3f %7.3f %7.3f %7.3f %d %d %d\n",igr,idv,dx,dy,fr0,fr1,fr2,inr0,inr1,inr2);

    chk=0;
    for(ind=0; ind<2*ngrid; ind++){chk+=div_array[idv].coef[ind];}
    if( fabs(chk)>.01 )fprintf(stderr,"make_div_ar: %d check: %lf\n",idv,chk);
	
    idv++;
  }
  fprintf(stderr,"make_div_ar exit\n");
  fclose(divfile);
}


struct tm *parse_date_str( long t_i)
{
  struct tm *t_o;
  char *tz;

  tz = getenv("TZ");
  setenv("TZ", "", 1);
  tzset();

  t_o=malloc(sizeof(struct tm));
  t_o->tm_year=(int)(t_i/1e8);
  t_o->tm_mon=(int)((t_i-1e8*t_o->tm_year)/1e6);
  t_o->tm_mday=(int)((t_i-1e8*t_o->tm_year-1e6*t_o->tm_mon)/1e4);
  t_o->tm_hour=(int)((t_i-1e8*t_o->tm_year-1e6*t_o->tm_mon-t_o->tm_mday*1e4)/1e2);
  t_o->tm_min=(int)(t_i-1e8*t_o->tm_year-1e6*t_o->tm_mon-t_o->tm_mday*1e4-t_o->tm_hour*1e2);
  t_o->tm_yday=mdays[t_o->tm_mon-1]+t_o->tm_mday;
  t_o->tm_wday=dayofweek(t_o->tm_mday,t_o->tm_mon,t_o->tm_year);
  if( IS_LEAPYEAR(t_o->tm_year) && t_o->tm_mon>2) t_o->tm_yday++;
  t_o->tm_sec=0;
  t_o->tm_year-=1900;
  t_o->tm_mon-=1;
  t_o->tm_isdst=0;
  return t_o;
}

time_t fname_to_time(char *fname)
{
  int datev;
  char datestr[10]="0";
  char hr_str[3]="0";
  char mn_str[3]="0";
  struct tm f_tm;
  char *tz;
  
  tz = getenv("TZ");
  setenv("TZ", "", 1);
  tzset();
  
  memset(datestr, '\0', sizeof datestr);
  memset(hr_str, '\0', sizeof hr_str);
  memset(mn_str, '\0', sizeof mn_str);

  strncpy(datestr,fname,8);
  strncpy(hr_str,fname+9,2);
  strncpy(mn_str,fname+11,2);
  datev=atoi(datestr);
  f_tm.tm_year=datev/10000;
  f_tm.tm_mon=(datev-10000*f_tm.tm_year)/100;
  f_tm.tm_mday=datev-10000*f_tm.tm_year-100*f_tm.tm_mon;
  f_tm.tm_yday=mdays[f_tm.tm_mon-1]+f_tm.tm_mday;
  if( IS_LEAPYEAR(f_tm.tm_year) && f_tm.tm_mon>2) f_tm.tm_yday++;
  f_tm.tm_hour=atoi(hr_str);
  f_tm.tm_min=atoi(mn_str);
  f_tm.tm_sec=0;
  f_tm.tm_isdst=0;
  f_tm.tm_wday=dayofweek(f_tm.tm_mday,f_tm.tm_mon,f_tm.tm_year);
  f_tm.tm_year-=1900; /* unix epoch year correction */
  f_tm.tm_mon-=1;     /* unix epoch month 0 to 11 */

  return mktime(&f_tm);
}

int split_dateline(char* date_line, MLDstr *mlDstr){

  char * token0 = strtok(date_line," ");
  char * token1 = strtok(NULL," ");

  sscanf(token0,"%d-%d-%d",&mlDstr->yr, &mlDstr->mo, &mlDstr->dy);
  sscanf(token1,"%d:%d:%d",&mlDstr->hr, &mlDstr->mt, &mlDstr->sc);

  return(0);
}


int read_ml_record(FILE *fp, MLDstr *mlDstr){

  double n_rmse,e_rmse;
  char date_line[256];
  int stat;
  int npts;
  int i;

  if( fgets(date_line, sizeof(date_line), fp)==NULL ){
    return(-1);
  }

  
    stat=split_dateline(date_line,mlDstr);
    fprintf(stderr,"date: %d %d %d\n  time: %d %d %d\n",mlDstr->yr,mlDstr->mo,mlDstr->dy,mlDstr->hr,mlDstr->mt,mlDstr->sc);

    fscanf(fp,"%f %f %f %f %f %f\n",&mlDstr->Bx,&mlDstr->By,&mlDstr->Bz,&mlDstr->Au,&mlDstr->Al,&mlDstr->v_sw);

    fscanf(fp,"%d",&npts);
    fprintf(stderr,"npts %d \n",npts);


    if( mlDstr->lats!=NULL )free(mlDstr->lats);
    mlDstr->lats=(double *)calloc(npts,sizeof(double));
    if( mlDstr->lons!=NULL) free(mlDstr->lons);
    mlDstr->lons=(double *)calloc(npts,sizeof(double));
    if( mlDstr->vn!=NULL) free(mlDstr->vn);
    mlDstr->vn=(double *)calloc(npts,sizeof(double));
    if( mlDstr->ve!=NULL) free(mlDstr->ve);
    mlDstr->ve=(double *)calloc(npts,sizeof(double));
    if( mlDstr->vn_cov!=NULL) free(mlDstr->vn_cov);
    mlDstr->vn_cov=(double *)calloc(npts,sizeof(double));
    if( mlDstr->ve_cov!=NULL) free(mlDstr->ve_cov);
    mlDstr->ve_cov=(double *)calloc(npts,sizeof(double));
    if( mlDstr->vmag!=NULL) free(mlDstr->vmag);
    mlDstr->vmag=(double *)calloc(npts,sizeof(double));
    if( mlDstr->vaz!=NULL) free(mlDstr->vaz);
    mlDstr->vaz=(double *)calloc(npts,sizeof(double));

   
    for( i=0; i<npts; i++ ){
      if((stat=fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf\n",
		      (mlDstr->lats+i),(mlDstr->lons+i),(mlDstr->vn+i),
		      (mlDstr->ve+i),&n_rmse,&e_rmse,
		      (mlDstr->vmag+i),(mlDstr->vaz+i)))==EOF)
	{return(-1);}
      mlDstr->vn_cov[i]=MIN(n_rmse*n_rmse,MAX_ML_COV);
      mlDstr->ve_cov[i]=MIN(e_rmse*e_rmse,MAX_ML_COV);


      /* mlDstr->vn_cov[i]=n_rmse; */
      /* mlDstr->ve_cov[i]=e_rmse; */
      mlDstr->lons[i]=(float)inv_MLTConvertYMDHMS_v2(mlDstr->yr,
						     mlDstr->mo,mlDstr->dy,
						     mlDstr->hr,mlDstr->mt,
						     mlDstr->sc,mlDstr->lons[i]);

      if( mlDstr->lons[i]<0 )mlDstr->lons[i]+=360;
    }
    mlDstr->npts=npts;
    return(ftell(fp));
}

int regrid_ml(MLDstr *mlDstr, GridMLDstr *rgMLDstr){

  double lat;
  double lon;
  int indx;
  int ll,lr,ul,ur;
  int j;
  int ml_lat_indx;
  
  double min_lon=1000;
  double max_lon=0;
  double a0, az0, b0, bz0;
  double a1, az1, b1, bz1;
  double a2, az2, b2, bz2;
  double a3, az3, b3, bz3;
  double area0, area1, area2, area3;
  double tot_area;
  double *vn,*ve,*vn_cov,*ve_cov;

  double dlon_l,dlat_l;
  
  vn=mlDstr->vn;
  ve=mlDstr->ve;
  vn_cov=mlDstr->vn_cov;
  ve_cov=mlDstr->vn_cov;

  for( j=0; j<2*ML_NLON; j++ ){
    if( mlDstr->lons[j]>max_lon )max_lon=mlDstr->lons[j];
    if( mlDstr->lons[j]<min_lon )min_lon=mlDstr->lons[j];
  }

  dlat_l=ML_DLAT;
  for( j=0; j<ngrid; j++ ){
    lat=grid[j].center_lat;
    lon=grid[j].center_lon;    
    
    if( fabs(lat) < ML_LAT_0 ){
      rgMLDstr->vn[j]=0;
      rgMLDstr->ve[j]=0;
      
      rgMLDstr->vn_cov[j]=MIN_ML_ERR;
      rgMLDstr->ve_cov[j]=MIN_ML_ERR;
      continue;
    }

    ml_lat_indx=(int)((fabs(lat)-ML_LAT_0)/ML_DLAT);
    indx=ml_lat_indx*ML_NLON;

    if( lon<mlDstr->lons[indx] )while( mlDstr->lons[indx]>mlDstr->lons[indx-1] )indx++;

    if( lon < min_lon ){
      indx=MAX(0,indx-ML_NLON);      
      while((mlDstr->lats[indx]<lat) && (mlDstr->lons[indx+1]>mlDstr->lons[indx]) && indx<mlDstr->npts){
	indx++;
      }
    }
    
    ll=indx;

    while(( mlDstr->lons[indx] < lon )&&( indx< mlDstr->npts )){
      ll=indx;
      indx++;
    }
    
    lr=ll+1;
    ul=ll+ML_NLON;
    ur=ul+1;

    if( lon>max_lon ){
      indx=MAX(0,ml_lat_indx*ML_NLON-1);
      while((mlDstr->lats[indx]<lat) && (mlDstr->lons[indx+1]>mlDstr->lons[indx]) && indx<mlDstr->npts){
	indx++;
      }
      ll=indx;
      ul=ll+ML_NLON;
      lr=indx+1;
      if(mlDstr->lats[lr]>lat)lr=MAX(0,lr-ML_NLON);
      ur=lr+ML_NLON;
    }
    
    if(mlDstr->lats[lr]>lat){
      lr=MAX(0,lr-ML_NLON);
      ur=lr+ML_NLON;
    }

    
    dlon_l=fabs(mlDstr->lons[lr]-mlDstr->lons[ll]);
    if( dlon_l>360 )dlon_l=fabs(dlon_l-360);
      
    if( ll+ML_NLON >= mlDstr->npts ){
      rgMLDstr->vn[j]=(vn[ll]*(mlDstr->lons[lr]-lon)+vn[lr]*(lon-mlDstr->lons[ll]))/dlon_l;
      rgMLDstr->ve[j]=(ve[ll]*(mlDstr->lons[lr]-lon)+ve[lr]*(lon-mlDstr->lons[ll]))/dlon_l;
      rgMLDstr->vn_cov[j]=(vn_cov[ll]*(mlDstr->lons[lr]-lon)+vn_cov[lr]*(lon-mlDstr->lons[ll]))/dlon_l;
      rgMLDstr->ve_cov[j]=(ve_cov[ll]*(mlDstr->lons[lr]-lon)+ve_cov[lr]*(lon-mlDstr->lons[ll]))/dlon_l;
      continue;
    }

    if( lat==mlDstr->lats[ll] ){
      rgMLDstr->vn[j]=(vn[ll]*(mlDstr->lons[lr]-lon)+vn[lr]*(lon-mlDstr->lons[ll]))/dlon_l;
      rgMLDstr->ve[j]=(ve[ll]*(mlDstr->lons[lr]-lon)+ve[lr]*(lon-mlDstr->lons[ll]))/dlon_l;
      rgMLDstr->vn_cov[j]=(vn_cov[ll]*(mlDstr->lons[lr]-lon)+vn_cov[lr]*(lon-mlDstr->lons[ll]))/dlon_l;
      rgMLDstr->ve_cov[j]=(ve_cov[ll]*(mlDstr->lons[lr]-lon)+ve_cov[lr]*(lon-mlDstr->lons[ll]))/dlon_l;
      continue;
    }

    if( lon==mlDstr->lons[lr] ){
      rgMLDstr->vn[j]=(vn[lr]*(mlDstr->lats[ur]-lat)+vn[ur]*(lat-mlDstr->lats[lr]))/dlat_l;
      rgMLDstr->ve[j]=(ve[lr]*(mlDstr->lats[ur]-lat)+ve[ur]*(lat-mlDstr->lats[lr]))/dlat_l;
      rgMLDstr->vn_cov[j]=(vn_cov[lr]*(mlDstr->lats[ur]-lat)+vn_cov[ur]*(lat-mlDstr->lats[lr]))/dlat_l;
      rgMLDstr->ve_cov[j]=(ve_cov[lr]*(mlDstr->lats[ur]-lat)+ve_cov[ur]*(lat-mlDstr->lats[lr]))/dlat_l;
      continue;
    }
    
    if( lon==mlDstr->lons[ll] ){
      rgMLDstr->vn[j]=(vn[ll]*(mlDstr->lats[ul]-lat)+vn[ul]*(lat-mlDstr->lats[ll]))/dlat_l;
      rgMLDstr->ve[j]=(ve[ll]*(mlDstr->lats[ul]-lat)+ve[ul]*(lat-mlDstr->lats[ll]))/dlat_l;
      rgMLDstr->vn_cov[j]=(vn_cov[ll]*(mlDstr->lats[ul]-lat)+vn_cov[ul]*(lat-mlDstr->lats[ll]))/dlat_l;
      rgMLDstr->ve_cov[j]=(ve_cov[ll]*(mlDstr->lats[ul]-lat)+ve_cov[ul]*(lat-mlDstr->lats[ll]))/dlat_l;
      continue;
    }

    if( (fabs(lon-mlDstr->lons[ll])<0.1)||(fabs(lat-mlDstr->lats[ll])<0.1) ){
      area0=0;
    }else{
      sub_sphazm(lon,lat,mlDstr->lons[ll],lat,&az0, &a0);
      sub_sphazm(lon,lat,lon,mlDstr->lats[ll],&bz0, &b0);
      area0=a0*b0;
    }
    
    if( (fabs(lon-mlDstr->lons[lr])<0.1)||(fabs(lat-mlDstr->lats[lr])<0.1) ){
      area1=0;
    }else{
      sub_sphazm(lon,lat,mlDstr->lons[lr],lat,&az1, &a1);
      sub_sphazm(lon,lat,lon,mlDstr->lats[lr],&bz1, &b1);
      area1=a1*b1;
    }
    
     if( (fabs(lon-mlDstr->lons[ul])<0.1)||(fabs(lat-mlDstr->lats[ul])<0.1) ){
      area2=0;
    }else{
       sub_sphazm(lon,lat,mlDstr->lons[ul],lat,&az2, &a2);
       sub_sphazm(lon,lat,lon,mlDstr->lats[ul],&bz2, &b2);
       area2=a2*b2;
     }

    if( (fabs(lon-mlDstr->lons[ur])<0.1)||(fabs(lat-mlDstr->lats[ur])<0.1) ){
      area3=0;
    }else{
      sub_sphazm(lon,lat,mlDstr->lons[ur],lat,&az3, &a3);
      sub_sphazm(lon,lat,lon,mlDstr->lats[ur],&bz3, &b3);
      area3=a3*b3;
    }

    tot_area=fabs(area0)+fabs(area1)+fabs(area2)+fabs(area3);
    
    rgMLDstr->vn[j]=(vn[ll]*area3+vn[lr]*area2+vn[ul]*area1+vn[ur]*area0)/tot_area;
    rgMLDstr->ve[j]=(ve[ll]*area3+ve[lr]*area2+ve[ul]*area1+ve[ur]*area0)/tot_area;
    

    rgMLDstr->vn_cov[j]=(vn_cov[ll]*area3+vn_cov[lr]*area2+vn_cov[ul]*area1+vn_cov[ur]*area0)/tot_area;
    rgMLDstr->ve_cov[j]=(ve_cov[ll]*area3+ve_cov[lr]*area2+ve_cov[ul]*area1+ve_cov[ur]*area0)/tot_area;

    if( !isfinite(rgMLDstr->vn_cov[j]) || !isfinite(rgMLDstr->ve_cov[j])||(fabs(rgMLDstr->vn[j])>1000)){
      fprintf(stderr,"not finite %d %d %d %d %d\n",j,ll,lr,ul,ur);
      fprintf(stderr,"ve %lf %lf %lf %lf %lf %lf\n",lat,lon,ve[ll],ve[lr],ve[ul],ve[ur]); 
      fprintf(stderr,"vn %lf %lf %lf %lf %lf %lf\n",lat,lon,vn[ll],vn[lr],vn[ul],vn[ur]); 
      fprintf(stderr,"ve_cov %lf %lf %lf %lf %lf %lf\n",lat,lon,ve_cov[ll],ve_cov[lr],ve_cov[ul],ve_cov[ur]); 
      fprintf(stderr,"vn_cov %lf %lf %lf %lf %lf %lf\n",lat,lon,vn_cov[ll],vn_cov[lr],vn_cov[ul],vn_cov[ur]); 
      fprintf(stderr,"a_s %lf %lf %lf %lf %lf %lf\n",lat,lon,a0,a1,a2,a3); 
      fprintf(stderr,"a_s %d %d %lf %lf %lf %lf\n",ll,lr,lat,lon,mlDstr->lons[ll],mlDstr->lons[ul]); 
      fprintf(stderr,"b_s %lf %lf %lf %lf %lf %lf\n",lat,lon,b0,b1,b2,b3);
      rgMLDstr->vn_cov[j]=(vn_cov[ll]+vn_cov[lr]+vn_cov[ul]+vn_cov[ur])/4;
      rgMLDstr->ve_cov[j]=(ve_cov[ll]+ve_cov[lr]+ve_cov[ul]+ve_cov[ur])/4;
    }

    if(!isfinite(rgMLDstr->ve[j]) || !isfinite(rgMLDstr->vn[j])){
      fprintf(stderr,"NAN: %d %f %f\n",j,rgMLDstr->ve[j],rgMLDstr->vn[j]);
      fprintf(stderr,"areas: %f %f %f %f\n",area0,area1,area2,area3);
      fprintf(stderr,"lats: %f %f %f %f %f\n",lat,mlDstr->lats[ll],mlDstr->lats[lr],mlDstr->lats[ul],mlDstr->lats[ur]);
      fprintf(stderr,"lons: %f %f %f %f %f\n",lon,mlDstr->lons[ll],mlDstr->lons[lr],mlDstr->lons[ul],mlDstr->lons[ur]);
      fprintf(stderr,"north values: %f %f %f %f\n",vn[ll],vn[lr],vn[ul],vn[ur]);
      fprintf(stderr,"east values: %f %f %f %f\n",ve[ll],ve[lr],ve[ul],ve[ur]);
      rgMLDstr->vn[j]=(vn[ll]+vn[lr]+vn[ul]+vn[ur])/4;
      rgMLDstr->ve[j]=(ve[ll]+ve[lr]+ve[ul]+ve[ur])/4;      
      rgMLDstr->vn_cov[j]=(vn_cov[ll]+vn_cov[lr]+vn_cov[ul]+vn_cov[ur])/4;
      rgMLDstr->ve_cov[j]=(ve_cov[ll]+ve_cov[lr]+ve_cov[ul]+ve_cov[ur])/4;
    }
  }
    
  return(0);
}


FILE_INFO *file=NULL;
FILE_INFO *select_file(char *radar, time_t time)
{
  char yr_str[5], mo_str[3], dy_str[3];
  char *raid_path;
  char dir_path[PATH_LEN];
  time_t ftime;
  time_t diftime;
  time_t mindif=100000;
  DIR *dp;
  struct dirent *ep;
  char *tz;

  tz = getenv("TZ");
  setenv("TZ", "", 1);
  tzset();

  raid_path=getenv("RAID_PATH");
  struct tm *in_time;
  in_time=gmtime(&time);
  int yr=in_time->tm_year+1900;
  int mo=in_time->tm_mon+1;
  int dy=in_time->tm_mday;
  /* int hr=in_time->tm_hour; */
  /* int mn=in_time->tm_min; */
  /* int sc=in_time->tm_sec; */

  CNV_TO_STR(yr,yr_str);
  CNV_TO_STR(mo,mo_str);
  CNV_TO_STR(dy,dy_str);
  
  sprintf(file->dir_path,"%s","");
  sprintf(file->fname,"%s","");

  sprintf(dir_path,"%s%s/%s.%s/",raid_path,yr_str,mo_str,dy_str);
  fprintf(stderr,"\n %s%s/%s.%s\n",raid_path,yr_str,mo_str,dy_str);
  fprintf(stderr,"\n dirpath: %s\n",dir_path);
  if((dp=opendir(dir_path))==NULL) {
    fprintf(stderr,"---- COULDN'T OPEN DATA DIRECTORY ----\n");
    closedir(dp);
    return(file);
  }

  while( (ep=readdir(dp)) )
    {
      if( strstr(ep->d_name,".gz") != NULL) continue;
      if( strstr(ep->d_name,".bz2") != NULL) continue;
      if( strstr(ep->d_name,radar) != NULL)
	{
	  ftime=fname_to_time(ep->d_name);
	  diftime=time-ftime;
	  if( diftime>=0 && diftime<mindif )
	  {
	    sprintf(file->dir_path,"%s",dir_path);
	    sprintf(file->fname,"%s",ep->d_name);
	    mindif=diftime;
	  }
	}
    }
  closedir(dp);
  return file;
}



time_t ftime(struct RadarParm *prm){
  int yr,mo,dy,hr,mt,sc;
  yr=prm->time.yr;
  mo=prm->time.mo;
  dy=prm->time.dy;
  hr=prm->time.hr;
  mt=prm->time.mt;
  sc=prm->time.sc;
  return (time_t)TimeYMDHMSToEpoch(yr,mo,dy,hr,mt,sc);
}

void get_pos_ar(struct RadarSite *site, struct RadarParm *prm, struct RadarPos *rdrpos){
  
  double lat,lon,alt=300;
  double magazm;
  int rn,bm,rsep,frang;
  int s;
  int rxrise,yr;
  int chisham=1;
  rsep=prm->rsep;
  frang=prm->frang;
  rxrise=prm->rxrise;
  yr=prm->time.yr;
  
  for( bm=0; bm<site->maxbeam; bm++) for( rn=0; rn<=prm->nrang; rn++){
      
      /* Calculate magnetic latitude, longitude, and azimuth of range/beam
       * position */
      s=RPosInvMag(bm,rn,yr,site,frang,rsep,rxrise,
		   alt,&lat,&lon,&magazm,chisham,0);
      
      
      rdrpos->lat[bm][rn]=lat;
      rdrpos->lon[bm][rn]=lon;      
      rdrpos->kazm[bm][rn]=magazm;
    }
}


void map_pos_to_grid(struct RadarPos rdrpos,int frang,int rsep,double bmsep, struct RadarMap *rmap){
  int ib,ir,ig,ig_last;
  int bsteps,ibstep;
  int count;
  double mlat,mlon,kazm,mlat_c,mlon_c;
  double bm_delta;

  bsteps=MAX_COUNT/2;

  for( ib=0; ib<MAX_BEAM; ib++)for( ir=0; ir<MAX_GATES; ir++){
      mlat_c=rdrpos.lat[ib][ir];
      mlon_c=rdrpos.lon[ib][ir];
      kazm=rdrpos.kazm[ib][ir];
      
      count=0;
      ig_last=-1;
      bm_delta=((double)frang+(double)rsep*(double)ir)*bmsep/(double)bsteps;
      /* fprintf(stderr,"\n map_pos_pt_grid: %d %d %d %f %f\n",ir,frang,rsep,bmsep,bm_delta); */
      for(ibstep=0; ibstep<bsteps; ibstep++){
	sub_sphcal(mlon_c, mlat_c, kazm-90, (double)ibstep*bm_delta, &mlon, &mlat);
	/* if(mlat_c !=0)fprintf(stderr,"map_pos_to_grid: %d %f %f %f %f %f\n",ir,mlon_c,mlat_c,kazm,mlon,mlat); */
	for(ig=0; ig<ngrid; ig++)if(cell_in_out(mlat,mlon,grid[ig])){
	    if( ig != ig_last){
	      rmap->cell[ib][ir][count]=ig;
	      ig_last=ig;
	      count++;
	    }
	  }
      }
      
      for(ibstep=0; ibstep<bsteps; ibstep++){
	sub_sphcal(mlon_c, mlat_c, kazm+90, (double)ibstep*bm_delta, &mlon, &mlat);
	for(ig=0; ig<ngrid; ig++)if(cell_in_out(mlat,mlon,grid[ig])){
	    if( ig != ig_last){
	      rmap->cell[ib][ir][count]=ig;
	      ig_last=ig;
	      count++;
	    }
	  }
      }
      rmap->cell_count[ib][ir]=count;      
    }
  /* sleep(1); */
}
  
double grid_lat(int ig){
  return((grid[ig].lat[0]+grid[ig].lat[1]+grid[ig].lat[2]+grid[ig].lat[3])/4);
}

double grid_lon(int ig){
  return((grid[ig].lon[0]+grid[ig].lon[1]+grid[ig].lon[2]+grid[ig].lon[3])/4);
}



void corotation_radar(struct RadarSite *site, double* cor_mag, double* cor_azm){
  double lat0,lon0,lat1,lon1,h;
  double azm,rng;
  int err;
  
  err = AACGM_v2_Convert(site->geolat,site->geolon,300, &lat0,&lon0, &h, 0);    
  err = AACGM_v2_Convert(site->geolat,site->geolon+1, 300, &lat1,&lon1, &h, 0);
  sub_sphazm(lon0,lat0,lon1,lat1,&azm, &rng);
  
  *cor_mag=2*PI*LRE*1000*cosd(site->geolat)/SECS_P_DAY;
  *cor_azm=azm;
}


struct RadarParm *prm;
struct FitData *fit;
struct FitIndex *inx;

struct RadarNetwork *network;  
struct Radar *radar;

struct RadarPos *rpos;


extern int cgls(int m,int n,double* a[],double* b, double* x);
MLDstr mlDstr;
GridMLDstr rgMLDstr;
	       
int main(int argc, char *argv[]){

  /* read instruction file
     file should have date, time, lat and lon of area corners, radars to contribute
     filtering instructions
  */

  FILE *fp;
  FILE *vf;
  FILE *mlf;
  FILE *rptrs[NRAD];
  FILE *dataout;
  FILE *losout;
  FILE *divfile;
  char file_name[128];
  int nbms;
  struct tm *t_start;
  struct tm *t_end; 
  time_t fit_time, esec;
  char *envstr;
  struct RadarSite *site=NULL;
  struct RadarMap *rmap;
  int *grid_data_count;
  gsl_matrix *coef;
  gsl_matrix *coefTrans;
  gsl_vector *data;
  gsl_matrix *A;
  gsl_matrix *CDI;  
  gsl_vector *GD;  
  int stat;
  double *solution_l;
  double ml_time;

  double **aa;
  double *bb;

  double min_lat=MIN_LAT;
  double max_lat=MAX_LAT;
  double min_lon=0;
  double max_lon=360;


  int fpos;
  int yr,mo,dy,hr,mt;
  /* int syr,smo,sdy,shr,smt; */
  /* int eyr,emo,edy,ehr,emt; */
  /* double sc,ssc,esc; */
  int j,jr,jc,jeqn;
  int neqn,ndata;
  clock_t t1,t2;
  double sc;
  
  int yr_st,mo_st,dy_st,hr_st,mt_st;
  int fit_yr,fit_mo,fit_dy,fit_hr,fit_mt;
  double sc_st,fit_sc;
  int jj,jjc,jjr,count;
  int jcmn,jcmx,jrmn,jrmx;
  int frang,rsep;
  int *rsepl;
  int *frangl;
  int *nrangl;
  double bmsep;
  double kazm;
  double avg_val;
  double cor_mag,cor_azm;
  int grid_cell;
  int try_count;
  int count_lim=20;
  int s,jrmin,jrmax;

  double **filt_ar;
  double **filt_ar1;
  double **filt_err;
  double **filt_cnt;
  double med_filt_ar[9],median,tmp,tp;
  int **good_ar;

  mlDstr.lats=NULL;
  mlDstr.lons=NULL;
  mlDstr.ve=NULL;
  mlDstr.vn=NULL;
  mlDstr.ve_cov=NULL;
  mlDstr.vn_cov=NULL;
  mlDstr.vmag=NULL;
  mlDstr.vaz=NULL;
  
  smooth=0;
  MODEL_SCALE=1.;

  if (argc <= 1){
    fprintf(stderr,"********NO INSTRUCTION FILE GIVEN********\n");
    exit(-1);
  }
  strcpy(file_name,argv[1]);
  if ((fp=fopen(file_name,"r")) == 0) {
    fprintf(stderr,"******FIT INSTRUCTION FILE %s NOT FOUND******\n",argv[1]);
    exit(-1);
  }
  /* parse file */
  parse_instructions(fp);
  fclose(fp);


  t_start=parse_date_str(start_time);
  t_end=parse_date_str(end_time);

  fit_time=mktime(t_start);
  esec=mktime(t_end);  

  
  TimeEpochToYMDHMS(fit_time,&yr,&mo,&dy,&hr,&mt,&sc);
  sc=0;
  fprintf(stderr,"%d %d %d %d %d %f\n",yr,mo,dy,hr,mt,sc);
  IGRF_SetDateTime(yr,mo,dy,hr,mt,(int)sc);
  AACGM_v2_SetDateTime(yr,mo,dy,hr,mt,(int)sc);  
  
  min_lon_ng=99999.9;
  max_lon_ng=-99999.9;
  for( j=0; j<nbp; j++ ){
    if( bp[j].lon<min_lon_ng )min_lon_ng=bp[j].lon;
    if( bp[j].lon>max_lon_ng )max_lon_ng=bp[j].lon;
  }
  
  /* create grid
   */
  
  int nlat=(max_lat-min_lat)/dlat;
  int jlat,nlon;
  double lat,dlon_l;
  start_lon=(double *)calloc(nlat,sizeof(double));

  srand(time(NULL)); // randomize seed
  for( jlat=0; jlat<nlat; jlat++ )
    {
      lat=min_lat+((double)jlat+.5)*dlat;      
      dlon_l=dlat/cosd(lat);
      nlon=(int)((max_lon-min_lon)/dlon_l);
      dlon_l=360/(double)nlon;
      start_lon[jlat]=get_random()*dlon_l/2;
    }
    
  get_grid_size(min_lat, max_lat);
  grid=(CELL *)calloc(ngrid,sizeof(CELL));
  
  grd_file=fopen("grid.dat","w");
  make_grid(min_lat, max_lat);


  if( strcmp(hemisphere,"south")==0 ){
    max_lat=-MIN_LAT;
    min_lat=-MAX_LAT;
    
    /* dlat*=-1; */
    /* dlat_ng*=-1; */
  }

  
  fprintf(grd_file,"min_lat: %f\n",min_lat);
  fprintf(grd_file,"max_lat: %f\n",max_lat);
  fprintf(grd_file,"min_lon: %f\n",min_lon);
  fprintf(grd_file,"max_lon: %f\n",max_lon);
  fprintf(grd_file,"ngrid: %d\n",ngrid);
  
  for( j=0; j<ngrid; j++ ){
    fprintf(grd_file,"%d %f  %f\n",j,grid[j].center_lat,grid[j].center_lon);
    fprintf(grd_file,"%f  %f  %f  %f\n",grid[j].lat[0],grid[j].lat[1],grid[j].lat[2],grid[j].lat[3]);
    fprintf(grd_file,"%f  %f  %f  %f\n",grid[j].lon[0],grid[j].lon[1],grid[j].lon[2],grid[j].lon[3]);    
  }
  fclose(grd_file);

  grid_data_count=malloc(ngrid*sizeof(int));
  
  rgMLDstr.vn=(double *)calloc(ngrid,sizeof(double));
  rgMLDstr.ve=(double *)calloc(ngrid,sizeof(double));
  rgMLDstr.vn_cov=(double *)calloc(ngrid,sizeof(double));
  rgMLDstr.ve_cov=(double *)calloc(ngrid,sizeof(double));
  
  neighbors=calloc(ngrid,sizeof(struct Neighbor));
  find_neighbors();
  
  /* create array of coefficients for calculating divergence */
  divfile=fopen("divergence.dat","w");
  make_div_ar();

  for( jr=0; jr<ngrid-nedge; jr++){
    fprintf(divfile,"%d ",jr);
    for( jc=0; jc<2*ngrid; jc++)if(div_array[jr].coef[jc] !=0 )fprintf(divfile,"   %d %f\n",jc,div_array[jr].coef[jc]);
    fprintf(divfile,"\n");
  }
  fclose(divfile);
  
  envstr=getenv("SD_RADAR");
  if (envstr==NULL) {
    fprintf(stderr,"Environment variable 'SD_RADAR' must be defined.\n");
    exit(-1);
  }

  fp=fopen(envstr,"r");

  if (fp==NULL) {
    fprintf(stderr,"Could not locate radar information file.\n");
    exit(-1);
  }

  network=RadarLoad(fp);
  fclose(fp); 
  if (network==NULL) {
    fprintf(stderr,"Failed to read radar information.\n");
    exit(-1);
  }

  envstr=getenv("SD_HDWPATH");
  if (envstr==NULL) {
    fprintf(stderr,"Environment variable 'SD_HDWPATH' must be defined.\n");
    exit(-1);
  }

  RadarLoadHardware(envstr,network);

  /* find fit files and read into model-array and data vector
     requires determining k-vectors
   */

  if(file)free(file);
  file=calloc(1,sizeof(FILE_INFO));

  int res;
  nrad--;
  for( jj=0; jj<nrad; jj++ )
    {
      rptrs[jj]=NULL;
      file=select_file(radar_list[jj],fit_time);
      /* if ( file->fname !=NULL){ */
      if ( (res=strcmp(file->fname,"")) !=0){
	fprintf(stderr,"%s%s\n",file->dir_path,file->fname);
	sprintf(file_name,"%s/%s",file->dir_path,file->fname);
	if(rptrs[jj])fclose(rptrs[jj]);
	if((rptrs[jj]=fopen(file_name,"r"))== NULL)fprintf(stderr,"could not open file: %s",file_name);
      }
    }
  fprintf(stderr,"first files open\n");
  
  if ((mlf=fopen(ml_file,"r")) == 0) {
    fprintf(stderr,"******ML FILE %s NOT FOUND******\n",ml_file);
    exit(-1);
  }

  fpos=read_ml_record(mlf, &mlDstr);  
  
  
  fprintf(stderr,"read ml_record\n");
  
  vf=fopen("vel_out","w");
  losout=fopen("los_out","w");
  fp=fopen("coef_array","w");

  nrangl=calloc(nrad, sizeof(int));
  rsepl=calloc(nrad, sizeof(int));
  frangl=calloc(nrad, sizeof(int));
  
  prm=RadarParmMake();
  fit=FitMake();
  rpos=calloc(nrad, sizeof(struct RadarPos));
  rmap=calloc(nrad, sizeof(struct RadarMap));
  
  
  while( fit_time<esec )
    {
      TimeEpochToYMDHMS(fit_time,&fit_yr,&fit_mo,&fit_dy,&fit_hr,&fit_mt,&fit_sc);      
      TimeEpochToYMDHMS(fit_time-(time_t)avg_ival/2,&yr_st,&mo_st,&dy_st,&hr_st,&mt_st,&sc_st);
      printf("%d %d %d %d %d %f\n",yr_st,mo_st,dy_st,hr_st,mt_st,sc_st);
      fprintf(stderr,"fit_time time: %d %d %d %d %d %f\n",fit_yr,fit_mo,fit_dy,fit_hr,fit_mt,fit_sc);
      fprintf(stderr,"Interval start time: %d %d %d %d %d %f\n",yr_st,mo_st,dy_st,hr_st,mt_st,sc_st);
      ndata=0;
      for(jj=0; jj<ngrid; jj++)grid_data_count[jj]=0;
      for( jj=0; jj<nrad; jj++){


	if( rptrs[jj]!=NULL ){
	  s=fseek(rptrs[jj],0L,SEEK_SET);
	  if( s==0 )s=FitFseek(rptrs[jj],yr_st,mo_st,dy_st,hr_st,mt_st,sc_st,NULL,inx);
	  if( s==0 )s=FitFread(rptrs[jj],prm,fit);
	  if( s==0 )fprintf(stderr,"%s file time: %d %d %d %d %d %d %d %ld\n",radar_list[jj],prm->time.yr,
			    prm->time.mo,prm->time.dy,prm->time.hr,prm->time.mt,prm->time.sc,prm->stid,
			    ftime(prm)-fit_time);
	  if( s!=0 ){
	    fprintf(stderr,"1. File does not contain the requested interval. %d:%d\n",hr_st,mt_st);
	    if(rptrs[jj]){
	      fprintf(stderr,"*******Closing %s file\n",radar_list[jj]);
	      fclose(rptrs[jj]);
	      rptrs[jj]=NULL;
	    }
	  }
	}

	if( rptrs[jj]==NULL ){
	  file=select_file(radar_list[jj],fit_time-(time_t)avg_ival/2);
	  if ( (res=strcmp(file->fname,"")) ==0)continue;
	  /* if( file->fname==NULL ) continue; */
	  sprintf(file_name,"%s/%s",file->dir_path,file->fname);
	  fprintf(stderr,"*******Opening file %s\n",file_name); 
	  rptrs[jj]=fopen(file_name,"r");
	  if( rptrs[jj]!=NULL ){
	    s=FitFseek(rptrs[jj],yr_st,mo_st,dy_st,hr_st,mt_st,sc_st,NULL,inx);
	    if( s==0 )s=FitFread(rptrs[jj],prm,fit);
	    /* TimeEpochToYMDHMS(fit_time+(time_t)avg_ival,&yr,&mo,&dy,&hr,&mt,&sc); */
	  }
	}
	  
	radar=RadarGetRadar(network,prm->stid);
	site=RadarYMDHMSGetSite(radar,yr_st,mo_st,dy_st,hr_st,mt_st,(int) 0);

	if( INERTIAL ){
	  corotation_radar(site,&cor_mag,&cor_azm);
	}
	
	if (site==NULL) {fprintf(stderr,"NULL site\n"); continue;}

	nbms=(int)site->maxbeam;
	filt_ar=malloc(nbms*sizeof(double*));
	for( jr=0; jr<nbms; jr++)filt_ar[jr]=(double*)calloc(prm->nrang,sizeof(double));
	filt_ar1=malloc(nbms*sizeof(double*));
	for( jr=0; jr<nbms; jr++)filt_ar1[jr]=(double*)calloc(prm->nrang,sizeof(double));
	filt_err=malloc(nbms*sizeof(double*));
	for( jr=0; jr<nbms; jr++)filt_err[jr]=(double*)calloc(prm->nrang,sizeof(double));
	filt_cnt=malloc(nbms*sizeof(double*));
	for( jr=0; jr<nbms; jr++)filt_cnt[jr]=(double*)calloc(prm->nrang,sizeof(double));
	good_ar=malloc(nbms*sizeof(int*));
	for( jr=0; jr<nbms; jr++)good_ar[jr]=(int*)calloc(prm->nrang,sizeof(int));
	nrangl[jj]=prm->nrang;
	
	try_count=0;
	if( rptrs[jj]!=NULL ){
	  s=FitFseek(rptrs[jj],yr_st,mo_st,dy_st,hr_st,mt_st,sc_st,NULL,inx);
	  while( ((ftime(prm)-fit_time)<avg_ival/2) && (try_count<count_lim)){
	    if( s!=0 ){
	      try_count++;
	      fprintf(stderr,"2. File does not contain the requested interval. %d:%d\n",hr_st,mt_st);
	      fprintf(stderr,"%s  %ld\n",radar_list[jj],fit_time+(time_t)avg_ival);
	      file=select_file(radar_list[jj],fit_time+(time_t)avg_ival);
	      sprintf(file_name,"%s/%s",file->dir_path,file->fname);
	      if(rptrs[jj]){
		fprintf(stderr,"*******Closing %s file\n",radar_list[jj]);
		fclose(rptrs[jj]);
	      }
	      rptrs[jj]=NULL;
	      /* if( file->fname == NULL )continue; */
	      if ( (res=strcmp(file->fname,"")) ==0)continue;
	      fprintf(stderr,"*******Opening file %s\n",file_name); 
	      if( (rptrs[jj]=fopen(file_name,"r")) ==NULL )continue;
	      /*	      rewind(rptrs[jj]);*/
	      /* TimeEpochToYMDHMS(fit_time+(time_t)avg_ival,&yr,&mo,&dy,&hr,&mt,&sc); */
	      if((s=FitFseek(rptrs[jj],fit_yr,fit_mo,fit_dy,fit_hr,fit_mt,fit_sc,NULL,inx))!=0 ||
		 (s=FitFread(rptrs[jj],prm,fit))!=0){
		fprintf(stderr,"problem reading %s\n",file_name);
		fclose(rptrs[jj]);
		rptrs[jj]=NULL;
		break;
	      }
	      fprintf(stderr,"main: %s\n",file->fname);
	    }


	    if( s==0 ){

	      rsep=prm->rsep;
	      frang=prm->frang;
	      if((rsep != rsepl[jj]) || (frang != frangl[jj])){
		fprintf(stderr,"recalculating position array %s %d %d  %d  %d\n",
			radar_list[jj],rsep,rsepl[jj],frang,frangl[jj]);
		radar=RadarGetRadar(network,prm->stid);
		site=RadarYMDHMSGetSite(radar,fit_yr,fit_mo,fit_dy,fit_hr,fit_mt,(int) fit_sc);
		if (site==NULL){fprintf(stderr,"site was null\n"); continue;}
		bmsep=(double)site->bmsep*(double)PI/180.0;
		get_pos_ar(site,prm,&rpos[jj]);
		map_pos_to_grid(rpos[jj],frang,rsep,bmsep,&rmap[jj]);
		rsepl[jj]=rsep;
		frangl[jj]=frang;
	      }
	      if(prm->nrang != nrangl[jj]){
		fprintf(stderr,"%d station: %s  nrang mismatch old: %d   new: %d...reallocating\n",jj,radar_list[jj],nrangl[jj],prm->nrang);

		for( jr=0; jr<nbms; jr++)if(filt_ar[jr] != NULL)free(filt_ar[jr]); 
		for( jr=0; jr<nbms; jr++)if(filt_ar1[jr] != NULL)free(filt_ar1[jr]); 
		for( jr=0; jr<nbms; jr++)if(filt_err[jr] != NULL)free(filt_err[jr]); 
		for( jr=0; jr<nbms; jr++)if(filt_cnt[jr] != NULL)free(filt_cnt[jr]); 
		for( jr=0; jr<nbms; jr++)if(good_ar[jr] != NULL)free(good_ar[jr]);


		for( jr=0; jr<nbms; jr++)filt_ar[jr]=(double*)calloc(prm->nrang,sizeof(double));
		for( jr=0; jr<nbms; jr++)filt_ar1[jr]=(double*)calloc(prm->nrang,sizeof(double));
		for( jr=0; jr<nbms; jr++)filt_err[jr]=(double*)calloc(prm->nrang,sizeof(double));
		for( jr=0; jr<nbms; jr++)filt_cnt[jr]=(double*)calloc(prm->nrang,sizeof(double));
		for( jr=0; jr<nbms; jr++)good_ar[jr]=(int*)calloc(prm->nrang,sizeof(int));
		
		if(prm->nrang == 0)for( jr=0; jr<nbms; jr++){
		    filt_ar[jr]=NULL;
		    filt_ar1[jr]=NULL;
		    filt_err[jr]=NULL;
		    filt_cnt[jr]=NULL;
		    good_ar[jr]=NULL;
		  }		  
		nrangl[jj]=prm->nrang;
	      }
	      jrmin=(MIN_RANGE-frang)/rsep;
	      jrmax=(MAX_RANGE-frang)/rsep;
	      if( labs(ftime(prm)-fit_time) < (double)avg_ival/2){
		yr=prm->time.yr;
		mo=prm->time.mo;
		dy=prm->time.dy;
		hr=prm->time.hr;
		mt=prm->time.mt;
		sc=(double)prm->time.sc;
		/* fprintf(stderr,"%d %d %d %d %d %d\n",yr,mo,dy,hr,mt,(int)sc); */
		/* fprintf(stderr,"fit %d %d %d %d %d %d\n",fit_yr,fit_mo,fit_dy,fit_hr,fit_mt,(int)fit_sc); */
		for( jr=jrmin; jr<jrmax; jr++ )
		  if(fit->rng[jr].qflg == 1 && fit->rng[jr].gsct == 0 && 
		     (fit->rng[jr-1].gsct == 0 && fit->rng[jr+1].gsct == 0) &&
		     fabs(fit->rng[jr].v_err) < MAX_V_ERR &&
		     fabs(fit->rng[jr].v) < MAX_V && fabs(fit->rng[jr].v) > MIN_V){
		    
		    filt_ar[prm->bmnum][jr]+=fit->rng[jr].v;
		    filt_err[prm->bmnum][jr]+=ERR_SCALE*fit->rng[jr].v_err*fit->rng[jr].v_err;
		    filt_cnt[prm->bmnum][jr]+=1.;
		    
		  }
	      }
	    }
	    s=FitFread(rptrs[jj],prm,fit);
	  }
	}

	for( jc=0; jc<site->maxbeam; jc++ )for( jr=1; jr<prm->nrang; jr++){
	    jcmn=MAX(jc-1,0);
	    jcmx=MIN(jc+1,site->maxbeam-1);
	    jrmn=MAX(jr-1,0);
	    jrmx=MIN(jr+1,prm->nrang-1);
	    count=0;
	    avg_val=0;
	    filt_ar1[jc][jr]=0;
	    good_ar[jc][jr]=0;
	    if( filt_cnt[jc][jr]==0 ) continue;
	    
	    for(jjc=jcmn; jjc<=jcmx; jjc++)for(jjr=jrmn; jjr<=jrmx; jjr++){
		if( filt_cnt[jjc][jjr]!=0 ){
		  med_filt_ar[count]=filt_ar[jjc][jjr]/filt_cnt[jjc][jjr];
		  count++;
		  avg_val+=(filt_ar[jjc][jjr]/filt_cnt[jjc][jjr]);
		}
	      }

	    if( count<MIN_COUNT)continue;
	    
	    for( jjc=0; jjc<count; jjc++ )for( jjr=jjc+1; jjr<count; jjr++ ){
		tp=med_filt_ar[jjc];
		if( tp>med_filt_ar[jjr] ){
		  tmp=med_filt_ar[jjr];
		  med_filt_ar[jjr]=tp;
		  med_filt_ar[jjc]=tmp;
		  tp=tmp;
		}
	      }
	    
	    if( count%2 ){
	      median=med_filt_ar[(int)(count/2)];
	    }else{
	      median=(med_filt_ar[(int)(count/2)]+med_filt_ar[(int)(count/2)+1])/2;
	    }
		
	    /* for(jjc=0; jjc<count; jjc++) */
	    /*   fprintf(stderr,"%d med_filt_ar %lf \n",jjc,med_filt_ar[jjc]); */
	    
	    /* fprintf(stderr,"The median was: %lf\n",median); */

	    
	    avg_val/=(double)count;
	    
	    good_ar[jc][jr]=1;

	    /* filt_ar1[jc][jr]=avg_val; */
	    filt_ar1[jc][jr]=median;
	    if(count<MIN_COUNT){
	      good_ar[jc][jr]=0;
	      filt_ar1[jc][jr]=0;
	      continue;
	    }
	  }
	for( jc=0; jc<nbms; jc++ )for( jr=1; jr<prm->nrang; jr++){
	    for( jjc=0; jjc<rmap[jj].cell_count[jc][jr]; jjc++){
	      if ((grid_cell=rmap[jj].cell[jc][jr][jjc]) != -1 && good_ar[jc][jr] == 1 &&
		  filt_cnt[jc][jr]>0 && fabs(rpos[jj].lat[jc][jr]) < MAX_DATA_LAT){	
		grid_data_count[grid_cell]++;
		    
		kazm=rpos[jj].kazm[jc][jr];
	      
		/* printf("%d %d %d %f %f %f %f %f %f %f\n",prm->stid,jc,jr,rpos[jj].lat[jc][jr], */
		/* 	     rpos[jj].lon[jc][jr],kazm,filt_ar1[jc][jr]/filt_cnt[jc][jr], */
		/* 	     grid_lat(grid_cell),grid_lon(grid_cell)); */
		kazm_array=realloc(kazm_array,(ndata+1)*sizeof(struct mod_array));
		kazm_array[ndata].coef=calloc(2*ngrid,sizeof(double));
		los_data=realloc(los_data,(ndata+1)*sizeof(double));
		in_los_data=realloc(in_los_data,(ndata+1)*sizeof(double));
		los_kazm=realloc(los_kazm,(ndata+1)*sizeof(double));
		los_lats=realloc(los_lats,(ndata+1)*sizeof(double));
		los_lons=realloc(los_lons,(ndata+1)*sizeof(double));
		los_err=realloc(los_err,(ndata+1)*sizeof(double));
		kazm_array[ndata].coef[grid_cell]=sind(kazm);
		kazm_array[ndata].coef[grid_cell+ngrid]=cosd(kazm);
		
		los_lats[ndata]=grid_lat(grid_cell);
		los_lons[ndata]=grid_lon(grid_cell);
		los_kazm[ndata]=kazm;
		in_los_data[ndata]=filt_ar1[jc][jr];
		los_data[ndata]=(-1.)*filt_ar1[jc][jr]; /* +cor_mag*cosd(kazm-cor_azm); */
		los_err[ndata]=filt_err[jc][jr]/filt_cnt[jc][jr];
		/* fprintf(stderr,"%d %d %f %f %f %f %f\n",ndata,prm->stid,los_lats[ndata],los_lons[ndata],los_data[ndata],los_err[ndata],cor_mag*cosd(kazm-cor_azm)); */
		
		ndata++;
	      }
	    }
	  }
	for( jr=0; jr<nbms; jr++)if(filt_ar[jr] != NULL)free(filt_ar[jr]); 
	free(filt_ar);
	for( jr=0; jr<nbms; jr++)if(filt_ar1[jr] != NULL)free(filt_ar1[jr]); 
	free(filt_ar1);
	for( jr=0; jr<nbms; jr++)if(filt_err[jr] != NULL)free(filt_err[jr]); 
	free(filt_err);
	for( jr=0; jr<nbms; jr++)if(filt_cnt[jr] != NULL)free(filt_cnt[jr]); 
	free(filt_cnt);
	for( jr=0; jr<nbms; jr++)if(good_ar[jr] != NULL)free(good_ar[jr]); 
	free(good_ar);
      }
    
    /* synchronize ML model with grid data */
    
    ml_time=TimeYMDHMSToEpoch(mlDstr.yr,mlDstr.mo,mlDstr.dy,mlDstr.hr,
			      mlDstr.mt,(double)mlDstr.sc);
    while( ml_time<fit_time ){
      if((fpos=read_ml_record(mlf, &mlDstr))==-1)break;
      ml_time=TimeYMDHMSToEpoch(mlDstr.yr,mlDstr.mo,mlDstr.dy,
				mlDstr.hr,mlDstr.mt,(double)mlDstr.sc);
      
    }

    fprintf(stderr,"grid: %d %d %d %d %d %lf\n",yr,mo,dy,hr,mt,sc);
    fprintf(stderr,"ML: %d %d %d %d %d %d\n",mlDstr.yr,mlDstr.mo,mlDstr.dy,mlDstr.hr,mlDstr.mt,mlDstr.sc);
    
    /* Interpolate ML model to grid */
    if( regrid_ml(&mlDstr,&rgMLDstr)!=-1){


      fprintf(stdout,"%d %d %d %d %d %d\n",yr,mo,dy,hr,mt,(int)sc);
      fprintf(stdout,"%d\n",ngrid);
      for( j=0; j<ngrid; j++ ){
	
	fprintf(stdout,"%5.2f %5.2f %8.3f %8.3f %8.3f %8.3f\n",grid[j].center_lat,grid[j].center_lon,rgMLDstr.ve[j],rgMLDstr.vn[j],rgMLDstr.ve_cov[j],rgMLDstr.vn_cov[j]);
      }
      
      neqn=3*ngrid-nedge+ndata;
      
      fprintf(stderr,"NEQN=%d\n",neqn);
      coef=gsl_matrix_calloc((size_t)neqn,(size_t)ngrid*2);
      coefTrans=gsl_matrix_calloc((size_t)ngrid*2,(size_t)neqn);
      CDI=gsl_matrix_calloc((size_t)neqn,(size_t)neqn);
      data=gsl_vector_calloc((size_t)neqn);
      
      
      fprintf(stderr,"***ALLOCATED***");
      
      
      jeqn=0;
      /* Edges are low-latitude boundary points, which are the first nedge points of the grid */
      for( jr=0; jr<ngrid; jr++){
	if(isnan(rgMLDstr.ve[jr])||isnan(rgMLDstr.vn[jr])){
	  fprintf(stderr,"NAN: %d %f %f\n",jr,rgMLDstr.ve[jr],rgMLDstr.vn[jr]);
	}
	gsl_matrix_set(coef,jeqn,jr,1);
	gsl_vector_set(data,jeqn,rgMLDstr.ve[jr]);	
	gsl_matrix_set(CDI,jeqn,jeqn,1/MAX(rgMLDstr.ve_cov[jr],MIN_ML_ERR));
	jeqn++;
	gsl_matrix_set(coef,jeqn,jr+ngrid,1);
	gsl_vector_set(data,jeqn,rgMLDstr.vn[jr]);
	gsl_matrix_set(CDI,jeqn,jeqn,1/MAX(rgMLDstr.vn_cov[jr],MIN_ML_ERR));
	jeqn++;
      }
      fprintf(stderr,"MODEL SET  %d***",jeqn);
      
      
      for( jr=0; jr<ngrid-nedge; jr++){
	for( jc=0; jc<2*ngrid; jc++)gsl_matrix_set(coef,jeqn,jc,div_array[jr].coef[jc]);
	gsl_vector_set(data,jeqn,0);
	gsl_matrix_set(CDI,jeqn,jeqn,1/DIV_ERR);
	jeqn++;
      }
      fprintf(stderr,"DIVERGENCE SET  %d***",jeqn);
      
      /* for( j=0; j<ngrid; j++ )grid_data_count[j]=0; */
      for( jr=0; jr<ndata; jr++){
	for( jc=0; jc<2*ngrid; jc++)gsl_matrix_set(coef,jeqn,jc,kazm_array[jr].coef[jc]);
	gsl_vector_set(data,jeqn,los_data[jr]);
	/* gsl_matrix_set(CDI,jeqn,jeqn,1/(MAX(los_err[jr]*los_err[jr],MIN_ERR))); */
	gsl_matrix_set(CDI,jeqn,jeqn,1/(MAX(los_err[jr],MIN_ERR)));
	jeqn++;
      }
      fprintf(stderr,"DATA SET***jeqn %d****\n",jeqn);
      
      
      
      gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,coef,CDI,0.,coefTrans); 
    
      GD=gsl_vector_calloc(2*ngrid);
      gsl_blas_dgemv(CblasNoTrans,1.,coefTrans,data,0.,GD);

      A=gsl_matrix_calloc(2*ngrid,2*ngrid);
      gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,coefTrans,coef,0.,A); 
    
      aa=malloc(2*ngrid*sizeof(double*));
      for( jr=0; jr<2*ngrid; jr++)aa[jr]=(double*)malloc(2*ngrid*sizeof(double));
      
      for( jr=0; jr<2*ngrid; jr++)for( jc=0; jc<2*ngrid; jc++)aa[jr][jc]=gsl_matrix_get(A,jr,jc);

    

      if( WRITE_COEF ){
    
	dataout=fopen("data_out","w");
	for( jr=0; jr<2*ngrid; jr++ ){
	  fprintf(dataout,"%lf\n",gsl_vector_get(GD,jr));
	}
	for( jr=0; jr<2*ngrid; jr++ )for(jc=0; jc<ngrid; jc++ )if( aa[jr][jc]!=0 )fprintf(dataout,"coefficient %d %d %lf\n",jr,jc,aa[jr][jc]); 
	fclose(dataout);
      }
      
      for( jr=0; jr<2*ngrid; jr++)for( jc=0; jc<2*ngrid; jc++)if( !isfinite(aa[jr][jc]) )fprintf(stderr,"coefficient %d %d %lf\n",jr,jc,aa[jr][jc]); 
      
      bb=(double*)calloc(2*ngrid,sizeof(double));
      for( jr=0; jr<2*ngrid; jr++ )bb[jr]=gsl_vector_get(GD,jr);
      for( jr=0; jr<2*ngrid; jr++)if( !isfinite(bb[jr]) )fprintf(stderr,"data %d %lf\n",jr,bb[jr]); 
    
      solution_l=calloc(2*ngrid,sizeof(double));
    
      t1=clock();
      stat=cgls(2*ngrid,2*ngrid,aa,bb,solution_l);
      t2=clock();
      
      fprintf(stderr,"time for %d by %d cgls: %f\n",neqn,ngrid*2,((double)t2-(double)t1)/CLOCKS_PER_SEC);
      
      fprintf(stderr,"grid: %d %d %d %d %d %lf\n",yr,mo,dy,hr,mt,sc);
      fprintf(stderr,"ML: %d %d %d %d %d %d\n",mlDstr.yr,mlDstr.mo,mlDstr.dy,mlDstr.hr,mlDstr.mt,mlDstr.sc);    
      TimeEpochToYMDHMS(fit_time+(time_t)avg_ival/2,&yr,&mo,&dy,&hr,&mt,&sc);
      fprintf(stderr,"fit_time: %d %d %d %d %d %lf\n",yr,mo,dy,hr,mt,sc);

      fprintf(vf,"%d %d %d %d %d %d\n",yr,mo,dy,hr,mt,(int)sc);
      fprintf(vf,"%f %f %f %f %f %f\n",mlDstr.Bx,mlDstr.By,mlDstr.Bz,mlDstr.v_sw,mlDstr.Au,mlDstr.Al);
      fprintf(vf,"%d\n",ngrid);
      for( jc=0; jc<ngrid; jc++)
	fprintf(vf,"%f %f %f %f %f %f %d\n",grid[jc].center_lat,grid[jc].center_lon,solution_l[jc],solution_l[jc+ngrid],rgMLDstr.ve[jc],rgMLDstr.vn[jc],grid_data_count[jc]);
    
      fflush(vf);

      fprintf(losout,"%d %d %d %d %d %d\n",yr,mo,dy,hr,mt,(int)sc);
      fprintf(losout,"%f %f %f %f %f %f\n",mlDstr.Bx,mlDstr.By,mlDstr.Bz,mlDstr.v_sw,mlDstr.Au,mlDstr.Al);
      fprintf(losout,"%d\n",ndata);
      for( jc=0; jc<ndata; jc++)
	fprintf(losout,"%f %f %f %f\n",los_lats[jc],los_lons[jc],los_data[jc],los_kazm[jc]);
      
      
      free(solution_l);
      
      for( jr=0; jr<2*ngrid; jr++)free(aa[jr]); 
      free(aa);
      free(bb);
    
      gsl_matrix_free(coef);
      gsl_matrix_free(coefTrans);
      gsl_matrix_free(CDI);
      gsl_vector_free(GD);
      gsl_vector_free(data);
      gsl_matrix_free(A);

      for( jr=0; jr<ndata; jr++)free(kazm_array[jr].coef);
    }
    
    fit_time+=(time_t)avg_ival;
    
  }
  free(div_array);
  free(grid);
  free(grid_data_count);
  free(neighbors);
  free(rgMLDstr.vn);
  free(rgMLDstr.ve);
  free(rgMLDstr.vn_cov);
  free(rgMLDstr.ve_cov);

  fclose(vf);
  fclose(fp); 
}
