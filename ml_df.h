#include "cnvgrid.h"

#define PATH_LEN 512
#define F_NAME_LEN 64
#define MXFILES 1000
#define NRAD 33
#define MAX_GATES 1000
#define MAX_BEAM 25
#define MAX_COUNT 10

#define CNV_TO_STR(num,str) if(num < 10){sprintf(str,"0%d",num);} else{sprintf(str,"%d",num);}
#define IS_LEAPYEAR(year) ((year%4)&&(!(year%100))||(year%400))

typedef struct cell{
  double lat[4];
  double lon[4];
  double center_lat;
  double center_lon;
  int lat_indx;
  int lon_indx;
  int edge;
}CELL;

typedef struct boundary{
  double lat;
  double lon;
}B_POINT;

typedef struct Corotation{
  double vc_mag;
  double vc_azm;
} COROTATION;

typedef struct Neighbor{
  int left[2];
  int right[2];
  int lower[3];
  int upper[3];
}NEIGHBOR;

typedef struct mod_array{
  double* coef;
}M_ARRAY;

typedef struct time{
  int yr;
  int mo;
  int dy;
  int hr;
  int mt;
  int sc;
  long us;
} TIME;

typedef struct finfo{
  char dir_path[PATH_LEN];
  char fname[F_NAME_LEN];
}FILE_INFO;

typedef struct RadarPos{
  double lat[MAX_BEAM][MAX_GATES];
  double lon[MAX_BEAM][MAX_GATES];
  double kazm[MAX_BEAM][MAX_GATES];
}RPOS;  

typedef struct RadarMap{
  int cell[MAX_BEAM][MAX_GATES][MAX_COUNT];
  int cell_count[MAX_BEAM][MAX_GATES];
}RMAP;  

typedef struct PotRec{
  double pot;
  int count;
  double glat;
  double glon;
  double mlat;
  double mlon;
  double mlt;
}POTEN;

/* void parse_instructions(); */
/* void find_fit_files(struct tm *,struct tm *,char *[], char *[], int *); */
struct tm *parse_date_str(long);

int make_local_pgrid(struct cell * , int, struct CnvGrid * );

typedef struct MLDSTR{
  int yr;
  int mo;
  int dy;
  int hr;
  int mt;
  int sc;
  float Bx,By,Bz,Au,Al,v_sw;
  int npts;
  double* lats;
  double* lons;
  double* vn;
  double* ve;
  double* vn_cov;
  double* ve_cov;
  double* vmag;
  double* vaz;
}MLDstr;


typedef struct GRDMLDSTR{
  double* vn;
  double* ve;
  double* vn_cov;
  double* ve_cov;
}GridMLDstr;
