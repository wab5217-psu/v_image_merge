#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <dirent.h>
#include <string.h>
#include <sys/stat.h>
#include <stdio.h>
#include <getopt.h>
#include <math.h>
#include "radar.h"
#include "rpos.h"
#include "geodtgc.h"
#include "v_image.h"
#include "../v_image_src/image.h"

#include "aacgmlib_v2.h"
#include "magcmp.h"

#define SPS_STID 22
#define MCM_STID 20
#define KOD_STID 7
#define BKS_STID 33

#define MIN_PWR .1


double myatan2d(double arg1, double arg2){
  double pi=acos(-1.);
  return(atan2(arg1,arg2)*180./pi);
}


void norm_vec(double *x, double *y, double *z);

/**
 * Converts from global geocentric spherical coordinates (r,theta,phi) to
 * global Cartesian coordinates (x,y,z), with input values in degrees.
 **/
void sphtocar(double r, double theta, double phi,
              double *x, double *y, double *z);


void glbthor(int iopt, double lat, double lon,
             double *rx, double *ry, double *rz,
             double *tx, double *ty, double *tz);

void geodtgc(int iopt, double *gdlat, double *gdlon,
             double *grho, double *glat,
             double *glon, double *del);

void ang_rb_filt(struct vIMAGE_STR *vstr){

  int jr,ja;
  int dcount;
  int count_thresh;

  count_thresh=0.33*vstr->nrang;
  
  for(ja=0; ja<MAX_nANG; ja++){
    dcount=0;
    for( jr=0; jr<MAX_NRANG; jr++ ){
      if( vstr->pwr[jr][ja]!=BAD_VALUE )dcount++;
    }
    if( dcount>=count_thresh )for( jr=0; jr<MAX_NRANG; jr++ ){
      vstr->pwr[jr][ja]=BAD_VALUE;
      vstr->vel[jr][ja]=BAD_VALUE;
      vstr->wid[jr][ja]=BAD_VALUE;
    }
  }

  double hold_pwr[MAX_NRANG][MAX_nANG];
  double hold_vel[MAX_NRANG][MAX_nANG];
  double hold_wid[MAX_NRANG][MAX_nANG];

  for(ja=0; ja<MAX_nANG; ja++){
    for( jr=1; jr<MAX_NRANG-1; jr++ ){
      if( vstr->pwr[jr-1][ja]==BAD_VALUE && vstr->pwr[jr+1][ja]==BAD_VALUE){
	hold_pwr[jr][ja]=BAD_VALUE;
	hold_vel[jr][ja]=BAD_VALUE;
	hold_wid[jr][ja]=BAD_VALUE;
      }else{
	hold_pwr[jr][ja]=vstr->pwr[jr][ja];
	hold_vel[jr][ja]=vstr->vel[jr][ja];
	hold_wid[jr][ja]=vstr->wid[jr][ja];
      }
    }
  }

  for(ja=0; ja<MAX_nANG; ja++)for( jr=1; jr<MAX_NRANG-1; jr++ ){
      vstr->pwr[jr][ja]=hold_pwr[jr][ja];
      vstr->vel[jr][ja]=hold_vel[jr][ja];
      vstr->wid[jr][ja]=hold_wid[jr][ja];
    }
}



int read_vImage(FILE *fp, struct vIMAGE_STR *vstr){

  int nread;
  int rng,aindx,nang;
  double ang,amp,vel,wid;
  int j0,j1;
  
  for( j0=0; j0<MAX_NRANG; j0++ )for( j1=0; j1<MAX_nANG; j1++ ){
      vstr->pwr[j0][j1]=BAD_VALUE;
      vstr->vel[j0][j1]=BAD_VALUE;
      vstr->wid[j0][j1]=BAD_VALUE;
    }
  
  if((nread=fscanf(fp,"%d %d %d %d %d %d",&(vstr->year),&(vstr->month),&(vstr->day),&(vstr->hour),&(vstr->minut),&(vstr->sec)))!=6 )return(-1);
  if((nread=fscanf(fp,"%d %d %d",&(vstr->nrang),&nang,&(vstr->smsep)))!=3 )return(-1);  
  if((nread=fscanf(fp,"%d %d %d",&(vstr->beam),&(vstr->freq),&(vstr->mpinc)))!=3 )return(-1);  
  while( (nread=fscanf(fp,"%d %d %lf %lf %lf %lf",&rng,&aindx,&ang,&amp,&wid,&vel))>0){
    vstr->pwr[rng][aindx]=amp*amp;
    vstr->wid[rng][aindx]=wid;
    vstr->vel[rng][aindx]=vel;
    vstr->angle[rng][aindx]=ang;

  }
  return(0);
}

void get_pos_ar(struct RadarSite *site, int year, int rsep, double *angs,  struct RadarImPos *rdrpos){
  
  double lat,lon,azm;
  double d,psi;
  int rn,bm,s;
  int chisham=1;
  double flat,flon,frho;
  double fx,fy,fz;
  
  double gx,gy,gz;
  double glat,glon;
  double gdlat,gdlon,gdrho;
  double gbx,gby,gbz; 
  double ghx,ghy,ghz;
  double bx,by,bz,b;
  double dummy;
  
  int rxrise=0;
  int frang=180;

  double height=350.;
  double range_edge=0;

  int center=0;

  if (center==0) {
    range_edge=-0.5*rsep*20/3;
  }

  gdlat=site->geolat;
  gdlon=site->geolon;

  
  for( bm=0; bm<MAX_nANG; bm++) for( rn=0; rn<=MAX_NRANG; rn++){      
      
      psi=angs[bm];
      d=slant_range(frang,rsep,rxrise,range_edge,rn+1);

      fldpnth(site->geolat,site->geolon,psi,site->boresite,
	      height,d,&frho,&lat,&lon,chisham);

      rdrpos->lat[rn][bm]=lat;
      rdrpos->lon[rn][bm]=lon;

      /* fldpnt_azm(llon,llat,lon,lat,&azm);  */

      flat=lat;
      flon=lon;
      
      sphtocar(frho,flat,flon,&fx,&fy,&fz);

      /* Convert radar site geodetic latitude/longitude (gdlat,gdlon) to
       * geocentric spherical coordinates (glat,glon) and distance from the
       * center to the surface of the oblate spheroid (gdrho) */
      geodtgc(1,&gdlat,&gdlon,&gdrho,&glat,&glon,&dummy);
      
      /* Convert radar geocentric coordinates (gdrho,glat,glon) to global
       * Cartesian coordinates (gbx,gby,gbz) */
      sphtocar(gdrho,glat,glon,&gbx,&gby,&gbz);
      
      /* Calculate vector from the radar to center of range/beam cell (gx,gy,gz) */
      gx=fx-gbx;
      gy=fy-gby;
      gz=fz-gbz;
      
      /* Normalize the vector from the radar to center of range/beam cell */
      norm_vec(&gx,&gy,&gz);

      glbthor(1,flat,flon,&gx,&gy,&gz,&ghx,&ghy,&ghz);

      /* Normalize the local horizontal radar-to-range/beam cell vector */

      norm_vec(&ghx,&ghy,&ghz);

      s=IGRFMagCmp(year,frho,flat,flon,&bx,&by,&bz,&b);
      if (s==-1) fprintf(stderr,"IGRF failed ****************************\n");


      /* Normalize the magnetic field vector */
      norm_vec(&bx,&by,&bz);
      
      /* Calculate a new local vertical component such that the radar-to-range/beam
       * vector becomes orthogonal to the magnetic field at the range/beam position
       * (gh dot b = 0) */
      ghz=-(bx*ghx+by*ghy)/bz;

      /* Normalize the new radar-to-range/beam vector (which is now orthogonal to B) */
      norm_vec(&ghx,&ghy,&ghz);

      
      /* Calculate the azimuth of the orthogonal radar-to-range/beam vector */
      azm=myatan2d(ghy,-ghx);
      
      rdrpos->kazm[rn][bm]=azm;

      /* llat=lat; */
      /* llon=lon; */

    }
}

struct RadarNetwork *network;  
struct Radar *radar;
static int filt_flag;

int main(int argc, char *argv[]){

  char *envstr;
  char dir[10];
  char date[24];
  DIR *dp;
  FILE *fp;
  struct dirent *entry;
  struct stat statbuf;
  struct vIMAGE_STR vstr;
  struct RadarImPos rpos;
  int stat;
  int jr,ja;
  double vel[MAX_NRANG][MAX_nANG];
  double vhold[MAX_NRANG][MAX_nANG];
  double pwr[MAX_NRANG][MAX_nANG];
  double wid[MAX_NRANG][MAX_nANG];
  double angs[MAX_nANG];
  int dcount[MAX_NRANG][MAX_nANG];
  double rsep;

  struct RadarSite *site=NULL;

  filt_flag=0;
  int c;
  while(1){
    static struct option long_options[]=
      {
	{"med_filt",no_argument, &filt_flag, 1}
      };
    int option_index=0;
    c=getopt_long(argc,argv,"abc",long_options,&option_index);
    if( c==-1 ) break;
  }

  char rchar[3];
  char rcode[4];
  strcpy(date,argv[optind]);  
  strcpy(rcode,argv[optind+1]);
  sprintf(dir,"%s",".");

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
  

  for( jr=0; jr<MAX_NRANG; jr++ )for( ja=0; ja<MAX_nANG; ja++ ){
      vel[jr][ja]=0;
      pwr[jr][ja]=0;
      wid[jr][ja]=0;
      dcount[jr][ja]=0;
    }

  if((dp = opendir(dir)) == NULL) {
    fprintf(stderr,"cannot open directory: %s\n", dir);
    exit(-1);
  }  
  while( (entry=readdir(dp)) != NULL ){
    lstat(entry->d_name,&statbuf);
    if(S_ISREG(statbuf.st_mode) && (strncmp(date,entry->d_name,12) == 0) && (strstr(entry->d_name,"vel_v_ang") != NULL)) {
      fprintf(stderr,"%s\n",entry->d_name);
      fp=fopen(entry->d_name,"r");
      if( (stat=read_vImage(fp,&vstr)) !=0 ) continue;

      for( jr=0; jr<MAX_NRANG; jr++ )for( ja=0; ja<MAX_nANG; ja++ ) if( vstr.pwr[jr][ja]!=BAD_VALUE ){
	    /* if( 10*log10(vstr.pwr[jr][ja])>80 ) fprintf(stderr,"%s %d %d %f %f\n",entry->d_name,jr,ja,vstr.vel[jr][ja],log10(vstr.pwr[jr][ja])*10); */
	    vel[jr][ja]+=vstr.vel[jr][ja];
	    pwr[jr][ja]+=vstr.pwr[jr][ja];
	    wid[jr][ja]+=vstr.wid[jr][ja];
	    dcount[jr][ja]++;
	  }
      fclose(fp);
    }
  }  

  for( jr=0; jr<MAX_NRANG; jr++ )for( ja=0; ja<MAX_nANG; ja++ ) if( dcount[jr][ja]!=0 ){
	vel[jr][ja]/=(double)dcount[jr][ja];
	pwr[jr][ja]/=(double)dcount[jr][ja];
	wid[jr][ja]/=(double)dcount[jr][ja];
      }else{
	vel[jr][ja]=BAD_VALUE;
	pwr[jr][ja]=BAD_VALUE;
	wid[jr][ja]=BAD_VALUE;
      }
  envstr=getenv("SD_HDWPATH");
  RadarLoadHardware(envstr,network);

  int stid;
  char rname[8];
  if( strcmp("kod",rcode)==0 ){
    stid=KOD_STID;
    sprintf(rchar,"d");
    sprintf(rname,"kod");
  }else if( strcmp("mcm",rcode)==0 ){
    stid=MCM_STID; /***** NEED TO FIX THIS *******/
    sprintf(rchar,"a");
    sprintf(rname,"mcm");
  }else if( strcmp("sps",rcode)==0 ){
    stid=SPS_STID; /***** NEED TO FIX THIS *******/
    sprintf(rchar,"a");
    sprintf(rname,"sps");
  }else if( strcmp("bks",rcode)==0 ){
    stid=BKS_STID; /***** NEED TO FIX THIS *******/
    sprintf(rchar,"a");
    sprintf(rname,"bks");
  }else{ 
    fprintf(stderr,"********** Invalid station code %s ***********\n",rcode);
    exit(1);
  }
    
  for( ja=0; ja<MAX_nANG; ja++ )angs[ja]=(double)ANG_MIN+(double)ja*DANG;

  rsep=(double)vstr.smsep*.15;
  radar=RadarGetRadar(network,stid);
  site=RadarYMDHMSGetSite(radar,vstr.year,vstr.month,vstr.day,vstr.hour,vstr.minut,vstr.sec);
  if( site==NULL )fprintf(stderr,"NULL site\n");
  get_pos_ar(site,vstr.year,rsep,angs,&rpos);


  /* median filter */
  
  double count;
  int jrr,jaa;
  if( filt_flag ){
    for( jr=1; jr<MAX_NRANG-1; jr++ )for( ja=1; ja<MAX_nANG-1; ja++){
	count=0.;
	vhold[jr][ja]=0;
	for( jrr=jr-1; jrr<=jr+1; jrr++ )for( jaa=ja-1; jaa<=ja+1; jaa++ ){
	    if( vel[jrr][jaa] != BAD_VALUE && pwr[jrr][jaa] >= MIN_PWR ){
	      count+=1.;
	      vhold[jr][ja]+=vel[jrr][jaa];
	    }
	  }
	if( count > 1. )vhold[jr][ja]/=count; else vhold[jr][ja]=BAD_VALUE;
      }
    
    for( jr=1; jr<MAX_NRANG-1; jr++ )for( ja=1; ja<MAX_nANG-1; ja++)vel[jr][ja]=vhold[jr][ja];
  }

  for( jr=0; jr<MAX_NRANG; jr ++ ){
    vel[jr][0]=BAD_VALUE;
    pwr[jr][0]=BAD_VALUE;
    wid[jr][0]=BAD_VALUE;
    vel[jr][1]=BAD_VALUE;
    pwr[jr][1]=BAD_VALUE;
    wid[jr][1]=BAD_VALUE;
    vel[jr][MAX_nANG-1]=BAD_VALUE;
    pwr[jr][MAX_nANG-1]=BAD_VALUE;
    wid[jr][MAX_nANG-1]=BAD_VALUE;
    vel[jr][MAX_nANG-2]=BAD_VALUE;
    pwr[jr][MAX_nANG-2]=BAD_VALUE;
    wid[jr][MAX_nANG-2]=BAD_VALUE;
  }
  
  FILE *vout;
  char outf[64];
  sprintf(outf,"%s.v_image",date);
  if( (vout=fopen(outf,"w")) == NULL ){
    fprintf(stderr,"--- FAILED TO OPEN OUTPUT FILE ---\n");
    exit(-1);
  }
  
  fprintf(vout,"%s\n",date);
  fprintf(vout,"%d %d\n",vstr.nrang,MAX_nANG);
  fprintf(vout,"%s\n",rname);
  for( jr=0; jr<vstr.nrang; jr++ )for( ja=0; ja<MAX_nANG; ja++ )
				    fprintf(vout,"%d %d %lf %lf %lf %lf %lf %lf\n",jr,ja,vel[jr][ja],pwr[jr][ja],wid[jr][ja],rpos.lat[jr][ja],rpos.lon[jr][ja],rpos.kazm[jr][ja]);

  fflush(vout);
  fclose(vout);  
}
