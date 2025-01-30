#define MAX_NRANG 225
#define MAX_nANG 60
#define BAD_VALUE -99999.

struct RadarNetwork *network;  
struct Radar *radar;
static int filt_flag;

struct vIMAGE_STR{
  int year;
  int month;
  int day;
  int hour;
  int minut;
  int sec;
  int nsec;
  int nrang;
  int beam;
  int mpinc;
  int smsep;
  int freq;
  char radar[8];
  double angle[MAX_NRANG][MAX_nANG];
  double pwr[MAX_NRANG][MAX_nANG];
  double vel[MAX_NRANG][MAX_nANG];
  double wid[MAX_NRANG][MAX_nANG];  
};

struct RadarImPos{
  double lat[MAX_NRANG+1][MAX_nANG];
  double lon[MAX_NRANG+1][MAX_nANG];
  double kazm[MAX_NRANG+1][MAX_nANG];
};  


void fldpnth(double gdlat, double gdlon, double psi, double bore,
             double fh, double r, double *frho, double *flat,
             double *flon, int chisham);

void fldpnt_azm(double mlat, double mlon, double nlat, double nlon, double *az);

