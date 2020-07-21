/*
*                 OpenATS
*          Open Auto Tracking System
*
*
*
*
*
*
*
*
*
*
*
*
*
*/
/*
/                       Open  Auto Tracking System
/
/  
/  Use this tool can tracking satellites's Freq and angle by 
/  the GPS module and the TLE data of satellites,output with Socket packet 
/  and COM port.this tool will very useful for the tracking sateliite receiver
/  or UAV.the tool read the data from GPS,so you need a GPS module connected befor
/  
/  Use: 
/  
/
/ 
/
/
//GPS module need support NMEA-0183 
/
/
/
/
/
/
/计算太阳角度函数
/
/
/
/
*/
#include <gps.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <getopt.h>
#include <assert.h>
#include <signal.h>
#include <ncurses.h>
#include <time.h>
#include <sys/select.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <sys/timeb.h>
#include <unistd.h>
#include <libwebsockets.h>
#include <mysql/mysql.h>
#include "socket.h"
#include "minIni.h"
#include "sgdp4h.h"
#include "serialib.h"
#include "gpsdecoder.h"
#include "cJSON.h"

#define BUFSIZE 256
#define ST_SIZE 128
#ifndef PI
#define PI 3.141592653589793
#endif
#ifndef M_PI
#define M_PI  3.141592653589793
#endif
#ifndef DEG
#define DEG(X) ((X)*180.0/PI)
#endif     

#define D2R M_PI/180.0
#define R2D 180.0/M_PI
#define XKMPER 6378.135 // Earth radius in km
#define XKMPAU 149597879.691 // AU in km
#define FLAT (1.0/298.257)

#define GPS_Y 70
#define GPS_X 10
#define SAT_Y 5
#define SAT_X 10
#define LOK_Y 5
#define LOK_X 10

#define SOCKET_BUF 200
#define EXAMPLE_RX_BUFFER_BYTES (10)

double mjd;
int fix_az = 0;           //fix az default disable,if find the '$GPHDT' will auto enable
int tracking = 1;
int tcp_error = 0;        //tcp receive monitor '1' is error
int usedb = 0;
long satno;
int ant_status;       //tracking mode report,1 is tracking satellite,2 is tracking target,3 is remote control
unsigned controlfreq = 200000;
unsigned netfreq = 1000000;
unsigned loopfreq = 100000;
int socket_fd,receive_fd,listen_fd,gps_result;
int ant_dev,gps_dev,ant_ret,gps_ret;
char tlefile[100], satname[20], control_freq[20],loop_freq[10],net_freq[10];
char gps_time[1];   //use GPS time or system time
char gpsd[2];       //use GPS module or GPSD data
char antennaid[5],receive_port[10];
orbit_t orb;
FILE *ini_file;
FILE *tlefile_fd;
MYSQL mysql;
MYSQL_RES *res = NULL;
MYSQL_ROW row;
struct tm *utc;
static struct gps_data_t gpsdata;
const char inifile[] = "options.ini";

typedef struct                 //get angle struct
{
    double az;
    double el;
    double range;
} s_azelcord;

typedef struct                 //local position
{
    double lat;
    double lng;
    double alt;
} s_geogcord;

s_geogcord location;
xyz_t sun;

gpsinfo gps_info = {           //gpsdata struct
    .lat = 0.0,
    .lng = 0.0,
    .alt = 0.0,
};     

s_geogcord target = {
    .lat = 0.0,
    .lng = 0.0,
    .alt = 0.0,
};

s_azelcord lookangles ={
    .az = 0.0,
    .el = 0.0,
    .range = 0.0,
};

typedef struct                 //antenna report
{
    int antennaid;
    float lat;
    float lng;
    float alt;
    float az;
    float el;
    int status;
    char ip[20];
} antenna;
antenna antinfo;

void help();
void gohome(int signo);
int sys_mjd();
int gps_mjd();
char * makejson(antenna *ant);
int parsejson(char * msg);
void read_gpsd(int test);
int ftime(struct timeb *tp);
void print_orb(orbit_t *orb);
void receive_server (void *ptr);
int send_message ();
static char *st_start(char *buf);
void compute_observer_and_solar_positions(void);
void sunpos_xyz(double mjd,xyz_t *pos,double *ra,double *de);
double altitude(orbit_t orb,double mjd,s_geogcord s);
static double d_read(char *str, int start, int stop);
static long i_read(char *str, int start, int stop);
int read_twoline(FILE *fp, long satno, orbit_t *orb,char *satname);
void *vector(size_t num, size_t size);
long search_satno(char *satname);
double date2mjd(int Year, int Month, double Day);
void fix_angle(s_azelcord *ant_angles, gpsinfo *gps_info);
int sat_angle(orbit_t *orb,double mjd, s_geogcord s,double *v, double *az,double *el);
int sat_tracking(orbit_t *orb);
void target_angle(s_geogcord *ground, s_geogcord *obs);
int target_tracking();
void hand_control();
void equatorial2horizontal(double mjd, s_geogcord s,double ra,double de,double *az,double *el);
int write_database(antenna *antinfo);
WINDOW *create_window(int,int,int,int);

void help ()
{
    printf("\n");
	printf("+===============================================================================+\n");
	printf("#                  ______                         _______ _____                 #\n");
	printf("#                 /  __  \\                     /\\|__   __/ ____|                #\n");
    printf("#                |  |  |  |____   ___ _ __    /  \\  | | | (___                  #\n");
	printf("#                |  |  |  |  _ \\ / _ \\ '_ \\  / /\\ \\ | |  \\___ \\                 #\n");
	printf("#                |  |__|  | |_) |  __/ | | |/ ____ \\| |  ____) |                #\n");
	printf("#                 \\______/| .__/ \\___|_| |_/_/    \\_\\_| |_____/                 #\n");
	printf("#                         | |                                                   #\n");
	printf("#                         |_|      Open Auto Tracking System                    #\n");
	printf("# _____             _____ _____ ______ _                                        #\n");
	printf("#|  __ \\     /\\    / ____|_   _|  ____| |                                       #\n");
	printf("#| |__) |   /  \\  | (___   | | | |__  | |                                       #\n");
	printf("#|  _  /   / /\\ \\  \\___ \\  | | |  __| | |                                       #\n");
	printf("#| | \\ \\  / ____ \\ ____) |_| |_| |____| |____                                   #\n");
	printf("#|_|  \\_\\/_/    \\_\\_____/|_____|______|______|       www.openats.cn             #\n");
	printf("+===============================================================================+\n");
	printf("\n");
	printf("# use:\n");
    printf("# '-h' for help\n");
    printf("# '-u' update/download TLE file\n");
    printf("# '-s' Tracking satellite by satellite name\n");
    printf("# '-n' Tracking satellite by NORAD number\n");
    printf("# '-r' Remote control\n");
    printf("# '-t' to tracking target forexample UAV or Airplane\n");
    printf("# '-i' Show satellite TLE info by NORAD number\n");
    printf("# '-g  show GPS module data\n");
	printf("\n");
}

WINDOW *create_window(int h,int w,int y,int x)
{
    WINDOW * local_win;
    local_win = newwin(h,w,y,x);
    //box(local_win,0,0);
    wrefresh(local_win);
    return local_win;
}

//print orbit
void print_orb(orbit_t *orb)
{
	printf("+===============================================================================+\n");
    printf("# SAT NAME: %s    NORAD ID: %ld\n", satname, (long)orb->satno);
    printf("---------------------------------------------\n");
    printf("# Satellite designation = %s\n",orb->desig);
    printf("# Epoch year = %d day = %.8f\n", orb->ep_year, orb->ep_day);
    printf("# Eccentricity = %.7f\n", orb->ecc);
    printf("# Equatorial inclination = %.4f deg\n", DEG(orb->eqinc));
    printf("# Argument of perigee = %.4f deg\n", DEG(orb->argp));
    printf("# Mean anomaly = %.4f deg\n", DEG(orb->mnan));
    printf("# Right Ascension of Ascending Node = %.4f deg\n", DEG(orb->ascn));
    printf("# Mean Motion (number of rev/day) = %.8f\n", orb->rev);
    printf("# First derivative of mean motion = %e\n",orb->ndot2);
    printf("# Second derivative of mean motion = %e\n",orb->nddot6);
    printf("# BSTAR drag = %.4e\n", orb->bstar);
    printf("# Orbit number = %ld\n", orb->norb);
	printf("+===============================================================================+\n");
}

void gohome(int signo)         // ctrl+c 
{
    //printf("\ninput ctrl+c now exit...\nPlease input 'enter' to exit OpenATS\n");
    tracking = 0;
    //printf("tracking is:%d\n",tracking);
}
void tcperror(int signo)
{
    tcp_error = 1;
}

static char *st_start(char *buf)
{
    if(buf == NULL) return buf;

    while(*buf != '\0' && isspace(*buf)) buf++;

    return buf;
}

void *vector(size_t num, size_t size)
{
    void *ptr;

    ptr = calloc(num, size);
    if(ptr == NULL)
        {
        printf("vector: Allocation failed %u * %u\n", num, size);
        }

    return ptr;
}

static long i_read(char *str, int start, int stop)
{
    long itmp=0;
    char *buf, *tmp;
    int ii;

    start--;    /* 'C' arrays start at 0 */
    stop--;

    tmp = buf = (char *)vector(stop-start+2, sizeof(char));

    for(ii = start; ii <= stop; ii++)
        {
        *tmp++ = str[ii];   /* Copy the characters. */
        }
    *tmp = '\0';            /* NUL terminate */

    itmp = atol(buf);       /* Convert to long integer. */
    free(buf);

    return itmp;
}

static double d_read(char *str, int start, int stop)
{
    double dtmp=0;
    char *buf, *tmp;
    int ii;

    start--;
    stop--;

    tmp = buf = (char *)vector(stop-start+2, sizeof(char));

    for(ii = start; ii <= stop; ii++)
        {
        *tmp++ = str[ii];   /* Copy the characters. */
        }
    *tmp = '\0';            /* NUL terminate */

    dtmp = atof(buf);       /* Convert to long integer. */
    free(buf);

    return dtmp;
}

/*
// Sun position
void compute_observer_and_solar_positions(void)
{
    int i;
    xyz_t obsvel;

    for (i=0;i<m.length;i++) {
    // Compute MJDs
    pt[i].mjd=m.mjd+(double) i/86400.0;
    mjd2date(pt[i].mjd,pt[i].nfd,0);

    // Observer position
    obspos_xyz(pt[i].mjd,&pt[i].obspos,&obsvel);

    // Solar position
    sunpos_xyz(pt[i].mjd,&pt[i].sunpos,&pt[i].sra,&pt[i].sde);
    equatorial2horizontal(pt[i].mjd,pt[i].sra,pt[i].sde,&pt[i].sazi,&pt[i].salt);
  }

  return;
}

void sunpos_xyz(double mjd,xyz_t *pos,double *ra,double *de)
{
    double jd,t,l0,m,e,c,r;
    double n,s,ecl;

    jd=mjd+2400000.5;
    t=(jd-2451545.0)/36525.0;
    l0=modulo(280.46646+t*(36000.76983+t*0.0003032),360.0)*D2R;
    m=modulo(357.52911+t*(35999.05029-t*0.0001537),360.0)*D2R;
    e=0.016708634+t*(-0.000042037-t*0.0000001267);
    c=(1.914602+t*(-0.004817-t*0.000014))*sin(m)*D2R;
    c+=(0.019993-0.000101*t)*sin(2.0*m)*D2R;
    c+=0.000289*sin(3.0*m)*D2R;

    r=1.000001018*(1.0-e*e)/(1.0+e*cos(m+c));
    n=modulo(125.04-1934.136*t,360.0)*D2R;
    s=l0+c+(-0.00569-0.00478*sin(n))*D2R;
    ecl=(23.43929111+(-46.8150*t-0.00059*t*t+0.001813*t*t*t)/3600.0+0.00256*cos(n))*D2R;

    *ra=atan2(cos(ecl)*sin(s),cos(s))*R2D;
    *de=asin(sin(ecl)*sin(s))*R2D;

    pos->x=r*cos(*de*D2R)*cos(*ra*D2R)*XKMPAU;
    pos->y=r*cos(*de*D2R)*sin(*ra*D2R)*XKMPAU;
    pos->z=r*sin(*de*D2R)*XKMPAU;

    return;
}
*/

long gps2mjd (long gpscycle, long gpsweek, long gpsseconds)
{
    long gpsdays;

    gpsdays = ((gpscycle * 1024) + gpsweek) * 7 + (gpsseconds / 86400);
    return date2mjd(1980, 1, 6) + gpsdays;
}

// nfd2mjd add "ms"
int sys_mjd()
{
    double ms, sec;
    double dday;
    //struct tm *ptm;
    struct timeb tp;

    ftime(&tp);
    //printf("%d\n",tp.time);
    utc = gmtime(&tp.time);
    //ms = (double)tp.millitm/1000;
    //printf("%f\n",ms);
    ms = utc->tm_sec + ((double)tp.millitm/1000);
    //printf("sec is %f\n",sec);
    //printf("year:%d  mon:%d  day:%d  hor:%d   min:%d  sec:%d\n",utc->tm_year+1900,utc->tm_mon+1,utc->tm_mday,utc->tm_hour,utc->tm_min,utc->tm_sec);
    dday=utc->tm_mday+(utc->tm_hour/24.0)+(utc->tm_min/1440.0)+(ms/86400.0);
    //printf("%d,%d,%.6f\n",utc.tm_year+1900,utc.tm_mon,dday);
    mjd=date2mjd(utc->tm_year+1900,utc->tm_mon+1,dday);
    //printf("mjd is %.6f\n",mjd);
    return 0;
}


double date2mjd(int year,int month,double day)
{
    int a,b;
    double jd;

    if (month<3) {
        year--;
        month+=12;
    }

    a=floor(year/100.);
    b=2.-a+floor(a/4.);

    if (year<1582) b=0;
    if (year==1582 && month<10) b=0;
    if (year==1582 && month==10 && day<=4) b=0;

    jd=floor(365.25*(year+4716))+floor(30.6001*(month+1))+day+b-1524.5;

    return jd-2400000.5;
}

double modulo(double x,double y)
{
    x=fmod(x,y);
    if (x < 0.0) x += y;

    return x;
}

// Greenwich Mean Sidereal Time
double gmst(double mjd)
{
    double t,gmst;
    t=(mjd-51544.5)/36525.0;
    gmst = modulo(280.46061837+360.98564736629*(mjd-51544.5)+t*t*(0.000387933-t/38710000),360.0);
    
    return gmst;
}

// Greenwich Mean Sidereal Time
double dgmst(double mjd)
{
    double t,dgmst;
    t = (mjd-51544.5)/36525.0;
    dgmst = 360.98564736629+t*(0.000387933-t/38710000);

    return dgmst;
}

//red the TLE data and return &orb
int read_twoline(FILE *fp, long search_satno, orbit_t *orb, char *satname)
{
    static long search;
    static char line0[ST_SIZE],line1[ST_SIZE],line2[ST_SIZE];
    char stTemp[ST_SIZE] = "\0";
    char *st0,*st1,*st2;
    int found = 0,nLen;
    double bm, bx;

    // Set defaults
    strcpy(orb->desig,"");

    st0 = line0;
    st1 = line1;
    st2 = line2;
  /*
    do {
    if(fgets(line0,SIZE-1, fp) == NULL)
      return -1;
    st0 = st_start(line0);
    } while(st0[0] != '0');
    sprintf(satname,"%s",st0+2);
    printf("the satname is[%s]\n",satname);
  */


    /*if(search_satno > 0) {        //search search_satno data
      found = 0;

    } else {                      //search all satno data
      found = 1;
      search_satno = atol(st1+2);
    }

  */
  //printf("the search_satno is %ld\n",search_satno);
  //sprintf(search, "%ld", search_satno);
  //printf("the search target is %s\n",search);

    do {
        if(search_satno == 0)
        {
            found = 1;
            search_satno = atol(st1+2);
        }
    //printf("the search_satno is %ld\n",search_satno);
    do {
        if(fgets(line0,ST_SIZE-1, fp) == NULL)
        {
            //printf("gets line error\n");
            return -1;
        }
    } while(line0[0] != '0');
    st0 = st_start(line0);
    //printf("st0 is :%s\n",st0);
    nLen = strlen(st0);
    for(int i=0,j=0; i<nLen; i++)
    {
        if(st0[i] != '\n')
        {
            stTemp[j] = st0[i];
            j++;
        }
    }
    memset(st0,'\0', nLen);
    strcpy(st0, stTemp);
    memset(stTemp,'\0',ST_SIZE);
    //printf("st0 is :%s\n",st0);
    sprintf(satname,"%s",st0+2);
    //printf("satname is :%s\n",satname);
    //printf("----------------\n");

    if(fgets(line1, ST_SIZE-1, fp) == NULL)
        return -1;
    st1 = st_start(line1);
    search = atol(st1+2);
    //printf("search is %ld\n",search);
    if(search == search_satno) {
        found = 1;
        break;
    }
    } while(found != 1);
  //printf("############################################################################\n");
  //printf("%s\n",satname);
  //} while(fgets(line1, SIZE-1, fp) != NULL); 

    /* 
    do {
        st1 = st_start(line1);
        if(strcmp(st1+2, search) == 0) {
        found = 1;
        printf("found 1\n");
        break;
    }
    } while(fgets(line1, SIZE-1, fp) != NULL);
*/
    //sprintf(search, "2 %05ld", search_satno);
  
    if(found) 
    {
        fgets(line2, ST_SIZE-1, fp);
        st2 = st_start(line2);
    }
    //printf("TLE data:\n");
    //printf("%s",st1);
    //printf("%s\n",st2);
  
    orb->ep_year = (int)i_read(st1, 19, 20);

    if(orb->ep_year < 57) orb->ep_year += 2000;
    else orb->ep_year += 1900; 
    orb->ep_day = d_read(st1, 21, 32);  
    orb->ndot2 = d_read(st1, 34, 43);
    bm = d_read(st1, 45, 50) * 1.0e-5;
    bx = d_read(st1, 51, 52);
    orb->nddot6 = bm * pow(10.0, bx);
    bm = d_read(st1, 54, 59) * 1.0e-5;
    bx = d_read(st1, 60, 61);
    orb->bstar = bm * pow(10.0, bx);
    orb->eqinc = RAD(d_read(st2,  9, 16));
    orb->ascn = RAD(d_read(st2, 18, 25));
    orb->ecc  = d_read(st2, 27, 33) * 1.0e-7;
    orb->argp = RAD(d_read(st2, 35, 42));
    orb->mnan = RAD(d_read(st2, 44, 51));
    orb->rev  = d_read(st2, 53, 63);
    orb->norb = i_read(st2, 64, 68);
    orb->satno = search;
    //sscanf(st1+9,"%s",orb->desig);
    strncpy(orb->desig,st1+9,8);
    orb->desig[8]='\0';
    return 0;
}

//search satno by satname
long search_satno (char *satname)
{
    char zeroname[46],zero[1],zero_one[15],zero_two[15],zero_three[15];
    char line0[ST_SIZE],line1[ST_SIZE],line2[ST_SIZE],stTemp[ST_SIZE] = "\0";    
    char *st0,*st1,*st2;  
    FILE * tle_fd;

    st0 = line0;
    st1 = line1;
    st2 = line2;

    //open TLE file
	if ((tle_fd = fopen(tlefile,"r")) == NULL)     //open with read only
	{
        printf("can't open TLE file! please check your cmd or file\n");
		exit(0);
    }

    do{ 
        if(fgets(line0,ST_SIZE,tle_fd) == NULL)
        {
            printf("can't get tledata from file.\n");
            return -1;
        }
        st0 = st_start(line0);
    }while(st0[0] != '0');
    st0 = st_start(line0);
    int nLen = strlen(st0);
    for(int i = 0,j = 0; i < nLen; i++)
    {
      if(st0[i] != '\n')
      {
        stTemp[j] = st0[i];
        j++;
      }
    }
    memset(st0,'\0', nLen);
    strcpy(st0, stTemp);
    sscanf(line0,"%s %s %s %s",zero,zero_one,zero_two,zero_three);
    //printf("zeroone is %s\n",zero_one);
    //printf("zerotwo is %s\n",zero_two);
    sprintf(zeroname,"%s %s %s",zero_one,zero_two,zero_three);          
    while (strncmp(zeroname,satname,strlen(satname)) != 0)
    {
        do{ 
            if(fgets(line0,ST_SIZE,tle_fd) == NULL)
            {
              printf("can't find the %s\n",satname);
              return -1;
            }
        }while(line0[0] != '0');
        st0 = st_start(line0);
        sscanf(line0,"%s %s %s %s",zero,zero_one,zero_two,zero_three);
        sprintf(zeroname,"%s %s %s",zero_one,zero_two,zero_three);
        //printf("searching the zeroname is %s\n",zeroname);
    }
    //printf("the satellite name is:[%s]\n",satname);
    if(fgets(line1,ST_SIZE,tle_fd) == NULL)
    {
        printf("can't find the %s TLE data\n",satname);
        return -1;
    }
    st1 = st_start(line1);
    long satno = atoi(st1+2);
    //printf("the satellite number is:[%ld] \n",satno);

    return satno;
}

int sat_angle (orbit_t *orb, double mjd, s_geogcord s, double *v, double *az, double *el)
{
    double dx,dy,dz,dvx,dvy,dvz,r;
    double ra,de;
    xyz_t satpos,obspos,satvel,obsvel;

    //get卫星和天线的地心XYZ坐标
    obspos_xyz(mjd,s,&obspos,&obsvel);
    satpos_xyz(mjd+2400000.5,&satpos,&satvel);


    dx = satpos.x-obspos.x;
    dy = satpos.y-obspos.y;
    dz = satpos.z-obspos.z;
    dvx = satvel.x-obsvel.x;
    dvy = satvel.y-obsvel.y;
    dvz = satvel.z-obsvel.z;

    r = sqrt(dx*dx+dy*dy+dz*dz);
    *v = (dvx*dx+dvy*dy+dvz*dz)/r;

    ra = modulo(atan2(dy,dx)*R2D,360.0);
    de = asin(dz/r)*R2D;

    equatorial2horizontal(mjd,s,ra,de,az,el);

    return 0;
}

void equatorial2horizontal(double mjd, s_geogcord s,double ra,double de,double *az,double *el)
{
    double h;
    double azi;

    h = gmst(mjd)+s.lng-ra;

    azi = modulo(atan2(sin(h*D2R),cos(h*D2R)*sin(s.lat*D2R)-tan(de*D2R)*cos(s.lat*D2R))*R2D,360.0);
    *el = asin(sin(s.lat*D2R)*sin(de*D2R)+cos(s.lat*D2R)*cos(de*D2R)*cos(h*D2R))*R2D;

    if(azi >= 180)
        *az = azi-180;
    else if(azi < 180)
        *az = azi+180;

    return;
}

// Observer position
void obspos_xyz(double mjd, s_geogcord s, xyz_t *pos, xyz_t *vel)
{
    double ff,gc,gs,theta,sl,dtheta;
    sl = sin(s.lat*D2R);
    ff = sqrt(1.0-FLAT*(2.0-FLAT)*sl*sl);
    gc = 1.0/ff+s.alt/XKMPER;
    gs = (1.0-FLAT)*(1.0-FLAT)/ff+s.alt/XKMPER;
    theta = gmst(mjd)+s.lng;
    dtheta = dgmst(mjd)*D2R/86400;
    pos->x = gc*cos(s.lat*D2R)*cos(theta*D2R)*XKMPER;
    pos->y = gc*cos(s.lat*D2R)*sin(theta*D2R)*XKMPER; 
    pos->z = gs*sin(s.lat*D2R)*XKMPER;
    vel->x = -gc*cos(s.lat*D2R)*sin(theta*D2R)*XKMPER*dtheta;
    vel->y = gc*cos(s.lat*D2R)*cos(theta*D2R)*XKMPER*dtheta; 
    vel->z = 0.0;

    return;
}


// SGDP4 algitude  暂时未用 
double altitude(orbit_t orb,double mjd,s_geogcord s)
{
    double dx,dy,dz,dvx,dvy,dvz,r;
    double ra,de,azi,alt;
    xyz_t satpos,obspos,satvel,obsvel;

    // Loop over data points
    obspos_xyz(mjd,s,&obspos,&obsvel);
    satpos_xyz(mjd+2400000.5,&satpos,&satvel);

    dx = satpos.x-obspos.x;  
    dy = satpos.y-obspos.y;
    dz = satpos.z-obspos.z;
    r = sqrt(dx*dx+dy*dy+dz*dz);

    ra = modulo(atan2(dy,dx)*R2D,360.0);
    de = asin(dz/r)*R2D;

    equatorial2horizontal(mjd,s,ra,de,&azi,&alt);

    return alt;
}

int write_database(antenna *antinfo)      //write into database
{
    char query_str[1000];
    int i,n;

    printf("the table name is :%d\n",antinfo->antennaid);
    //printf("the satename is :%s",satname);
    //printf("satname leth is :%ld\n",strlen(satname));
    res = mysql_list_tables(&mysql,antinfo->antennaid);

    if((row=mysql_fetch_row(res)) != NULL)   //already have table
    {

        //写入数据库格式：
        //          天线编号，天线维度，天线经度，天线高度，天线方位角，天线仰角，天线状态
        //字段名字：antenna, lat, lng, alt, az, el, status


        printf("the table :[%s] already have\n",antinfo->antennaid);
        sprintf(query_str,"INSERT INTO `%d` (antenna,lat,lng,alt,az,el,status) VALUES (\
%d,\
%.6f,\
%.6f,\
%.2f,\
%.4f,\
%.4f,\
%d,\
%d\
);",antinfo->antennaid,antinfo->antennaid,antinfo->lat,antinfo->lng,antinfo->alt,antinfo->az,antinfo->el,antinfo->status,antinfo->ip);
        //printf("the query_str is :%s\n",query_str);
        i = mysql_real_query(&mysql,query_str,strlen(query_str));
        if(i != 0)
        {
            printf("write fail\n");
            printf("mysql_real_query :%s\n",mysql_error(&mysql));
            return -1;
        }
        printf("write data success!\n");
    }
    else if(row == NULL)    //don't have the satellite-id table
    {
        printf("the table :%d not have,now new-built\n",antinfo->antennaid);
        sprintf(query_str,"CREATE TABLE `%s` (\
Update_Time timestamp NOT NULL ON UPDATE CURRENT_TIMESTAMP DEFAULT CURRENT_TIMESTAMP,\
antenna int(10),\
lat float(10,6),\
lng float(10,6),\
alt float(8,2),\
az float(10,4),\
el float(10,4),\
status int(10),\
ip char(20)\
);",antinfo->antennaid);
        //printf("the mysqlcmd is:\n%s\n",query_str);
        n = mysql_real_query(&mysql,query_str,strlen(query_str));
        if(n != 0)
        {
            printf("mysql_real_query :%s\n",mysql_error(&mysql));
            return -1;          
        }
        //printf("the table :%s already CREATE success\n",table);

        sprintf(query_str,"INSERT INTO `%d` (antenna,lat,lng,alt,az,el,status) VALUES (\
%d,\
%.6f,\
%.6f,\
%.2f,\
%.4f,\
%.4f,\
%d,\
%d\
);",antinfo->antennaid,antinfo->antennaid,antinfo->lat,antinfo->lng,antinfo->alt,antinfo->az,antinfo->el,antinfo->status,antinfo->ip);
        printf("write data success!\n");
        i = mysql_real_query(&mysql,query_str,strlen(query_str));
        if(i != 0)
        {
            printf("mysql_real_query :%s\n",mysql_error(&mysql));
            return -1;
        }
    }
    return 0;
}

char * makejson(antenna *ant)
{

    cJSON * pJsonRoot = NULL;

    pJsonRoot = cJSON_CreateObject(); //新建一个JSON项目：pSubJson
    if(NULL == pJsonRoot)
    {
        //error happend here
        return NULL;
    }

    //add 字符串、数字和bool变量
    cJSON_AddNumberToObject(pJsonRoot, "antennaid", ant->antennaid);
    cJSON_AddNumberToObject(pJsonRoot, "lat", ant->lat);
    cJSON_AddNumberToObject(pJsonRoot, "lng", ant->lng);
    cJSON_AddNumberToObject(pJsonRoot, "alt", ant->alt);
    cJSON_AddNumberToObject(pJsonRoot, "az", ant->az);
    cJSON_AddNumberToObject(pJsonRoot, "el", ant->el);
    cJSON_AddNumberToObject(pJsonRoot, "status", ant->status);
    cJSON_AddStringToObject(pJsonRoot, "ip", ant->ip);

    char * p = cJSON_Print(pJsonRoot);
    //else use : 
    //char * p = cJSON_PrintUnformatted(pJsonRoot);
    if(NULL == p)
    {
        //convert json list to string faild, exit
        //because sub json pSubJson han been add to pJsonRoot, so just delete pJsonRoot, if you also delete pSubJson, it will coredump, and error is : double free
        cJSON_Delete(pJsonRoot);
        return NULL;
    }
    //free(p);

    cJSON_Delete(pJsonRoot);

    return p;
}


int parsejson(char * msg)
{

    if(NULL == msg)
    {
        return 0;
    }
    cJSON * pJson = cJSON_Parse(msg);
    if(NULL == pJson)                                                                                         
    {
        // parse faild, return
        return 0;
    }
   
    // get number from json
    cJSON * pSub = cJSON_GetObjectItem(pJson, "antennaid");
    if(NULL == pSub)
    {
        //get number from json faild
        return 0;
    }
    printf("lat : %d\n", pSub->valueint);

    // get number from json
    pSub = cJSON_GetObjectItem(pJson, "target");
    if(NULL == pSub)
    {
        //get number from json faild
        return 0;
    }
    printf("lat : %d\n", pSub->valueint);

    cJSON_Delete(pJson);
    return 1;
}


void read_gpsd(int test)
{
    double sec;
    double dday;
    char readgpsd_freq[10];

    ini_gets("", "readgpsd_freq","dummy", readgpsd_freq, sizeof(readgpsd_freq), inifile);

    if (gps_open(GPSD_SHARED_MEMORY, NULL, &gpsdata) != 0)
    {
	    printf("gpsd not running\n");
        tracking = 0;
        return;
    }
    while(tracking)
    {
        gps_read(&gpsdata,NULL,0);

        location.lat = gpsdata.fix.latitude;
        location.lng = gpsdata.fix.longitude;
        location.alt = gpsdata.fix.altHAE;

        if(atoi(gps_time))
        {
            //printf("%d\n",gpstime.tv_sec);
            //printf("%ld\n",gpstime.tv_nsec);
            utc = gmtime(&gpsdata.fix.time.tv_sec);
            sec = utc->tm_sec + (double)gpsdata.fix.time.tv_nsec/1000000;
            //printf("%d %d %d\n",utc->tm_year+1900,utc->tm_mon +1, utc->tm_mday);
            dday=utc->tm_mday+(utc->tm_hour/24.0)+(utc->tm_min/1440.0)+(sec/86400.0);
            //printf("%d,%d,%.6f\n",utc.tm_year+1900,utc.tm_mon,dday);
            mjd=date2mjd(utc->tm_year+1900,utc->tm_mon+1,dday);

            gps_info.year = utc->tm_year + 1900;
            gps_info.month = utc->tm_mon +1;
            gps_info.day = utc->tm_mday;
            gps_info.hour = utc->tm_hour;
            gps_info.min = utc->tm_min;
            gps_info.sec = utc->tm_sec;

        }

        gps_info.speed = gpsdata.fix.speed;
        gps_info.lat = gpsdata.fix.latitude;
        gps_info.lng = gpsdata.fix.longitude;
        gps_info.alt = gpsdata.fix.altHAE;

        if(test)
        {
            printf("latitude: %.7f°, longitude: %.7f°, alt:%.4f speed: %.6f km/h\n",gpsdata.fix.latitude,gpsdata.fix.longitude,gpsdata.fix.altHAE, 3.6 * gpsdata.fix.speed);
            //printf("latitude: %.7f°, longitude: %.7f°, alt:%.4f \n",location.lat,location.lng,location.alt);
            printf("mjd:%.6f\n",mjd);
            usleep(100000);
        }
        usleep(atoi(readgpsd_freq)*1000);
    }
    gps_close(&gpsdata);
}


//get GPS data from module 
void read_gps(int test)
{
	serial *gps_s;
	int listen_fd, socket_fd, rc;
	char gpgga_array[10] = "$GPGGA";
    char gprmc_array[10] = "$GPRMC";
    char gpgsv_array[10] = "$GPGSV";
    char gpgll_array[10] = "$GPGLL";
    char gphdt_array[10] = "$GPHDT";    //fix az angle
	char buffer[BUFSIZE],readgps_freq[10];
    char gps_port[50], gps_baudrate[10];
	struct timeval tv;
    fd_set rset; 

	bzero(buffer,sizeof(buffer));

    ini_gets("", "gps_port","dummy", gps_port, sizeof(gps_port), inifile);
    ini_gets("", "gps_baudrate","dummy", gps_baudrate, sizeof(gps_baudrate), inifile);
    ini_gets("", "readgps_freq","dummy", readgps_freq, sizeof(readgps_freq), inifile);
    unsigned int readgps = atoi(readgps_freq)*1000;

	if (serial_open(&gps_s, gps_port, atoi(gps_baudrate)) != 0){
		printf("open gps port error\n");
        //结束此进程
	}
    //fcntl(gps_s->fd, F_SETFL, FNDELAY);
    //printf("gps port opened\n");

    tv.tv_sec = 0;
    tv.tv_usec = 0;

    while(tracking){
        FD_ZERO(&rset);
        FD_SET(gps_s->fd,&rset);
        rc = select(gps_s->fd+1,&rset,NULL,NULL,&tv);
        if (rc > 0)
        {
            bzero(buffer,sizeof(buffer));
            serial_read(gps_s, buffer, '\n', 500);
            
            if(test == 1)      //-g show GPS data
            {

                
                //struct timeb tp;
                //ftime(&tp);
                //printf("start:%d\n",tp.millitm);
                printf("%s",buffer);
                if(strlen(buffer) > 20)
                {
                    if(strncmp(buffer,gpgga_array,6) == 0)
		            {
                        printf("get gpgga:%s\n",buffer);
                        gpgga(&gps_info, buffer);
                        //printf("=========================================================================================\n");            
			            //printf("           GPS lat:%.7f   GPS lng:%.7f   GPS alt:%.2f  GPS lock:%d\n",gps_info.lat, gps_info.lng, gps_info.alt, gps_info.sat);
			            //printf("=========================================================================================\n");
                    }
                    if(strncmp(buffer,gprmc_array,6) == 0)
                    {
                        printf("get gprmc:%s\n",buffer);
                        gprmc(&gps_info, buffer);
                    }
                    //ftime(&tp);
                    //printf("stop:%d\n",tp.millitm);
                    /*
                    printf("==============================================================================\n");
                    printf("gps time is :%02d/%02d/%02dT%02d:%02d:%02d\n",gps_info.year,gps_info.month,gps_info.day,gps_info.hour,gps_info.min,gps_info.sec);
                    printf("gps lat is:%f    lng:%f    alt:%f\n",gps_info.lat,gps_info.lng,gps_info.alt);
                    printf("gps valide is:%d\n",gps_info.valid);
                    printf("gps sepped is:%f\n",gps_info.speed);
                    printf("gps sat is:%d\n",gps_info.sat);
                    printf("==============================================================================\n");
                    */
                }
                continue;
                
            }
		    if(strncmp(buffer,gpgga_array,6) == 0)
		    {
                
                //printf("gpgga:\n%s\n",buffer);
                gpgga(&gps_info, buffer);
                //printf("=========================================================================================\n");            
			    //printf("           GPS lat:%.7f   GPS lng:%.7f   GPS alt:%.2f  GPS lock:%d\n",gps_info.lat, gps_info.lng, gps_info.alt, gps_info.sat);
			    //printf("=========================================================================================\n");
			    
            }
            if(strncmp(buffer,gprmc_array,6) == 0)
            {
                gprmc(&gps_info, buffer);
            }
            /*
            if(strncmp(buffer,gpggsv_array,6) == 0)
            {
                //printf("gpgsv is:%s\n",buffer);
                //gpggsv(&gps_info, buffer);
            }
            if(strncmp(buffer,gpgll_array,6) == 0)
            {
                //printf("gpgll is:%s\n",buffer);
                //gpgll(&gps_info, buffer);
            }
            */
            if(strncmp(buffer,gphdt_array,6) == 0)
            {
                //printf("gphdt is:%s\n",buffer);
                gphdt(&gps_info, buffer);
                fix_az = 1;
                //printf("fix angle is:%.2f\n",gps_info.offset);
            }

        }
        usleep(readgps);
	};

    //return;
}

enum protocols
{
	PROTOCOL_HTTP = 0,
	PROTOCOL_EXAMPLE,
	PROTOCOL_COUNT
};

static int callback_http( struct lws *wsi, enum lws_callback_reasons reason, void *user, void *in, size_t len )
{
	switch( reason )
	{
		case LWS_CALLBACK_HTTP:
			lws_serve_http_file( wsi, "example.html", "text/html", NULL, 0 );
			break;
		default:
			break;
	}

	return 0;
}

struct payload
{
	unsigned char data[LWS_SEND_BUFFER_PRE_PADDING + EXAMPLE_RX_BUFFER_BYTES + LWS_SEND_BUFFER_POST_PADDING];
	size_t len;
} received_payload;

static int callback_example( struct lws *wsi, enum lws_callback_reasons reason, void *user, void *in, size_t len )
{
	switch( reason )
	{
		case LWS_CALLBACK_RECEIVE:
			memcpy( &received_payload.data[LWS_SEND_BUFFER_PRE_PADDING], in, len );
			received_payload.len = len;
			lws_callback_on_writable_all_protocol( lws_get_context( wsi ), lws_get_protocol( wsi ) );
			break;

		case LWS_CALLBACK_SERVER_WRITEABLE:
			lws_write( wsi, &received_payload.data[LWS_SEND_BUFFER_PRE_PADDING], received_payload.len, LWS_WRITE_TEXT );
			break;

		default:
			break;
	}

	return 0;
}


static struct lws_protocols protocols[] =
{
	/* The first protocol must always be the HTTP handler */
	{
		"http-only",   /* name */
		callback_http, /* callback */
		0,             /* No per session data. */
		0,             /* max frame size / rx buffer */
	},
	{
		"example-protocol",
		callback_example,
		0,
		EXAMPLE_RX_BUFFER_BYTES,
	},
	{ NULL, NULL, 0, 0 } /* terminator */
};

//mark 大改，改成打包json给前端展示
//send data to web
int send_message ()
{
    struct lws_context_creation_info info;
	memset( &info, 0, sizeof(info) );

	info.port = 8000;
	info.protocols = protocols;
	info.gid = -1;
	info.uid = -1;

	struct lws_context *context = lws_create_context( &info );

	while(tracking)
	{
		lws_service( context, /* timeout_ms = */ 1000000 );
	}
    printf("hello,websockets done\n");
	lws_context_destroy( context );

    return 0;
    pthread_exit((void*) 0);
}

//receive antenna info from client station(TCP server mode)
void receive_message (void *ptr)
{
    int rc;
    char sendbuf[SOCKET_BUF],recvbuf[SOCKET_BUF];
    fd_set rset,wset;
    struct timeval tv;
    int port,len;
    char ip[20];
    int socket_fd = (int)ptr;

    ini_gets("", "net_freq","dummy", net_freq, sizeof(net_freq), inifile);
    netfreq = atoi(net_freq)*1000;

    tv.tv_sec = 0;
    tv.tv_usec = 0;

    while(tracking)
    {
        bzero(recvbuf,sizeof(recvbuf));
        /*
        rc = select(socket_fd+1,&rset,&wset,NULL,&tv);

        if(rc > 0)
        {   
            if(FD_ISSET(socket_fd,&wset));
            {
                //sprintf(sendbuf,"%s %s %s",antennaid,lookangles.az,lookangles.el);
                //printf("%s\n",sendbuf);
                //sprintf(sendbuf,"ANTENNAID:%s LAT:%.4F LNG:%.4F ALT:%.4F AZ:%.4f EL:%.4f STATUS:%d",antennaid,location.lat,location.lng,location.alt,lookangles.az,lookangles.el,ant_status);
                //sprintf(sendbuf,"hello");
                //tcp_send(socket_fd, sendbuf, sizeof(sendbuf));
                //printf("Send:%s\n",sendbuf);
            }
            if(FD_ISSET(socket_fd,&rset));
            {
                */
                //printf("socket receive ready\n");
        len = tcp_recv(socket_fd, recvbuf, 200);
        if(len == 0)
            break;
        printf("receive from server:%s\n",recvbuf);
        sscanf(recvbuf,"ANTENNAID:%d LAT:%f LNG:%f ALT:%f AZ:%f EL:%f STATUS:%d",&antinfo.antennaid,&antinfo.lat,&antinfo.lng,&antinfo.alt,&antinfo.az,&antinfo.el,&antinfo.status);
        printf("antid:%d,antlat:%f,antlng:%f,antalt:%f,antaz:%f,antel:%f,antstatus:%d\n",antinfo.antennaid,antinfo.lat,antinfo.lng,antinfo.alt,antinfo.az,antinfo.el,antinfo.status);
        if(usedb)
        {
            write_database(&antinfo); 
        }
        char *json = makejson(&antinfo);
        printf("json:%s\n",json);


    
    };
    tcp_close(socket_fd);
    pthread_exit((void*) 0);
}

//control antenna
void control_antenna(void *ptr)
{
	serial *ant_s;
	int rc;
	char sendbuf[BUFSIZE],recvbuf[BUFSIZE];
    char ant_port[ST_SIZE], ant_baudrate[ST_SIZE], minel[10], maxel[10];
    //char minel[10], maxel[10];
    s_azelcord ant_angles;
	struct timeval tv;
    fd_set rset,wset;

	bzero(recvbuf,sizeof(recvbuf));
    bzero(sendbuf,sizeof(sendbuf));

    ini_gets("", "ant_port","dummy", ant_port, sizeof(ant_port), inifile);
    ini_gets("", "ant_baudrate","dummy", ant_baudrate, sizeof(ant_baudrate), inifile);
    //ini_gets("", "minel", "dummy", minel, sizeof(minel), inifile); 
    //ini_gets("", "maxel", "dummy", maxel, sizeof(maxel), inifile); 
    ini_gets("", "control_freq", "dummy", control_freq, sizeof(control_freq), inifile);
    ini_gets("", "minel", "dummy", minel, sizeof(minel), inifile); 
    ini_gets("", "maxel", "dummy", maxel, sizeof(maxel), inifile); 

    controlfreq = atoi(control_freq)*1000;

	if(serial_open(&ant_s, ant_port, atoi(ant_baudrate)) != 0)
    {
		printf("open antenna port error\n");
        //报错并结束程序
        return;
	}
    //fcntl(ant_s->fd, F_SETFL, FNDELAY);
    //serial_write(ant_s, data);  //antenna reset to home

    tv.tv_sec = 0;
    tv.tv_usec = 0;

	while(tracking){

        if(lookangles.el <= atoi(minel) || lookangles.el >= atoi(maxel))
        {
            //printf("satellite not in range\n");
            //tracking = 0;
            break;
        }
	    bzero(recvbuf,sizeof(recvbuf));
        bzero(sendbuf,sizeof(sendbuf));

        FD_ZERO(&rset);
        FD_ZERO(&wset);
        FD_SET(ant_s->fd,&rset);
        FD_SET(ant_s->fd,&wset);
                
        ant_angles.az = lookangles.az;
        ant_angles.el = lookangles.el;

        //方位角修正函数，采集高端GPS的天线朝向角度，做天线角度修正，如果获取角度信息失败，则默认不使用
        //检测GPS模块输出的$GPHDT语句中的朝向角度，默认为NULL，注意，GPS模块需要支持双天线测向功能
        //仅动中通状态下使用，一般为高端GPS模块所具有的精准测向，普通GPS模块的航向角度不准确，不作为参考数据

        if(fix_az)
        {
            //printf("find fix_angle and fix az\n");
            fix_angle(&ant_angles, &gps_info);
        }
        rc = select(ant_s->fd+1,&rset,&wset,NULL,&tv);
        if(rc >0)
        {
            
            if(FD_ISSET(ant_s->fd,&rset));
            {
                //printf("ant read already\n");
                serial_read(ant_s,recvbuf,'\n',50);
                //printf("receive : %s\n",recvbuf);
            }
            
            if(FD_ISSET(ant_s->fd,&wset));
            {
                //printf("ant write already\n");
                sprintf(sendbuf,"%.2f %.2f",ant_angles.az, ant_angles.el);
                //printf("antenna angle %s\n lenis %d",sendbuf,sizeof(sendbuf));
                serial_write(ant_s, sendbuf);
                //printf("antenna:%s\n",sendbuf);
            }

            //printf("buffer form server :%s\n",buffer);
        }
        usleep(controlfreq);
	};

    //printf("stop control antenna\n");
    strcpy(sendbuf, "0");
    if(serial_write(ant_s, sendbuf) == 0)  //antenna reset to home
    {
        //printf("control antenna go home ok\n");
        //return;
    }else{
        //printf("can't control antenna go home\n");
    }
    sleep(1);
    //return;
}

//import function for OpenATS server
int remote_control()
{
    int socket_fd,listen_fd;
    int rc;
    int first = 1;  //first tracking to read TLE file
    int sat_tracking = 0;   //tracking satellite true or not
    char receive_port[10];
    char sendbuf[SOCKET_BUF], recvbuf[SOCKET_BUF];
    char dbhost[50],dbuser[20],dbpasswd[30],dbname[10],dbuse[1];
    struct timeval tv;
    fd_set rset, wset;
    pthread_t gpsd_tid, gps_tid, receive_tid, send_tid;
    WINDOW * track_win;
    void *status;
    int port;
    char ip[128];    

    ini_gets("", "receive_port","dummy", receive_port, sizeof(receive_port), inifile);
    ini_gets("", "net_freq","dummy", net_freq, sizeof(net_freq), inifile);
    ini_gets("", "database","dummy", dbuse, sizeof(dbuse), inifile);    
    ini_gets("", "gpsd","dummy", gpsd, sizeof(gpsd), inifile);

    netfreq = atoi(net_freq)*1000;
    usedb = atoi(dbuse);

    //create read_gps course
    if(atoi(gpsd))
    {
        pthread_create(&gpsd_tid, NULL, (void*) &read_gpsd, NULL);
        pthread_detach(&gpsd_tid);
    }else{
        pthread_create(&gps_tid, NULL, (void*) &read_gps, NULL);
        pthread_detach(&gps_tid);
    }
    pthread_create(&send_tid, NULL, (void*) &send_message, NULL);
    pthread_detach(&send_tid);

    do{   
        listen_fd = tcp_server("0.0.0.0", receive_port); 
    }while(listen_fd < 0);

    if(usedb)
    {
        ini_gets("", "database_name","dummy", dbname, sizeof(dbname),inifile);
        ini_gets("", "address","dummy", dbhost, sizeof(dbhost), inifile);
        ini_gets("", "username","dummy", dbuser, sizeof(dbuser), inifile);
        ini_gets("", "password","dummy", dbpasswd, sizeof(dbpasswd), inifile);
        //connect the database
        if(mysql_init(&mysql) == NULL)
        {
            printf("mysql init error: %s \n",mysql_error(&mysql));
            return -1;
        }
        if(mysql_real_connect(&mysql,
                        dbhost,         //database address
                        dbuser,         //database user
                        dbpasswd,       //database password
                        dbname,         //database name
                         0,
                        NULL,
                         0) == NULL)
        {
            printf("mysql_real_connect error:%s \n",mysql_error(&mysql));
        }
        printf("mysql database connected\n");
    }

    tv.tv_sec = 0;
    tv.tv_usec = 0;
    
    while(tracking)
    {
        FD_ZERO(&rset);
        FD_ZERO(&wset);
        FD_SET(socket_fd,&rset);
        FD_SET(socket_fd,&wset);
    	bzero(recvbuf,sizeof(recvbuf));
        bzero(sendbuf,sizeof(sendbuf));

        socket_fd = tcp_server_accept(listen_fd, ip, &port);
        if(socket_fd > 0)
        {
            printf("ip from client is %s\n",ip);
            strcpy(antinfo.ip, ip);
            printf("port from client is %d\n",port);
            pthread_create(&receive_tid, NULL, (void*) &receive_message, (void*) socket_fd);
            pthread_detach(&receive_tid);
        }
        //json结构



        usleep(loopfreq);
    };
}

//antenna tracking function
int sat_tracking(orbit_t *orb)
{
    int x,y;
    double vel;
    char ant_data[BUFSIZE], buf[BUFSIZE], buffer[BUFSIZE];
    char latitude[10], longtitude[10], altitude[10];
    pthread_t gps_tid, gpsd_tid, ant_tid, send_tid, key_tid;
    WINDOW * track_win;
    void *status;

	bzero(buffer,sizeof(buffer));
    ant_status = satno;     //tracking satellite report

    ini_gets("", "latitude", "dummy", latitude, sizeof(latitude), inifile);
    ini_gets("", "longtitude","dummy", longtitude, sizeof(longtitude), inifile);
    ini_gets("", "altitude", "dummy", altitude, sizeof(altitude), inifile);
    ini_gets("", "loop_freq", "dummy", loop_freq, sizeof(loop_freq), inifile);
    ini_gets("", "gps_time", "dummy", gps_time, sizeof(gps_time), inifile);
    ini_gets("", "antenna_id","dummy", antennaid, sizeof(antennaid), inifile);
    ini_gets("", "gpsd","dummy", gpsd, sizeof(gpsd), inifile);

    loopfreq = atoi(loop_freq)*1000;

    location.lat = atof(latitude);
    location.lng = atof(longtitude);
    location.alt = atof(altitude);
    //mvwprintw(track_win,LOK_Y+6,LOK_X+1,"location lat:%.4f  lng:%.4f   alt:%.2f\n",location.lat,location.lng,location.alt);
    
    if (init_sgdp4(orb) == SGDP4_ERROR)           //check sgp4
    {
        //printf("sgdp4 error\n");
        return -1;  
    }
    //create read_gps or read_gpsd course
    if(atoi(gpsd))
    {
        pthread_create(&gpsd_tid, NULL, (void*) &read_gpsd, NULL);
        pthread_detach(&gpsd_tid);
    }else{
        pthread_create(&gps_tid, NULL, (void*) &read_gps, NULL);
        pthread_detach(&gps_tid);
    }
    //create control_antenna course
    pthread_create(&ant_tid, NULL, (void*) &control_antenna, NULL);
    pthread_detach(&ant_tid);
    //create send_message course
    pthread_create(&send_tid, NULL, (void*) &send_message, NULL);
    pthread_detach(&send_tid);
    //getmaxyx(stdscr,y,x);
    track_win = create_window(50,200,SAT_Y,SAT_X);

    ////////////////////////////////////START OF SAT TRACKING LOOP///////////////////////////////////

    while(tracking)
    {

        if(tcp_error)       //connect error exit stop this thread
        {
            usleep(500000);
            //pthread_join(&send_tid, &status);
            //create send_message course
            pthread_create(&send_tid, NULL, (void*) &send_message, NULL);
            tcp_error = 0;
            continue;
        }
        if(!atoi(gps_time))
            sys_mjd();
        
        //caculate satellite AZ/EL
        sat_angle(orb, mjd, location,&lookangles.range,&lookangles.az,&lookangles.el);

        mvwprintw(track_win,LOK_Y,LOK_X," =========================================================================================\n");
        mvwprintw(track_win,LOK_Y+2,LOK_X,"                                        OpenATS");
        mvwprintw(track_win,LOK_Y+3,LOK_X,"                                    www.openats.cn");
        mvwprintw(track_win,LOK_Y+5,LOK_X," [SYSTEM]");        
        mvwprintw(track_win,LOK_Y+7,LOK_X," UTC TIME:   %d/%d/%d  %.2d:%.2d:%.2d\n",utc->tm_year+1900,utc->tm_mon+1,utc->tm_mday,utc->tm_hour,utc->tm_min,utc->tm_sec);        
        mvwprintw(track_win,LOK_Y+8,LOK_X," SYSTEM MJD TIME:  %f\n",mjd);
        mvwprintw(track_win,LOK_Y+9,LOK_X," SYSTEM LOCATION:");
        mvwprintw(track_win,LOK_Y+10,LOK_X," LAT:%.7f   LON:%.7f   ALT:%.2f\n",location.lat,location.lng,location.alt);
        mvwprintw(track_win,LOK_Y+11,LOK_X,"------------------------------------------------------------------------------------------");        
        mvwprintw(track_win,LOK_Y+13,LOK_X," [GPS]");	    
        mvwprintw(track_win,LOK_Y+15,LOK_X," GPS LAT:%.7f   GPS LNG:%.7f   GPS ALT:%.2f  GPS LOCK:%d\n",gps_info.lat,gps_info.lng,gps_info.alt,gps_info.sat);
        mvwprintw(track_win,LOK_Y+16,LOK_X," GPS TIME:%02d/%02d/%02d  %02d:%02d:%02d    SPEED:%f\n",gps_info.year,gps_info.month,gps_info.day,gps_info.hour,gps_info.min,gps_info.sec,gps_info.speed);
        mvwprintw(track_win,LOK_Y+17,LOK_X,"------------------------------------------------------------------------------------------\n");            
        mvwprintw(track_win,LOK_Y+19,LOK_X," [SATELLITE]");		
        mvwprintw(track_win,LOK_Y+21,LOK_X," SAT:%s        NORAD:%ld       ANTENNA:%s\n\n\n",satname,satno,antennaid);   
        attron(A_BOLD);
        //mvwprintw(track_win,LOK_Y+22,LOK_X,"\n"); 
        mvwprintw(track_win,LOK_Y+22,LOK_X," [ AZ : %0.7f  EL : %0.7f ]\n",lookangles.az,lookangles.el);
        attroff(A_BOLD);
        mvwprintw(track_win,LOK_Y+23,LOK_X,"==========================================================================================\n\n\n\n\n\n\n");
        wrefresh(track_win);
        usleep(loopfreq);
    };
    ////////////////////////////////////END OF SAT TRACKING LOOP///////////////////////////////////
    //pthread_join(ant_tid, &status);
    //pthread_join(send_tid, &status);

    return 0;
}

int target_tracking()
{
    int x,y;
    char buffer[BUFSIZE];
    char latitude[10], longtitude[10], altitude[10];
    char targetlat[10], targetlng[10], targetalt[10];
    pthread_t gps_tid, gpsd_tid, ant_tid, receive_tid, send_tid;
    WINDOW * track_win;
    void *status;

	memset(buffer,0,sizeof(buffer));
    ant_status = 0;     //tracking target report
    ini_gets("", "latitude", "dummy", latitude, sizeof(latitude), inifile);
    ini_gets("", "longtitude","dummy", longtitude, sizeof(longtitude), inifile);
    ini_gets("", "altitude", "dummy", altitude, sizeof(altitude), inifile);
    ini_gets("", "antenna_id","dummy", antennaid, sizeof(antennaid), inifile);
    ini_gets("", "targetlat", "dummy", targetlat, sizeof(targetlat), inifile);
    ini_gets("", "targetlng","dummy", targetlng, sizeof(targetlng), inifile);
    ini_gets("", "targetalt", "dummy", targetalt, sizeof(targetalt), inifile);   
    ini_gets("", "loop_freq", "dummy", loop_freq, sizeof(loop_freq), inifile);   
    ini_gets("", "gpsd","dummy", gpsd, sizeof(gpsd), inifile);
    ini_gets("", "gps_time", "dummy", gps_time, sizeof(gps_time), inifile);

    loopfreq = atoi(loop_freq)*1000;

    //get default antenna location
    location.lat = atof(latitude);
    location.lng = atof(longtitude);
    location.alt = atof(altitude);
    //get default target location
    target.lat = atof(targetlat);
    target.lng = atof(targetlng);
    target.alt = atof(targetalt);
    //create read_gps or read_gpsd course
    if(atoi(gpsd))
    {
        pthread_create(&gpsd_tid, NULL, (void*) &read_gpsd, NULL);
        pthread_detach(&gpsd_tid);
    }else{
        pthread_create(&gps_tid, NULL, (void*) &read_gps, NULL);
        pthread_detach(&gps_tid);
    }
    //create control_antenna course
    pthread_create(&ant_tid, NULL, (void*) &control_antenna, NULL);
    //pthread_detach(&ant_tid);
    //create receive_message course
    pthread_create(&receive_tid, NULL, (void*) &receive_message, NULL);
    pthread_detach(&receive_tid);
    //create send_message course
    pthread_create(&send_tid, NULL, (void*) &send_message, NULL);
    pthread_detach(&send_tid);
    getmaxyx(stdscr,y,x);
    track_win = create_window(40,200,SAT_Y,SAT_X);

    /////////////////////////////START OF TARGET TRACKING LOOP//////////////////////////////////
    //target_tracking loop
    while (tracking)
    {

        if(tcp_error)       //connect error exit stop this thread
        {
            
            usleep(500000);
            //pthread_join(&send_tid, &status);
            pthread_join(&receive_tid, &status);
            //create send_message course
            pthread_create(&send_tid, NULL, (void*) &send_message, NULL);
            pthread_create(&receive_tid, NULL, (void*) &receive_message, NULL);
            tcp_error = 0;
            continue;
        }

        if(!atoi(gps_time))
            sys_mjd();
        
        target_angle(&location,&target);

        mvwprintw(track_win,LOK_Y,LOK_X,"==========================================================================================\n");
        mvwprintw(track_win,LOK_Y+2,LOK_X,"                                        OpenATS");
        mvwprintw(track_win,LOK_Y+3,LOK_X,"                                    www.openats.cn");
        mvwprintw(track_win,LOK_Y+5,LOK_X," [SYSTEM]");  
        mvwprintw(track_win,LOK_Y+7,LOK_X," UTC TIME:   %d/%d/%d  %.2d:%.2d:%.2d\n",utc->tm_year+1900,utc->tm_mon+1,utc->tm_mday,utc->tm_hour,utc->tm_min,utc->tm_sec);
        mvwprintw(track_win,LOK_Y+8,LOK_X," SYSTEM MJD TIME:  %f\n",mjd);
        mvwprintw(track_win,LOK_Y+9,LOK_X," LOCATION LAT:%.7f   LON:%.7f   ALT:%.2f\n",location.lat,location.lng,location.alt);
        mvwprintw(track_win,LOK_Y+10,LOK_X,"------------------------------------------------------------------------------------------");        
        mvwprintw(track_win,LOK_Y+12,LOK_X," [GPS]");
        mvwprintw(track_win,LOK_Y+14,LOK_X," GPS LAT:%.7f   GPS LNG:%.7f   GPS ALT:%.2f  GPS LOCK:%d\n",gps_info.lat,gps_info.lng,gps_info.alt,gps_info.sat);
        mvwprintw(track_win,LOK_Y+15,LOK_X," GPS TIME:%02d/%02d/%02d  %02d:%02d:%02d    SPEED:%f\n",gps_info.year,gps_info.month,gps_info.day,gps_info.hour,gps_info.min,gps_info.sec,gps_info.speed);
        mvwprintw(track_win,LOK_Y+16,LOK_X,"------------------------------------------------------------------------------------------\n");            
        mvwprintw(track_win,LOK_Y+18,LOK_X," [TARGET]");	
        mvwprintw(track_win,LOK_Y+20,LOK_X," TARGET LAT:%.7f    LNG:%.7f    ALT:%.2f      ANTENNA:%s\n\n\n",target.lat,target.lng,target.alt,antennaid);   
        //mvwprintw(track_win,LOK_Y+21,LOK_X,"\n"); 
        attron(A_BOLD);
        mvwprintw(track_win,LOK_Y+21,LOK_X," [ AZ : %0.7f  EL : %0.7f  RANGE : %0.7f]\n",lookangles.az,lookangles.el,lookangles.range);
        attroff(A_BOLD);
		mvwprintw(track_win,LOK_Y+22,LOK_X,"==========================================================================================\n\n\n\n\n\n\n");
        wrefresh(track_win);
        //printf("target.az: %d el:%d\n",lookangles.az, lookangles.el);
        usleep(loopfreq);
    };
    /////////////////////////////END OF TARGET TRACKING LOOP/////////////////////////////////
    pthread_join(ant_tid, &status);
    //pthread_join(send_tid, &status);
    return 0;
}


void target_angle(s_geogcord *ground,s_geogcord *obs)
{
    
    /* WGS 84 Constants */
    double radius_e = 6378137;
    double f = 1/298.257223563;
    double ecc = 8.1819190842621E-2;

    /* Convert from decimal degrees to radians */
    double lat_ground = (M_PI/180)*ground->lat;
    double lon_ground = (M_PI/180)*ground->lng;
    double lat_obs = (M_PI/180)*obs->lat;
    double lon_obs = (M_PI/180)*obs->lng;

    /* Add radius of earth to altitudes */
    double r_ground = radius_e + ground->alt;
    double r_obs = radius_e + obs->alt;

    /* WGS 84 Geoid */
    double N = radius_e / sqrt(1 - pow(ecc,2)*pow(sin(lat_ground),2));

    /* Convert ground to Earth Centered Rotational (ECR) coordinates */
    double x_ground = (N+ground->alt)*cos(lat_ground)*cos(lon_ground);
    double y_ground = (N+ground->alt)*cos(lat_ground)*sin(lon_ground);
    double z_ground = (N*(1 - pow(ecc,2)) + ground->alt)*sin(lat_ground);

    /* Convert obv station to Earth Centered Rotational (ECR) coordinates */
    double x_obs = (N+obs->alt)*cos(lat_obs)*cos(lon_obs);
    double y_obs = (N+obs->alt)*cos(lat_obs)*sin(lon_obs);
    double z_obs = (N*(1 - pow(ecc,2)) + obs->alt)*sin(lat_obs);

    /* Calculate the range vector */
    double range_v_x = x_obs - x_ground;
    double range_v_y = y_obs - y_ground;
    double range_v_z = z_obs - z_ground;

    /* Transform range vector to Topocenteric Horizon */
    double rot_s = sin(lat_ground)*cos(lon_ground)*range_v_x + sin(lat_ground)*sin(lon_ground)*range_v_y - cos(lat_ground)*range_v_z;
    double rot_e = -1*sin(lon_ground)*range_v_x + cos(lon_ground)*range_v_y;
    double rot_z = cos(lat_ground)*cos(lon_ground)*range_v_x + cos(lat_ground)*sin(lon_ground)*range_v_y + sin(lat_ground)*range_v_z;

    double range = sqrt(pow(rot_s, 2) + pow(rot_e, 2) + pow(rot_z, 2));

    /* Calculate elevation and take care of divide by zero if they're the same point */
    double el = 0;
    if(range == 0){
        el = (M_PI)/2;
    }else{
        el = asin(rot_z/range);
    }

    /* Calculate the azmuth and take care of divide by zero */
    double az = 0;
    if(rot_s == 0){
        az = (M_PI)/2;
    }else{
        az = atan(-1*(rot_e/rot_s));
    }

    if(az < 0){
        az = az+(2*M_PI);
    }

    lookangles.az = az*(180/M_PI);
    lookangles.el = el*(180/M_PI);
    lookangles.range = range;
    
    return;
}


/* fix_angle function
/
注意，此函数用于动中通技术，需要采集测向板卡方位角数据或者是传感器数据来对天线运行角度做修正，目前仅修正方位角，后续添加仰角修正
天线方位角修正逻辑：
天线当前方位角：now
卫星方位角：target
角度差：△A
天线运行角度：△B
前提：
如果获取GPS方位角失败（天线是固定站，没有方位角传感器），则天线实际运行角度△B = target（不做修正）,上级函数体现

如果：target ≥ now  则：角度差△A = target - now       
如果：target < now  则：角度差△A = -（now - target）

如果：角度差△A < -180°
则天线运行角度△B = 360° + △A 
如果：角度差△A > 180°
则天线运行角度△B = -（360-△A）
以上都不满足：
则：△B = △A

场景一：
天线当前方位角50°，卫星方位角度为160°
△A = 160-50 = 110°
△B = △A = 110°

场景二：
天线当前方位角290°，卫星方位角30°
△A = -（290-30）= -260°
△B = 360 + （-260） = 100°
*/
void fix_angle (s_azelcord *ant_angles, gpsinfo *gpsdata)
{
    float a;

    if (ant_angles->az < gpsdata->offset)
        a = -(gpsdata->offset - ant_angles->az);
    else a = (ant_angles->az - gpsdata->offset);
    
    if (a < -180) ant_angles->az = (360 + a);
    else if (a > 180) ant_angles->az = -(360 - a);
    else ant_angles->az = a;
    return;
}


void hand_control ()
{
    serial *ant_s;
    float az,el;
    char ant_port[20],ant_baudrate[10],c,input[50];
    char input_az[15],input_el[15];
    char minel[10], maxel[10];

    ini_gets("","ant_port","dummy",ant_port,sizeof(ant_port),inifile);
    ini_gets("","ant_baudrate","dummy",ant_baudrate,sizeof(ant_baudrate),inifile);
    ini_gets("", "maxel", "dummy", maxel, sizeof(maxel), inifile); 

    printf("Ant port is%s\n",ant_port);

    if (serial_open(&ant_s, ant_port, atoi(ant_baudrate)) != 0){
		printf("open antenna port error\n");
        //报错并结束程序
	}

    while(tracking)
    {
        strcpy(input_az,"0");
        strcpy(input_el,"0");
        printf("--------------------------------------------------\n");
        printf("     input 'AZ EL' forexample : 123.49 56.73\n");
        printf("         input 'q' or 'exit' to exit\n");
        printf("        input 'X' 'Y' to fix antenna\n\n\n");
        printf("---> ");
        fgets(input,50,stdin);
        sscanf(input,"%s %s",input_az,input_el);
        if( strncmp(input,"exit",4) == 0 || strncmp(input,"EXIT",4) == 0 ||
            strncmp(input,"q",1) == 0 || strncmp(input,"Q",1) == 0)
        {
            printf("tracking stop\n");
            tracking = 0;
            break;
        }
        else if(strncmp(input,"X",1) == 0 || strncmp(input,"Y",1) == 0 )
        {
            serial_write(ant_s, input);
            printf("\n<FIX CMD>: %s\n",input);
        }
        else if(atof(input_az) > -150 && atof(input_az) < 450 && atof(input_el) > -10 && atof(input_el) < atoi(maxel))
        {
            printf("\n-->[AZ:%s  EL:%s]\n\n\n\n\n",input_az,input_el);
            serial_write(ant_s, input);
        }
        else{
            printf("\n[Are you kidding me?]\n");
        }

    };
    
    strcpy(input, "0");
    if(serial_write(ant_s, input) == 0)  //antenna reset to home
    {
        printf("control antenna go home ok\n");
        return;
    }else{
        printf("can't control antenna go home\n");
    }

    return;
}



/////////////////////////////////STRAT OF MAIN/////////////////////////////////

int main (int argc, char * const *argv)
{
    double az, el;
    int mode, arg, x, y;
    char buf[BUFSIZE],str,str1,str2,str3;
    char address[100],fixangle[1];
    getmaxyx(stdscr,y,x);

/////    ./OpenATS -s satellite     /////
//
// -t : target tracking
// -u : update tle data
// -s : tracking satellite name
// -n : tracking satellite number
// -i : show satellite information
// -g : show GPS module data
// -h : help page
// 
////////////////////////////////////////

    while((arg = getopt(argc,argv,"ts:n:uhi:mrg")) != -1) 
    {
        switch(arg) {
        case 'f':         //fix angle

            break;
        case 't':         //target tracking 
            mode = 1;
            break;
        case 's':         //tracking satellite by satellite name
            if((argv[optind+3]) == NULL)
                sprintf(satname,"%s %s %s %s",optarg,argv[optind],argv[optind+1],argv[optind+2]);
            else if((argv[optind+2]) == NULL)
                sprintf(satname,"%s %s %s",optarg,argv[optind],argv[optind+1]);
            else if((argv[optind+1]) == NULL)
                sprintf(satname,"%s %s",optarg,argv[optind]);
            else if(argv[optind] == NULL)
                sprintf(satname,"%s",optarg);
            mode = 2;
            break;
        case 'n':        //tracing satellite by satellite number
            satno = atoi(optarg);
            mode = 3;
            break;
        case 'u':        // update tle data
            mode = 4;
            break;
        case 'm':        //hand control
            mode = 5;   
            break;
        case 'h':        //help page
            help();
            return 0;
            break;
        case 'g':        //show GPS module data
            read_gps(1);
            break;
        case 'r':        //remote control
            mode = 6;
            break;
        case 'i':        //show satellite info
            satno = atoi(optarg);
            mode = 7;
            break;
        default:        
            help();
            return 0;
            break;
        }
    };

    if(argc == 0)
    {
        help();
        exit(0);
    }

    signal(SIGINT, gohome);
    signal(SIGTERM, gohome);
    signal(SIGQUIT, gohome);
    signal(SIGINT, gohome);
    signal(SIGPIPE,tcperror);

    /*  
    initscr();
    start_color();
    //getmaxyx(stdscr,y,x);
    init_pair(1,COLOR_WHITE,COLOR_BLUE);
    init_pair(2,COLOR_BLUE,COLOR_WHITE);
    init_pair(3,COLOR_RED,COLOR_WHITE);
 
    box(window1,ACS_VLINE,ACS_HLINE);
    box(window2,ACS_VLINE,ACS_HLINE);
    box(window3,ACS_VLINE,ACS_HLINE);
    
    int i = 250;

    attron(A_CHARTEXT);
    //int box(WINDOW *window1,2,2);
    wprintw(window1,"i:%d   X:%d\n",i,i++);
    attroff(A_CHARTEXT);
    wrefresh(window1);
    //refresh();
    attron(COLOR_PAIR(1));
    //int box(WINDOW *window1,2,2);
    wprintw(window2,"i:%d   X:%d\n",i,++i);
    wrefresh(window2);
    attroff(COLOR_PAIR(1));

    attron(COLOR_PAIR(2));
    //int box(WINDOW *window1,2,2);
    wprintw(window2,"i:%d   X:%d\n",i,++i);
    wrefresh(window2);
    attroff(COLOR_PAIR(2));
    //attron(A_BOLD);
    //mvprintw(5,5,"OpenATS");
    //mvprintw(6,5,"www.openats.cn");
    //attroff(A_BOLD);
    refresh();
    */

    ini_gets("", "tlefile", "dummy", tlefile, sizeof(tlefile), inifile);
    ini_gets("", "fixangle", "dummy", fixangle, sizeof(fixangle), inifile);
    fix_az = atoi(fixangle);
    //printf("fix_az is:%d\n",fix_az);

    if(mode == 1)           //target tracking
    {
        initscr();
        start_color();
        refresh();

        if (target_tracking() < 0)
        {
            printf("target_tracking error\n");
            exit(0);
        }
	    endwin();
    }
    else if(mode == 2)         //tracking by satname
	{
        //open TLE file
	    if ((tlefile_fd = fopen(tlefile,"r")) == NULL)     //open with read only
	    {
            printf("can't open TLE file! please check your cmd or file\n");
		    exit(0);
        }
	    //printf("open TLE file success\n");
        if ((satno = search_satno(satname)) < 0)
        {
            printf("can't find the satno by name: %s\n",satname);
            exit(0);
        }
        //printf("SATNAME:%s   NORAD:%ld\n",satname,satno);
        if (read_twoline(tlefile_fd,satno,&orb,satname) < 0)
        {
            printf("can't find the satellite: %s from TLE file\n",satname);
            exit(0);
        }
        //print_orb(&orb);

        initscr();
        start_color();
        refresh();
        
        if (sat_tracking(&orb) < 0)
        {
            printf("sat_tracking error\n");
            exit(0);
        }
	    endwin();
	}
    else if(mode == 3)        //tracking by satno
    {
        //open TLE file
	    if((tlefile_fd = fopen(tlefile,"r")) == NULL)     //open with read only
	    {
            printf("can't open TLE file! please check your cmd or file\n");
		    exit(0);
        }	
	    //printf("open TLE file success\n");
        if(read_twoline(tlefile_fd,satno,&orb,satname) < 0)
        {
            printf("can't find the satellite from TLE file\n");
            exit(0);
        }
		fclose(tlefile_fd);

        initscr();
        start_color();
        refresh();

        //print_orb(&orb);
        if(sat_tracking(&orb) < 0)
        {
            printf("sat_tracking error\n");
            exit(0);
        }
        endwin();
    }
    else if(mode == 4)     //update TLE
    {
        ini_gets("", "update_address", "dummy", address, sizeof(address), inifile);
        printf("address is %s\n",address);
        printf("downloading the TLE file from %s..\n",address);
        printf("please wait...\n");
        system(address);
        exit(0);
    }
    else if(mode == 5)
    {
        hand_control();
        exit(0);
    }
    else if(mode == 6)
    {
        remote_control();
        //send_message();
        exit(0);
    }
    else if(mode == 7)     //print satinfo
    {
	    if((tlefile_fd = fopen(tlefile,"r")) == NULL)     //open with read only
	    {
            printf("can't open TLE file! please check your cmd or file\n");
		    exit(0);
        }	
	    //printf("open TLE file success\n");
        if(read_twoline(tlefile_fd,satno,&orb,satname) < 0)
        {
            printf("can't find the satellite from TLE file\n");
            exit(0);
        }
		fclose(tlefile_fd);
        print_orb(&orb);
    }
    else{
        help();
        exit(0);
    }
    //int sgdp4(double tsince, int withvel, kep_t *kep);
    //void kep2xyz(kep_t *K, xyz_t *pos, xyz_t *vel);
    //int satpos_xyz(double jd, xyz_t *pos, xyz_t *vel)
}
/////////////////////////////////END OF MAIN///////////////////////////////




/* End of OpenATS.c*/