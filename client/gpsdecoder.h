#ifndef _GPSDECODER_H
#define _GPSDECODER_H

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>

typedef struct 
{
    double lat;
    double lng;
    double alt;
    int    sat;
    double offset;
    float  speed;
    int   hour, min, day, year, month, sec;
    float   bearing;	    // GPRMC
    float   dilution;       // GPGGA
    int valid;              //GPS data valid

} gpsinfo;

int gpgga(gpsinfo *gps, char *packet);
int gphdt(gpsinfo *gps, char *packet);
int gprmc(gpsinfo *gps, char *packet);

#endif