#include "gpsdecoder.h"


    /* $GPGGA
     * $GPGGA,hhmmss.ss,llll.ll,a,yyyyy.yy,a,x,xx,x.x,x.x,M,x.x,M,x.x,xxxx*hh
     * WORDS:
     *  0    = $GPGGA - Gloabl Positioning System Fixed Data
     *  1    = UTC of Position
     *  2    = Latitude
     *  3    = N or S
     *  4    = Longitude
     *  5    = E or W
     *  6    = GPS quality indicator (0=invalid; 1=GPS fix; 2=Diff. GPS fix)
     *  7    = Number of satellites in use [not those in view]
     *  8    = Horizontal dilution of position
     *  9    = Antenna altitude above/below mean sea level (geoid)
     *  10   = Meters  (Antenna height unit)
     *  11   = Geoidal separation (Diff. between WGS-84 earth ellipsoid and mean sea level.
     *      -geoid is below WGS-84 ellipsoid)
     *  12   = Meters  (Units of geoidal separation)
     *  13   = Age in seconds since last update from diff. reference station
     *  14   = Diff. reference station ID#
     *  15   = Checksum
     */


    /* $GPRMC
     * note: a siRF chipset will not support magnetic headers.
     * $GPRMC,hhmmss.ss,A,llll.ll,a,yyyyy.yy,a,x.x,x.x,ddmmyy,x.x,a*hh
     * ex: $GPRMC,230558.501,A,4543.8901,N,02112.7219,E,1.50,181.47,230213,,,A*66,
     *
     * WORDS:
     *  0    = $GPRMC - Recommended Minimum Specific GNSS Data
     *  1    = UTC of position fix
     *  2    = Data status (V=navigation receiver warning)
     *  3    = Latitude of fix
     *  4    = N or S
     *  5    = Longitude of fix
     *  6    = E or W
     *  7    = Speed over ground in knots
     *  8    = Track made good in degrees True, Bearing This indicates the direction that the device is currently moving in,
     *       from 0 to 360, measured in “azimuth”.
     *  9    = UT date
     *  10   = Magnetic variation degrees (Easterly var. subtracts from true course)
     *  11   = E or W
     *  12   = Checksum
     */


//get GPS speed and time
int gprmc(gpsinfo *gps, char *packet)
{
    int i;
    char* endptr;
    float time_temp;
    int date_temp;

    i = 7; //dump $GPRMC
    //hhmmss.sss
    //234567.890 
    //printf("packeti %s\n",&packet[i]);
    time_temp = atof(&packet[i])/10000;
    //printf("gprmc time_temp is:%6f\n",time_temp);
    gps->hour = trunc(time_temp);       //23
    //printf("gps hour is %d\n",gps->hour);
    gps->min = trunc((time_temp-gps->hour)*100);    //45
    //printf("gps min is %d\n",gps->min);
    gps->sec = time_temp*10000 - (gps->hour*10000 + gps->min*100);
    //printf("gps sec is :%f\n",gps->sec);

    while(packet[i++] != ',');  //2     
        if(packet[i] == 'A')            //gps data not true
            gps->valid = 1; 
        else{
            gps->valid = 0;
            return -1;
        }
    while(packet[i++] != ',');  //3
    while(packet[i++] != ',');  //4
    while(packet[i++] != ',');  //5
    while(packet[i++] != ',');  //6
    while(packet[i++] != ',');  //7
        gps->speed = atof(&packet[i]);
        gps->speed *= 1.15078;          //covert to mph
    while(packet[i++] != ',');  //8
    while(packet[i++] != ',');  //9     //UTC date
        //ddmmyy
        //234567
        date_temp = atoi(&packet[i]);
        gps->day = trunc(date_temp/10000);
        gps->month = trunc((date_temp/100)-gps->day*100);
        gps->year = date_temp - (gps->day*10000 + gps->month*100);

    return 0;
}


int gpgga(gpsinfo *gps, char *packet)
{
    //$GPGGA,045711.00,3958.89213,N,11629.10691,E,2,08,1.64,110.4,M,-8.2,M,,0000*45
    int i;
    char* endptr;
	double lat_temp,lon_temp;
    double degress,minutes;
    //int c = gpgga_checksum((packet+1));
    i = 7; //dump $GPGGA
    //get UTC time [hhmmss.sss]
    //GpsInfo.PosLLA.TimeOfFix.f = strtod(&packet[i], &endptr);
    //gps->time = strtod(&packet[i], &endptr);
    //printf("gps time is %.0f\n",gps->time);

    while(packet[i++] != ',');              // next field: latitude
     
    // get latitude [ddmm.mmmmm]
    lat_temp = strtod(&packet[i],&endptr)/100;
    degress = trunc(lat_temp);
    minutes = (lat_temp-degress)/0.6;
    //printf("gpgga lat_temp %f\n",lat_temp);

	while(packet[i++] != ',');              // next field: N/S indicator
	if(packet[i] == 'S')
        gps->lat = -(degress + minutes);
	else 
        gps->lat = degress + minutes;
	//printf("gps latitude is:%f\n",gps->lat);

    while(packet[i++] != ',');              // next field: longitude

    // get longitude [ddmm.mmmmm]
    lon_temp = strtod(&packet[i],&endptr)/100;
    degress = trunc(lon_temp);
    minutes = (lon_temp-degress)/0.6;
    //printf("lon_temp is: %f",lon_temp);
    while(packet[i++] != ',');              // next field: E/W indicator
 
    // correct latitute for E/W
	if(packet[i] == 'W')
	 	gps->lng = -(degress + minutes);
	else
        gps->lng = degress + minutes;

	//printf("gps longitude is:%f\n",gps->lng);

    while(packet[i++] != ',');              // next field: position fix status
 
    // position fix status
    // 0 = Invalid, 1 = Valid SPS, 2 = Valid DGPS, 3 = Valid PPS
    // check for good position fix
    if( (packet[i] != '0') && (packet[i] != ',') )
         //GpsInfo.PosLLA.updates++;

    while(packet[i++] != ',');              // next field: satellites used
     
    // get number of satellites used in GPS solution
    gps->sat = atoi(&packet[i]);
    //printf("sat %d\n",gps->sat);

    while(packet[i++] != ',');              // next field: HDOP (horizontal dilution of precision)
    while(packet[i++] != ',');              // next field: altitude
     
    // get altitude (in meters)
    gps->alt = strtod(&packet[i], &endptr);
    //printf("Altitude %.1f\n",gps->alt);
	while(packet[i++] != ',');              // next field: altitude units, always 'M'
    while(packet[i++] != ',');              // next field: geoid seperation
    while(packet[i++] != ',');              // next field: seperation units
    while(packet[i++] != ',');              // next field: DGPS age
    while(packet[i++] != ',');              // next field: DGPS station ID
    while(packet[i++] != '*');              // next field: checksum

    return 0;
}

/* 
int GPS_resolve_GPRMC(gpsinfo *gps, char *packet)  
{  
//$GPRMC,092427.604,V,4002.1531,N,11618.3097,E,0.000,0.00,280812,,E,N*08  
  
      
    if(&packet==',')  
    {  
        ++SectionID;  
        i=0;  
    }  
    else  
    {  
        switch(SectionID)  
        {  
            case 1://02:48:13         
                    gps->time[i++]=data;          
                    if(i==2 || i==5)  
                    {         
                        gps->time[i++]=':';       
                    }                 
                    gps->time[8]='\0';  
                break;  
            case 2:  
                if(data=='A')  
                    GPS_DATA.GPS_sv='>';  
                else  
                    GPS_DATA.GPS_sv='<';  
                break;  
            case 3://3158.4608  
                    GPS_DATA.GPS_wd[++i]=data;    
                    GPS_DATA.GPS_wd[12]='\0';                     
                break;  
                  
            case 4:  
                if(data=='N')  
                    GPS_DATA.GPS_wd[0]='N';  
                else if(data=='S')  
                    GPS_DATA.GPS_wd[0]='S';  
              
                break;  
            case 5://11848.3737,E  
                  
                    GPS_DATA.GPS_jd[++i]=data;    
                    GPS_DATA.GPS_jd[12]='\0';  
                break;  
                  
            case 6:  
                if(data=='E')  
                    GPS_DATA.GPS_jd[0]='E';  
                else if(data=='W')  
                    GPS_DATA.GPS_jd[0]='W';  
                  
                break;  
            case 7://10.05  
                    gps->speed[i++]=data;  
                    gps->speed[4]='\0';                       
                break;  
            case 9://15-07-06 -> 06-07-15  
                    GPS_DATA.GPS_date[i++]=data;      
                    if(i==2 || i==5)                          
                    {  
                        GPS_DATA.GPS_date[i++]='-';  
                    }                                 
                    GPS_DATA.GPS_date[8]='\0';                    
                break;  
        }  
    }         
  
}  

*/
//decode $GPHDT az fix angle
int gphdt(gpsinfo *gps, char *packet)
{
    char fixangle[10];

    sscanf(packet,"$GPHDT,%s,T",fixangle);
    //printf("angle:%s \n",fixangle);
    if(atof(fixangle) != 0)
    {
        printf("fix angle is %.2f\n",atof(fixangle));
        gps->offset = atof(fixangle);
        return 0;
    }else{
        return -1;
    }
}