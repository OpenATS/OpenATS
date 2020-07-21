#include <stdio.h>
#include "serialib.h"
#include <unistd.h>
#include <pthread.h>
#include <stdlib.h>
#include <inttypes.h>
#include <time.h>
#include "gpsdecoder.h"
#include "socket.h"
#include <pthread.h>

#define BUFFER 128

void readgps(void *vargp);



gpsinfo gps_info = {
	.time = 0.0,
	.lat  = 38.95,
	.lng  = 116.12,
	.alt  = 98,
	.NorS = 'N',
	.EorW = 'E',
	.sat  =  1,
};


int main(void)
{

    int listen_fd,socket_fd,port;	
	char buffer[BUFFER];
	char buf[BUFFER];
	char data[BUFFER], ip[128];
	pthread_t tid;

	//listen_fd = tcp_server("127.0.0.1", "5000");
	pthread_create(&tid, NULL, (void*) &readgps,  NULL);
	//pid_t pid = fork();
	while(1)
	{
		printf(" rasiel lat:%f\n  lon:%f\n  alt:%f\n",gps_info.lat, gps_info.lng, gps_info.alt);
		//sleep(1);
	}
}

void readgps(void *ptr)
{

	serial *s;
	int listen_fd,socket_fd,port;
	char array[10] = "$GPGGA";	
	char buffer[BUFFER];
	char buf[BUFFER];

	if (serial_open(&s, "/dev/ttyUSB0", 9600) != 0){
		printf("Problem with port opening\n");
	}

	memset(buffer,0,sizeof(buffer));

		do{

			serial_read(s, buffer, '\n', 1280);
			if(strncmp(buffer,array,6)==0)
			{
				sleep(1);
				serial_read(s,buffer,'\n',1280);
				printf("$GPGGA%s",buffer);
				sprintf(buf,"$GPGGA%s",buffer);
				//decode(buf);
				//printf("time is:%s \n",gps_info.time);
				//find_gps_str(buf,&gps1);
				gpgga(&gps_info, buf);
				printf("child lat:%f   child lng:%f    child alt:%f\n",gps_info.lat, gps_info.lng, gps_info.alt);
				printf("=============================================================================================\n");
				/*
				printf("----------------------------\n");
				printf("%d %d\n",gps->lat,gps->lon);
				printf("----------------------------\n");
				*/
				/*
				sprintf(data,"%s %s", gps_info.lat, gps_info.lng);

				socket_fd = tcp_server_accept(listen_fd, ip, &port);
				printf("")
				tcp_send(socket_fd, data, sizeof(data));
				tcp_recv(socket_fd, data, sizeof(data));
				printf("recv from client :%s\n",data);
				*/
			}
		
		}while(1);
}
