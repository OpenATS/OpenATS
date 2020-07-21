#ifndef _SERIALIB_H
#define _SERIALIB_H

#include <sys/time.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/shm.h>
#include <termios.h>
#include <string.h>
#include <fcntl.h> 
#include <stdio.h>

typedef struct serial_s {
	speed_t baud;
	char* port;
	int fd;
} serial;

int serial_open(serial **s, char* port, int baud);
int serial_read_char(serial *s, char *p);
int serial_read(serial *s, char *buf, char eol, unsigned int len);
int serial_write_char(serial *s, char p);
int serial_write(serial *s, char* str);

typedef struct timeval timer;

void timer_init(timer **t);
unsigned long int timer_elapsed(timer *t);


#endif