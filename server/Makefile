CC = gcc
objects=serialib.o socket.o minIni.o deep.o sgdp4.o gpsdecoder.o cJSON.o server.o
server:$(objects)
	gcc -o server $(objects) -lm -lpthread -lncurses -lwebsockets -lgps -I/usr/include/mysql/ $(shell mysql_config --libs) -g 
deep.o:deep.c sgdp4h.h
sgdp4.o:sgdp4h.h
cJSON.o:cJSON.h
gpsdecoder.o:gpsdecoder.h
serialib.o:serialib.h
socket.o:socket.h
server.o:server.c 
minIni.o:minIni.h

.PHONY:clean
clean:
	-rm openats $(objects)