CC = gcc
objects=serialib.o socket.o minIni.o deep.o sgdp4.o gpsdecoder.o OpenATS.o
openats:$(objects)
	gcc -o openats $(objects) -lm -lpthread -lncurses -lgps -lwebsockets -g
deep.o:deep.c sgdp4h.h 
sgdp4.o:sgdp4h.h
gpsdecoder.o:gpsdecoder.h
serialib.o:serialib.h
socket.o:socket.h
minIni.o:minIni.h
OpenATS.o:OpenATS.c

.PHONY:clean
clean:
	-rm openats $(objects)