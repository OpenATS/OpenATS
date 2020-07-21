#ifndef _SOCKET_H
#define _SOCKET_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <netdb.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <fcntl.h>

int inet_generic_ntop(struct sockaddr *src, char *ip, int size);
int inet_generic_port(struct sockaddr *src, int *port);
int inet_resolve(const char *host, const char *service, 
                 struct addrinfo **res, int ai_passive);
int tcp_server_socket(const char *host, const char *service, 
                      struct sockaddr *addr, int *addr_len);
int tcp_client_socket(const char *host, const char *service,
                      struct sockaddr *server_addr, int *addr_len);
int set_socket_reuseaddr(int socket_fd);
int set_socket_nonblocking(int socket_fd);
int tcp_server_bind(int socket_fd, const struct sockaddr *addr,
                    int addr_len);
int tcp_server_listen(int socket_fd, int backlog);
int __tcp_server_accept(int socket_fd, struct sockaddr *addr,
                      int *addr_len);
int tcp_server_accept(int socket_fd, char *ip, int *port);
int tcp_client_connect(int socket_fd, const struct sockaddr *server_addr,
                       int addr_len);
int tcp_client(const char *host, const char *service);
int tcp_server(const char *host, const char *service);
void tcp_close(int socket_fd);
#endif
