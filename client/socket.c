#include "socket.h"

/**
 * convert numberic ip in struct src to char *
 */
int inet_generic_ntop(struct sockaddr *src, char *ip, int size)
{
    switch(src->sa_family)
    {
        case AF_INET:
        {
            struct sockaddr_in *inet4 = (struct sockaddr_in *)src;
            if(inet_ntop(AF_INET, (void *)&inet4->sin_addr, ip, size) == NULL)
            {
                return -1;
            }
            return 0;
        }
        case AF_INET6:
        {
            struct sockaddr_in6 *inet6 = (struct sockaddr_in6 *)src;
            if(inet_ntop(AF_INET6, (void *)&inet6->sin6_addr, ip, size) == NULL)
            {
                return -1;
            }
            return 0;
        }
        default:
            return -1;
    }
}

/**
 * convert port in struct src to local machine byte order (little-endian or big-endian)
 */
int inet_generic_port(struct sockaddr *src, int *port)
{
    switch(src->sa_family)
    {
        case AF_INET:
        {
            struct sockaddr_in *inet4 = (struct sockaddr_in *)src;
            if(port)
            {
                *port = ntohs(inet4->sin_port);
            }
            return 0;
        }
        case AF_INET6:
        {
            struct sockaddr_in6 *inet6 = (struct sockaddr_in6 *)src;
            if(port)
            {
                *port = ntohs(inet6->sin6_port);
            }
            return 0;
        }
        default:
            return -1;
    }
}

/**
 * create a listen socket 
 */
int tcp_server_socket(const char *host, const char *service,
               struct sockaddr *addr, int *addr_len)
{
    int listen_fd;
    int ret;
    struct addrinfo *res, *res_save;

    ret = inet_resolve(host, service, &res, 1);
    if(ret != 0)
    {
        //log error;
        return -1;
    }

    res_save = res;
    
    assert(res != NULL);

    do
    {
        listen_fd = socket(res->ai_family, res->ai_socktype, res->ai_protocol);
        if(listen_fd != -1)
            break;
        perror("socket create");
    }while(res = res->ai_next);

    if(res == NULL)
        return -1;

    if(addr)
        *addr = *(res->ai_addr);

    if(addr_len)
        *addr_len = res->ai_addrlen;

    freeaddrinfo(res_save);
    
    return listen_fd;
}

/**
 * resolve host and service to addrinfo *
 */
int inet_resolve(const char *host, const char *service, 
                 struct addrinfo **res, int ai_passive)
{
    int ret;
    struct addrinfo hints;
    struct addrinfo *tmp;

    memset(&hints, 0, sizeof(struct addrinfo));
    if(ai_passive)
        hints.ai_flags = AI_PASSIVE; //used for listen fd in server end
    hints.ai_family = AF_UNSPEC;
    hints.ai_socktype = SOCK_STREAM;

    ret = getaddrinfo(host, service, &hints, &tmp);
    if(ret != 0)
    {
        //log error;
        printf("%s\n", gai_strerror(ret));
        return -1;
    }

    *res = tmp;
    return 0;
}

/**
 * create a client socket
 */
int tcp_client_socket(const char *host, const char *service,
                      struct sockaddr *server_addr, int *addr_len)
{
    struct addrinfo *res, *res_save;
    int ret;
    int socket_fd;

    ret = inet_resolve(host, service, &res, 0);
    if(ret == -1)
        return -1;

    res_save = res;

    do
    {
        socket_fd = socket(res->ai_family, res->ai_socktype, res->ai_protocol);
        if(socket_fd != -1)
            break;
    }while(res = res->ai_next);

    if(server_addr)
        *server_addr = *(res->ai_addr);

    if(addr_len)
        *addr_len = res->ai_addrlen;

    freeaddrinfo(res_save);
    
    return socket_fd;
}

/**
 * set socket SO_REUSEADDR
 */
int set_socket_reuseaddr(int socket_fd)
{
    int on = 1;
    
    if(setsockopt(socket_fd, SOL_SOCKET, SO_REUSEADDR, &on, sizeof(on)) == -1)
        return -1;
    return socket_fd;
}

int set_socket_nonblocking(int socket_fd)
{
    if(fcntl(socket_fd, F_SETFL, O_NONBLOCK) == -1)
        return -1;
    return socket_fd;
}

/**
 * bind a socket fd 
 */
int tcp_server_bind(int socket_fd, const struct sockaddr *addr, int addr_len)
{
    if(bind(socket_fd, addr, addr_len) == -1)
        return -1;
    return 0;
}

/**
 * listen a socket fd
 */
int tcp_server_listen(int socket_fd, int backlog)
{
    if(listen(socket_fd, backlog) == -1)
        return -1;
    return 0;
}

/**
 * internal accept
 */
int __tcp_server_accept(int socket_fd, struct sockaddr *addr, 
                      int *addr_len)
{
    int new_socket_fd;

    new_socket_fd = accept(socket_fd, addr, addr_len);
    if(new_socket_fd == -1)
    {
        //log error
        return -1;
    }
    return new_socket_fd;
}

/**
 * public interface for accept
 */
int tcp_server_accept(int socket_fd, char *ip, int *port)
{
    struct sockaddr addr;
    int addr_len = sizeof(struct sockaddr);
    int new_socket_fd;

    new_socket_fd = __tcp_server_accept(socket_fd, &addr, &addr_len);
    if(new_socket_fd == -1)
        return -1;

    if(inet_generic_ntop(&addr, ip, 128) == -1)
        return -1;
    if(port)
    {
        *port = inet_generic_port(&addr, port);
        if(*port == -1)
            return -1;
    }

    return new_socket_fd;
}

/**
 * connect
 */
int tcp_client_connect(int socket_fd, const struct sockaddr *server_addr,
                       int addr_len)
{
    int new_socket_fd;
    
    new_socket_fd = connect(socket_fd, server_addr, addr_len);
    if(new_socket_fd == -1)
    {
        //printf("newsocketfd -1\n");
        //log
        return -1;
    }
    return 0;
}

/**
 * create a socket fd and connect to (host, service)
 */
int tcp_client(const char *host, const char *service)
{
    int socket_fd;
    struct sockaddr addr;
    int addr_len;
    int ret;

    socket_fd = tcp_client_socket(host, service, &addr, &addr_len);
    if(socket_fd == -1)
        return -1;
    ret = tcp_client_connect(socket_fd, &addr, addr_len);;
    if(ret == -1)
        return -1;

    return socket_fd;
}

/**
 * 1.create a AI_PASSIVE socket fd
 * 2.bind the socket
 * 3.listen the socket
 */
int tcp_server(const char *host, const char *service)
{
    struct sockaddr addr;
    int addr_len;
    int listen_fd;
    int ret;

    listen_fd = tcp_server_socket(host, service, &addr, &addr_len);
    if(listen_fd == -1)
        return -1;

    ret = set_socket_reuseaddr(listen_fd);
    if(ret == -1)
        return -1;

    ret = tcp_server_bind(listen_fd, &addr, addr_len);
    if(ret == -1)
        return -1;

    ret = tcp_server_listen(listen_fd, 1024);
    if(ret == -1)
        return -1;

    return listen_fd;
}

/**
 * send size bytes
 */
int tcp_send(int socket_fd, const void *buf, int size)
{
    int count, total_write = 0;
    
    while(total_write != size)
    {
        count = write(socket_fd, buf, size - total_write);
        if(count == 0)
        {
            return total_write;
        }
        if(count == -1)
        {
            return -1;
        }

        total_write += count;
        buf += count;
    }
    
    return total_write;
}

/**
 * recv size bytes
 */
int tcp_recv(int socket_fd, const void *buf, int size)
{
    int count, total_read = 0;
    
    
    while(total_read != size)
    {
        count = read(socket_fd, buf, size - total_read);
        if(count == 0)
        {
            return total_read;
        }
        if(count == -1)
        {
            return -1;
        }

        total_read += count;
        buf += count;
    }

    return total_read;
}

/**
 * close the socket fd
 */
void tcp_close(int socket_fd)
{
    close(socket_fd);
}
