#ifndef UDP_COMM_H
#define UDP_COMM_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <string>
#include <chrono>

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>


void send_data(std::string & s_addr, std::string & port_num, float data[], const int & data_size);

#endif