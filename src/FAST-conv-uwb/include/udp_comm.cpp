#include "udp_comm.h"

void send_data(std::string & s_addr, std::string & port_num, float data[], const int & data_size)
// maximum Byte number: 65507
{
    int st = socket(AF_INET, SOCK_DGRAM, 0);
    if (st == -1)
    {
        printf("Create socket failed ! error message: %s\n", strerror(errno));
        return;
    }

    struct sockaddr_in addr;
    memset(&addr, 0, sizeof(addr));
    addr.sin_family = AF_INET;
    addr.sin_port = htons(atoi(port_num.c_str()));
    
    if (s_addr == "BROADCAST")
    {
        addr.sin_addr.s_addr = htonl(INADDR_BROADCAST);
        const int opt = 1; // BROADCAST
        if(setsockopt(st, SOL_SOCKET, SO_BROADCAST, (char *)&opt, sizeof(opt)) == -1)
        {
            printf("Set socket error.\n");
            return;
        }
     
    }
    else
    { 
        addr.sin_addr.s_addr = inet_addr(s_addr.c_str());
    }

    if (sendto(st, data, sizeof(float) * data_size, 0, (struct sockaddr *) &addr, sizeof(addr)) == -1)
    {
        printf("Data transfer failed! Error message: %s\n", strerror(errno));
        close(st);
        return;
    }

    close(st);

    return;
}