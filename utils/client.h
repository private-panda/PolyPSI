#pragma once

#include "channel.h"

class Client{
public:
    Channel channel;
    // std::string dest_ip_address;
    // int dest_port;
public:
    Client(){
        setup();
    }
    void setup(){
        channel.init();
    }
    void setDestIpAddress(std::string ip_addr){
        channel.setDestIpAddress(ip_addr);
    }
    void setDestPort(int port){
        channel.setDestPort(port);
    }
    void connect(){
        channel.cConnect(); 
    }
    void connect(std::string dest_ip_address, int dest_port){
        channel.setDestIpAddress(dest_ip_address);
        channel.setDestPort(dest_port);
        channel.cConnect();
    }
};