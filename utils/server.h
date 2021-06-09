#pragma once

#include <vector>
#include "channel.h"

class Server{
private:
    Channel master_channel;
public:
    // std::vector<Channel> channels;
    Channel channel;
    

public:
    Server(){
        setup();
    }
    ~Server(){
    }
    void setup(){
        master_channel.init();
        master_channel.cBind();
        master_channel.cListen();
        // channels.push_back(master_channel.cAccept()); // for every connection with a client, create a communication channel
        // channels.push_back(Channel());
        // master_channel.cAccept(channels[0]);
        // master_channel = channel;
        master_channel.cAccept(channel);
    }
};