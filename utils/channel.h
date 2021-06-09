#pragma once

#include <unistd.h>   //close  
#include <arpa/inet.h>    //close  
#include <sys/types.h>  
#include <sys/socket.h>  
#include <netinet/in.h>  
#include <sys/time.h> //FD_SET, FD_ISSET, FD_ZERO macros 
#include <string>
#include <iostream>
#include <cassert>
// #include <stringstream>

#include "seal/seal.h"

// const uint64_t NET_MAGIC_HELLO = 0x5052495643415453ull; // 'PRIVCATS'
// const uint32_t NET_MAGIC_VECTOR_UINT64 = 0x76756938ul; // 'vui8'
// const uint32_t NET_MAGIC_CIPHERTEXT = 0x63697074ul; // 'cipt'
// const uint32_t NET_MAGIC_VECTOR_CIPHERTEXT = 0x76636970ul; // 'vcip'
// const uint32_t NET_MAGIC_PUBLIC_KEY = 0x706b6579ul; // 'pkey'
// const uint32_t NET_MAGIC_SECRET_KEY = 0x706b65f7ul; // 'pkey'
// const uint32_t NET_MAGIC_RELIN_KEYS = 0x72656c6eul; // 'reln'
// const uint32_t NET_MAGIC_PARMS = 0x72656c7ful; // 'parms'

// #define MAX_BUFFER_SIZE 2097152

// #define MAX_BUFFER_SIZE 3145728

//when the poly_modulus_degree = 16384, the large relinear_key is as large as 8.014M, the ciphertext size is 1.78M.
//Therefore, we use a buffer with size 9M.
// #define MAX_BUFFER_SIZE 9437184
//  #define MAX_BUFFER_SIZE 8388608
#define MAX_BUFFER_SIZE 10485760



//A channel connects two ports
class Channel{
public:
    int socket_fd; //a socket
private:
    std::string src_ip_addr;
    int src_port;
    struct sockaddr_in src_sockaddr;

    std::string dest_ip_addr;
    int dest_port;
    struct sockaddr_in dest_sockaddr;

    int data_sent;
    int data_received;

    // seal::seal_byte buffer[MAX_BUFFER_SIZE];
    seal::seal_byte *buffer=nullptr;


    // std::shared_ptr<seal::SEALContext> context;
    // std::stringstream send_stream;
    // std::stringstream recv_stream;

public:
    Channel(){
        socket_fd = -1;
        data_sent = 0;
        data_received = 0;
        src_ip_addr = "127.0.0.1";
        src_port = 5555;
        dest_ip_addr = "127.0.0.1";
        dest_port = src_port;
        dest_sockaddr.sin_family = AF_INET;
        inet_pton(AF_INET, dest_ip_addr.c_str(), &dest_sockaddr.sin_addr);
        dest_sockaddr.sin_port = htons(dest_port);

        buffer = (seal::seal_byte *)malloc(MAX_BUFFER_SIZE*sizeof(seal::seal_byte));
    }
    ~Channel(){
        if(buffer != nullptr){
            free(buffer);
        }
        close(socket_fd);
    }
    // void cClose(){
    //     free(buffer);
    //     close(socket_fd);
    // }
public:
    //get the data size sent in unit MB, precision(3)
    double getDataSizeSent(){
        return data_sent*1.0/(1<<20);
    };
    //get the data size received
    double getDataSizeReceived(){
        return data_received*1.0/(1<<20);
    };

    void setSrcIpAddress(std::string &ip_addr_){
        src_ip_addr = ip_addr_;
        inet_pton(AF_INET, src_ip_addr.c_str(), &src_sockaddr.sin_addr);
    }
    void setSrcPort(int port_){
        src_port = port_;
        src_sockaddr.sin_port = htons(src_port);
    }

    void setDestIpAddress(std::string &ip_addr_){
        dest_ip_addr = ip_addr_;
        inet_pton(AF_INET, dest_ip_addr.c_str(), &dest_sockaddr.sin_addr);
    }
    void setDestPort(int port_){
        dest_port = port_;
        dest_sockaddr.sin_port = htons(dest_port);
    }
    void setDestSockaddr(struct sockaddr_in &dest_sockaddr_){
        dest_sockaddr = dest_sockaddr_;
    }

    bool init(){
        socket_fd = socket(AF_INET, SOCK_STREAM, 0);
        if(socket_fd < 0){
            std::cerr<<"Fail to get socket"<<std::endl;
            exit(-1);
        }
        return true;
    }

    bool cBind(){
        // int opt = (1024*1024);
        int opt = 1; //SO_RCVBUF 
        if (setsockopt(socket_fd, SOL_SOCKET, SO_REUSEADDR | SO_REUSEPORT, (const unsigned char *)&opt, sizeof(opt))){ 
            std::cerr<<"Set sock opt failed"<<std::endl;
            exit(EXIT_FAILURE); 
        } 
        // if (setsockopt(socket_fd, SOL_SOCKET, SO_SNDBUF, (const char *)&opt, sizeof(opt))){ 
        //     std::cerr<<"Set sock opt failed"<<std::endl;
        //     exit(EXIT_FAILURE); 
        // } 
        src_sockaddr.sin_family = AF_INET;
        src_sockaddr.sin_addr.s_addr = INADDR_ANY;
        // inet_pton(AF_INET, src_ip_addr.c_str(), &src_sockaddr);
        src_sockaddr.sin_port = htons(src_port);
        if(bind(socket_fd, (struct sockaddr *) &src_sockaddr, sizeof(src_sockaddr)) < 0){
            std::cerr<<"Fail to bind the IP address."<<std::endl;
            return false;
        }
        return true;
    }
    
    bool cListen(){
        if(listen(socket_fd, 5) < 0){ //set the maximum size for the backlog queue to 5
            std::cerr<<"Network listen failed"<<std::endl;
            return false;
        }
        return true;
    }

    //accept to get a new socket, must return by reference
    void cAccept(Channel &chl){
        // Channel chl;
        // struct sockaddr_in tmp_sockaddr;
        int dest_sockaddr_len = sizeof(dest_sockaddr);
        chl.socket_fd = accept(socket_fd, (sockaddr *)&dest_sockaddr, (socklen_t *)(&dest_sockaddr_len));
        if(chl.socket_fd < 0){
            std::cerr<<"Accept error"<<std::endl;
            exit(-1);
        }
        chl.setDestSockaddr(dest_sockaddr);
        // std::cout<<"the accept socket in server is: "<<chl.socket_fd<<std::endl;

        // char buffer[100];
        // int recv_len = recv(chl.socket_fd, buffer, sizeof(buffer), 0);
        // std::cout<<"Recv in accept: "<<recv_len<<", "<<buffer<<std::endl;

        // std::cout<<"Accept to get a socket: "<<chl.socket_fd<<std::endl;
        // char ip_addr[32];
        // inet_ntop(AF_INET, &(tmp_sockaddr.sin_addr), ip_addr, INET_ADDRSTRLEN);
        // dest_ip_addr = std::to_string(ip_addr);
        // dest_port = (int)ntohs(tmp_sockaddr.sin_port);
        // return std::move(chl);
    }

    //connect with the dest host
    bool cConnect(){
        if(connect(socket_fd, (sockaddr *)&dest_sockaddr, sizeof(dest_sockaddr)) < 0){
            std::cerr<<"Connection error!"<<std::endl;
            exit(-1);
        }
        // std::cout<<"connection succeed, the current socket is "<<socket_fd<<std::endl;
        return true;
    }

    int cSend(uint8_t *str, const uint32_t &len){
        // std::cout<<"the current socket to send is "<<socket_fd<<std::endl;
        // send(socket_fd, str, len, 0);
        // std::cout<<"value: str="<<str<<std::endl;
        int send_len = send(socket_fd, str, len, 0);
        data_sent += send_len;
        // std::cout<<send_len<<", ";
        return send_len;
    }
    int cRecv(uint8_t *str, const uint32_t &len){
        // std::cout<<"the current socket to receive is "<<socket_fd<<std::endl;
        int recv_len = recv(socket_fd, str, len, MSG_WAITALL);
        data_received += recv_len;
        // std::cout<<recv_len<<", ";
        return recv_len;
    }
    int cSend(uint32_t &num){
        int send_len = send(socket_fd, &num, sizeof(num), 0);
        data_sent += send_len;
        return send_len;
    }
    int cRecv(uint32_t &num){
        int recv_len = recv(socket_fd, &num, sizeof(num), 0);
        data_received += recv_len;
        
        return recv_len;
    }

    int cSend(seal::EncryptionParameters &parms){
        // assert(context);
        //send the size of the ctxt
        uint32_t parms_size = parms.save(buffer, MAX_BUFFER_SIZE);
        cSend(parms_size);

        // char buffer[parms_size];
        // send_stream.read(buffer, parms_size);
        int send_size = send(socket_fd, buffer, parms_size, 0);
        data_sent += send_size;
        // std::cout<<"Send: byte len "<<buffer1<<std::endl;
        assert(parms_size == send_size);

        return parms_size+sizeof(uint32_t);
        // return send_len1+send_len2;

    }

    int cRecv(seal::EncryptionParameters &parms){
        // assert(context);
        //send the size of the ctxt
        uint32_t parms_size = 0;
        cRecv(parms_size);
        
        int recv_size = recv(socket_fd, buffer, parms_size, 0);
        data_received += recv_size;

        auto load_size = parms.load(buffer, MAX_BUFFER_SIZE);

        assert(parms_size == load_size);
        // std::cout<<"Send: byte len "<<buffer1<<std::endl;

        return parms_size+sizeof(uint32_t);
        // return send_len1+send_len2;
    }



    // void cSetSEALContext(std::shared_ptr<seal::SEALContext> context){
    //     this->context = context;
    // }

    int cSend(seal::Serializable<seal::RelinKeys> &relin_keys){
        //send the size of the ctxt
        uint32_t relin_keys_size = relin_keys.save(buffer, MAX_BUFFER_SIZE);
        cSend(relin_keys_size);

        int send_size = send(socket_fd, buffer, relin_keys_size, 0);
        
        // std::cout<<"Send: byte len "<<buffer1<<std::endl;
        // std::cout<<"send relin_keys size: "<<send_size<<std::endl;
        assert(relin_keys_size == send_size);

        int tmp_data_sent = relin_keys_size+sizeof(uint32_t);
        data_sent += tmp_data_sent;
        return tmp_data_sent;
        // return send_len1+send_len2;

    }

    int cRecv(seal::SEALContext &context, seal::RelinKeys &relin_keys){
        //send the size of the ctxt
        // assert(context);

        uint32_t relin_keys_size = 0;
        cRecv(relin_keys_size);
        // assert(relin_keys_size <= MAX_BUFFER_SIZE);
        // std::cout<<"relin_keys size: "<<relin_keys_size<<std::endl;
        
        //the recv buffer may not large enough to receive the bytes
        uint32_t tmp = relin_keys_size;
        seal::seal_byte *p = buffer;
        while(tmp > 0){
            int recv_size = recv(socket_fd, p, tmp, 0);
            tmp -= recv_size;
            p += recv_size;
        }
        // recv_stream.write(buffer, relin_keys_size);
        // std::cout<<"before loading the relin_keys: "<<relin_keys_size<< std::endl;
        int relin_keys_size1 = relin_keys.load(context, buffer, relin_keys_size);

        // std::cout<<"after loading the relin_keys: "<<relin_keys_size1<<std::endl;
        assert(relin_keys_size == relin_keys_size1);
        // std::cout<<"Send: byte len "<<buffer1<<std::endl;

        int tmp_data_received = relin_keys_size+sizeof(uint32_t);
        data_received += tmp_data_received;
        return tmp_data_received;
    }

    int cSend(seal::Ciphertext &ctxt){
        //send the size of the ctxt
        uint32_t ctxt_size = ctxt.save(buffer, MAX_BUFFER_SIZE);
        cSend(ctxt_size);
        int send_size = send(socket_fd, buffer, ctxt_size, 0);

        assert(ctxt_size == send_size);

        int tmp_data_sent = ctxt_size+sizeof(uint32_t);
        data_sent += tmp_data_sent;
        return tmp_data_sent;
    }

    int cRecv(seal::SEALContext &context, seal::Ciphertext &ctxt){
        //send the size of the ctxt
        // assert(context);

        uint32_t ctxt_size = 0;
        cRecv(ctxt_size);
        // assert(ctxt_size <= MAX_BUFFER_SIZE);
        
        //the recv buffer may not large enough to receive the bytes
        uint32_t tmp = ctxt_size;
        seal::seal_byte *p = buffer;
        while(tmp > 0){
            int recv_size = recv(socket_fd, p, tmp, 0);
            tmp -= recv_size;
            p += recv_size;
        }
        // recv_stream.write(buffer, relin_keys_size);
        // std::cout<<"before loading the ctxt: "<<ctxt_size<< std::endl;
        int ctxt_size1 = ctxt.load(context, buffer, ctxt_size);

        // std::cout<<"after loading the ctxt: "<<ctxt_size1<<std::endl;
        assert(ctxt_size == ctxt_size1);
        // std::cout<<"Send: byte len "<<buffer1<<std::endl;

        int tmp_data_received = ctxt_size+sizeof(uint32_t);
        data_received += tmp_data_received;
        return tmp_data_received;
    }

    int cSend(seal::Serializable<seal::Ciphertext> &seri_ctxt){
        //send the size of the serializable ctxt
        uint32_t ctxt_size = seri_ctxt.save(buffer, MAX_BUFFER_SIZE);
        cSend(ctxt_size);
        int send_size = send(socket_fd, buffer, ctxt_size, 0);

        assert(ctxt_size == send_size);

        int tmp_data_sent = ctxt_size+sizeof(uint32_t);
        data_sent += tmp_data_sent;
        return tmp_data_sent;
    }

    // int cRecv(seal::Ciphertext &ctxt){
    //     //send the size of the ctxt
    //     assert(context);

    //     uint32_t ctxt_size = 0;
    //     cRecv(ctxt_size);
    //     assert(ctxt_size <= MAX_BUFFER_SIZE);
        
    //     //the recv buffer may not large enough to receive the bytes
    //     uint32_t tmp = ctxt_size;
    //     seal::SEAL_BYTE *p = buffer;
    //     while(tmp > 0){
    //         int recv_size = recv(socket_fd, p, tmp, 0);
    //         tmp -= recv_size;
    //         p += recv_size;
    //     }
    //     // recv_stream.write(buffer, relin_keys_size);
    //     // std::cout<<"before loading the ctxt: "<<ctxt_size<< std::endl;
    //     int ctxt_size1 = ctxt.load(context, buffer, ctxt_size);

    //     // std::cout<<"after loading the ctxt: "<<ctxt_size1<<std::endl;
    //     assert(ctxt_size == ctxt_size1);
    //     // std::cout<<"Send: byte len "<<buffer1<<std::endl;

    //     int tmp_data_received = ctxt_size+sizeof(uint32_t);
    //     data_received += tmp_data_received;
    //     return tmp_data_received;
    // }
};