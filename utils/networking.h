#pragma once

#include <cstdint>
#include <vector>

#include "boost/asio.hpp"
#include "seal/seal.h"

using namespace std;
using namespace boost::asio;

class Networking
{
public:
    Networking(ip::tcp::socket &socket);

    void set_seal_context(shared_ptr<seal::SEALContext> new_context);

    uint32_t read_uint32();
    void write_uint32(uint32_t value);

    uint64_t read_uint64();
    void write_uint64(uint64_t value);

    void read_hello();
    void write_hello();

    void read_uint64s(vector<uint64_t> &values);
    void write_uint64s(vector<uint64_t> &values);

    void read_HEParms(seal::EncryptionParameters &parms);
    void write_HEParms(seal::EncryptionParameters &parms);

    void read_ciphertext(seal::Ciphertext &ciphertext);
    void write_ciphertext(seal::Ciphertext &ciphertext);

    void read_ciphertexts(vector<seal::Ciphertext> &ciphertexts);
    void write_ciphertexts(vector<seal::Ciphertext> &ciphertexts);

    void read_public_key(seal::PublicKey &public_key);
    void write_public_key(seal::PublicKey &public_key);

    void read_secret_key(seal::SecretKey &secret_key);
    void write_secret_key(seal::SecretKey &secret_key);

    void read_relin_keys(seal::RelinKeys &relin_keys);
    void write_relin_keys(seal::Serializable<seal::RelinKeys> &relin_keys);

private:
    ip::tcp::socket &socket;
    boost::asio::streambuf read_buffer;
    std::istream read_stream;
    boost::asio::streambuf write_buffer;
    std::ostream write_stream;

    shared_ptr<seal::SEALContext> seal_context;
};