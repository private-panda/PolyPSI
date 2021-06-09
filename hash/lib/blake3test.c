#include "./BLAKE3/c/blake3.h"
// #include "blake3.h"
#include <stdio.h>
#include <unistd.h>

int main() {
  // Initialize the hasher.
  blake3_hasher hasher;
  blake3_hasher_init(&hasher);

  // // Read input bytes from stdin.
  // unsigned char buf[65536];
  // ssize_t n;
  // while ((n = read(STDIN_FILENO, buf, sizeof(buf))) > 0) {
  //   blake3_hasher_update(&hasher, buf, n);
  // }
  blake3_hasher_update(&hasher, "56789", 5);

  // Finalize the hash. BLAKE3_OUT_LEN is the default output length, 32 bytes.
  uint8_t output[BLAKE3_OUT_LEN];
  blake3_hasher_finalize(&hasher, output, BLAKE3_OUT_LEN);


  // Print the hash as hexadecimal.
  for (size_t i = 0; i < BLAKE3_OUT_LEN; i++) {
    printf("%02x", output[i]);
  }
  printf("\n");


  blake3_hasher_update(&hasher, "56789", 5);

  blake3_hasher_finalize(&hasher, output, BLAKE3_OUT_LEN);
  // Print the hash as hexadecimal.
  for (size_t i = 0; i < BLAKE3_OUT_LEN; i++) {
    printf("%02x", output[i]);
  }
  printf("\n");

  printf("finalize: \n");
  blake3_hasher_finalize(&hasher, output, BLAKE3_OUT_LEN);
  // Print the hash as hexadecimal.
  for (size_t i = 0; i < BLAKE3_OUT_LEN; i++) {
    printf("%02x", output[i]);
  }
  printf("\n");

  uint32_t output1=0;
  uint8_t *p = (uint8_t *)&output1;
  blake3_hasher_finalize(&hasher, (uint8_t *)&output1, 4);
  printf("a 32-bit output: %x\n", output1);
  printf("a 32-bit output: %d\n", output1);
  printf("a 32-bit output: %x, %x, %x, %x, %x\n", *p, *(p+1), *(p+2), *(p+3), output1);


  // uint8_t output2[4];
  // blake3_hasher_finalize(&hasher, output2, 4);
  // for (size_t i = 0; i < 4; i++) {
  //   printf("%02x", output2[i]);
  // }
  // printf("\n");

  // uint32_t *pp = (uint32_t *)output2;
  // printf("a 32-bit output: %x\n", *pp);

  return 0;
}
