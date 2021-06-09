import sympy
import math

# sigma_vbfs = [50, 46, 42, 38]
# ells = [16, 15, 14, 12]

# for i in range(1):
#     for j in range(5000):
#         prime = sympy.ntheory.generate.nextprime(1<<ells[i], j+1)
#         if prime**3 > (1<<sigma_vbfs[i]):
#             print(sigma_vbfs[i], prime, math.log2(prime**3))
#             break


# print(list(sympy.sieve.primerange(1<<16, 1<<16+20)))

#to do batching, a prime in SEAL should be in the form of k*polynomial_modulus*2+1
# the maximum number of ny=2^12, which corresponds to the polynomial modulus 2^14. Then we have plain_modulus=65537. 65537 already meet lambda=40 in VBF.
#  for ny=2^10, we have poly_modulus=2^12. Then we have plain_modulus=40961. However, 40961 cannot meet lambda>=40 when nx=2^24.
#  for ny=2^8, we have poly_modulus=2^10, which needs not to use relinear keys to do key switching. Then we have plain_modulus=12289.
# one more thing, to support nx=2^28, plain_modulus=163841, which is enough for batching for ny=2^12.
poly_modulus = 1<<12
for k in range(100):
    prime = k*poly_modulus*2+1
    if sympy.ntheory.isprime(prime):
        print(k, prime)

# for i in range(10):
#     for k in range(100):
#         prime = k*2**(i+8)*2+1
#         if sympy.ntheory.isprime(prime):
#             print(i, k, prime)