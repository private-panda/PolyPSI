# we use virtual bloom filter to defactor a long item with length sigma to k short items with length sigma_prime

from mpmath import exp, power, mp, mpf, nstr, log

#m is the length of the hash table, and m=2^(sigma_prime); n is number of items; k is the hash numbers
def func(m, n, k):
    # return power(1-exp(-k*n/m), k)
    return power(1-power(1-1.0/m, k*n), k)

# set n = 2**24+2**16
# 1. k = 2, sigma_prime = 46
# 2. k = 3, sigma_prime = 39
# 3. k = 4, sigma_prime = 37
# 4. k = 5, sigma_prime = 35
# 5. k = 6, sigma_prime = 34

# set n = 2**24+2**12
# 1. k = 2, sigma_prime = 46
# 2. k = 3, sigma_prime = 39
# 3. k = 4, sigma_prime = 36
# 4. k = 5, sigma_prime = 35
# 5. k = 6, sigma_prime = 34

# set n = 2**20+2**12
# 1. k = 2, sigma_prime = 42

# set n = 2**16+2**12
# 1. k = 2, sigma_prime = 38

k = 2
lambda_ = 40
# n = 2**24+2**12
n = 2**28+2**12
# n = 2**16+2**12
# n = 2**12+2**12
# n = 2**12+2**11
# n = 2**12+2**11

sigma_prime = 30

mp.dps = 1000
# p = mpf(1.0/(1<<lambda_))
# # nstr(p, 10)
# # print(p)

# ans = func(1<<sigma_prime, n, k)

# while  p < ans:
#     sigma_prime += 1
#     ans = func(1<<sigma_prime, n, k)
    

# print(sigma_prime)
# nstr(sigma_prime)

ans = func(40961**3, n, k)
# ans = func(163841**3, n, k)

print(log(ans, 2))