from math import comb, log2
from mpmath import mp, mpf, power, log, nprint

#insert n items into b slot, the maximum bin size is max_b. Now we get the probability.
def max_bin_prob(n, b, max_b):
    i = max_b
    ans = mpf(0.0)
    #decimal precision
    mp.dps = 1000
    while i<=n:
        # ans += (comb(n, i)*(1/b**i)*((1-1/b)**(n-i)))
        tmp = mpf(comb(n, i))
        for j in range(i):
            tmp /= b
        for j in range(n-i):
            tmp *= (1-1/b)
        # ans += comb(n, i)*power(1/b, i)*power(1-1/b, n-i)
        ans += tmp
        i += 1
    ans *= b
    return -log(ans, 2)

def max_bin_prob1(n, b, max_b):
    i = 0
    ans = 0.0
    #decimal precision
    mp.dps = 1000
    
    while i < max_b:
        # ans += (comb(n, i)*(1/b**i)*((1-1/b)**(n-i)))
        # tmp = comb(n, i)
        tmp = mpf(comb(n, i))
        for j in range(i):
            tmp /= b
        for j in range(n-i):
            tmp *= (1-1/b)
        ans += tmp
        nprint(tmp, 20)
        # ans += comb(n, i)*power(1/b, i)*power(1-1/b, n-i)
        i += 1
        # print(ans)
    ans = b*(1-ans)
    return -log2(ans)
# n=3*(1<<24)
# b=int(8192)
# max_b = 7000
n=(1<<10)
b=int(1024)
max_b = 17

p_len = max_bin_prob1(n, b, max_b)
# # print(p_len)
nprint(p_len, 5)
# print(comb(100,10))