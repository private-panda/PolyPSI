from math import log2, floor
i=1
j=0
# ell= 3
# alpha = 64
# ell= 2
# alpha = 128
ell= 3
alpha = 8
B = 74
B_prime = B/alpha
print(B_prime)
ans = []
while i < 2**ell:
    j = 0
    while j<= floor(log2(B_prime)/ell):
        # print(i*2**j, end=', ')
        ans.append(i*2**j)
        j += 1
    i += 1
    # print()

print(sorted(ans))
print(set(ans))