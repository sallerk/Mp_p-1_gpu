import sys
def primes_upto(n):
    s=[True]*(n+1); s[0]=s[1]=False
    for i in range(2,int(n**0.5)+1):
        if s[i]:
            for j in range(i*i,n+1,i): s[j]=False
    return [i for i,v in enumerate(s) if v]

def stage1_exponent(B1, p):
    E = 2*p
    for q in primes_upto(B1):
        pw=q
        while pw <= B1//q: pw*=q
        E*=pw
    return E

def make_fold(p):
    M = (1<<p) - 1
    def fold(x):
        neg = x < 0
        if neg: x = -x
        while x.bit_length() > p:
            x = (x & M) + (x >> p)
        if x >= M: x -= M
        if neg: x = (M - x) % M
        return x
    return M, fold

def lucasV(seed, n, p):
    M, fold = make_fold(p)
    A, B = seed % M, fold(seed*seed - 2)
    for b in bin(n)[3:]:
        T = fold(A*B - seed)
        if b == '1':
            B = fold(B*B - 2); A = T
        else:
            A = fold(A*A - 2); B = T
    return A

# sanity: ladder vs plain recurrence, small modulus
def slow(seed, n, M):
    a, b = 2, seed
    for _ in range(n-1): a, b = b, (seed*b - a) % M
    return b % M
def ladder_mod(seed, n, M):
    A, B = seed % M, (seed*seed-2) % M
    for b in bin(n)[3:]:
        T = (A*B - seed) % M
        if b=='1': B=(B*B-2)%M; A=T
        else:      A=(A*A-2)%M; B=T
    return A
for n in (1,2,3,10,97,1000,4095):
    assert ladder_mod(3,n,10**9+7) == slow(3,n,10**9+7), n
print("ladder == plain recurrence for n in 1..4095  OK")

for p, B1, seed in [(859433, 100, 3), (859433, 500, 5)]:
    E = stage1_exponent(B1, p)
    v = lucasV(seed, E, p)
    print("p=%-8d B1=%-5d seed=%d  E=%4d bits  V_E res64=%016x"
          % (p, B1, seed, E.bit_length(), v & 0xFFFFFFFFFFFFFFFF))
    sys.stdout.flush()
