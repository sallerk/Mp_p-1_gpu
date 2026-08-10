#!/usr/bin/env python3
"""Exact CPU reference residues for the G1 gate in src/main.cpp.

Computes 3^(2^n) mod (2^p - 1) with exact integer arithmetic -- no floating
point, no FFT -- so it is an independent check on the GPU transform. Reduction
uses Mersenne folding (x & M) + (x >> p) rather than division.
"""
import sys, time

def expexp2(p, n, base=3):
    M = (1 << p) - 1
    x = base
    for _ in range(n):
        x = x * x
        x = (x & M) + (x >> p)
        if x >= M:
            x -= M
    return x

# (p, n) pairs matching the PARTIAL[] table in src/main.cpp.
CASES = [(859433, 1000), (1257787, 1000)]

if __name__ == '__main__':
    for p, n in CASES:
        t = time.time()
        x = expexp2(p, n)
        print("p=%-8d n=%-5d res64=%016x   (%.1fs)"
              % (p, n, x & 0xFFFFFFFFFFFFFFFF, time.time() - t))
        sys.stdout.flush()
