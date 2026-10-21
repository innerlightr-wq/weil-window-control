"""Minimal pure-Python linear algebra (float).  Used only for [diagnostic]
exploration; every sign-critical conclusion is re-certified in iv.py arithmetic."""
import math

def zeros(n, m=None):
    m = n if m is None else m
    return [[0.0]*m for _ in range(n)]

def eye(n):
    A = zeros(n); 
    for i in range(n): A[i][i] = 1.0
    return A

def matmul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    C = zeros(n, m)
    for i in range(n):
        Ai = A[i]; Ci = C[i]
        for t in range(k):
            a = Ai[t]
            if a == 0.0: continue
            Bt = B[t]
            for j in range(m): Ci[j] += a*Bt[j]
    return C

def transpose(A): return [list(r) for r in zip(*A)]

def jacobi_eig(A0, tol=1e-14, maxsweep=200):
    """Symmetric eigendecomposition by cyclic Jacobi. Returns (vals, V) with
    A = V diag(vals) V^T."""
    n = len(A0)
    A = [list(r) for r in A0]
    V = eye(n)
    for _ in range(maxsweep):
        off = sum(A[i][j]**2 for i in range(n) for j in range(n) if i != j)
        if off <= tol*tol: break
        for p in range(n-1):
            for q in range(p+1, n):
                if abs(A[p][q]) < 1e-300: continue
                theta = (A[q][q]-A[p][p])/(2*A[p][q])
                t = (1 if theta >= 0 else -1)/(abs(theta)+math.sqrt(theta*theta+1))
                c = 1/math.sqrt(t*t+1); s = t*c
                for k in range(n):
                    akp, akq = A[k][p], A[k][q]
                    A[k][p] = c*akp - s*akq
                    A[k][q] = s*akp + c*akq
                for k in range(n):
                    apk, aqk = A[p][k], A[q][k]
                    A[p][k] = c*apk - s*aqk
                    A[q][k] = s*apk + c*aqk
                for k in range(n):
                    vkp, vkq = V[k][p], V[k][q]
                    V[k][p] = c*vkp - s*vkq
                    V[k][q] = s*vkp + c*vkq
    return [A[i][i] for i in range(n)], V

def psd_sqrt(M):
    vals, V = jacobi_eig(M)
    D = zeros(len(M))
    for i, v in enumerate(vals): D[i][i] = math.sqrt(max(v, 0.0))
    return matmul(matmul(V, D), transpose(V))

def spec_of_DM(Dvec, M):
    """Eigenvalues of diag(Dvec) @ M  (M PSD) = spectrum of M^{1/2} D M^{1/2}."""
    S = psd_sqrt(M)
    n = len(M)
    DS = [[Dvec[i]*S[i][j] for j in range(n)] for i in range(n)]
    X = matmul(S, DS)                      # S D S, symmetric
    for i in range(n):
        for j in range(i+1, n):
            m = 0.5*(X[i][j]+X[j][i]); X[i][j] = X[j][i] = m
    vals, _ = jacobi_eig(X)
    return sorted(vals)

def frob2(A): return sum(x*x for r in A for x in r)
def trace(A):  return sum(A[i][i] for i in range(len(A)))
