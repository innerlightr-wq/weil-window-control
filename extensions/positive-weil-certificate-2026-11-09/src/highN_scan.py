import sys, time
sys.path.insert(0,'/home/elias/GitHub/weil-window-control/src')
import mpmath as mp, zeta_window as zw
for L,name in [((4,5),'4/5'), ((81,100),'81/100')]:
    for N,dps in [(20,80),(24,100),(28,120)]:
        t0=time.time(); mp.mp.dps=dps
        M=zw.build_matrix(mp.mpf(L[0])/L[1], N, dps, sector=zw.EVEN)
        E=mp.eigsy(M, eigvals_only=True); lam=min(E)
        print("L=%-8s N=%-3d dps=%-4d lambda_min = %-22s (%.0fs)"%(name,N,dps,mp.nstr(lam,12),time.time()-t0), flush=True)
