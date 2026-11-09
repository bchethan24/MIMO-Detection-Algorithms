import math
import random
from itertools import product

DETECTOR = "MMSE-SIC"
MAX_ML_CANDIDATES = 100000

def constellation(name):
    if name == "QPSK":
        return [complex(1,1)/math.sqrt(2), complex(-1,1)/math.sqrt(2),
                complex(1,-1)/math.sqrt(2), complex(-1,-1)/math.sqrt(2)]
    levels = [-3,-1,1,3]
    return [complex(i,q)/math.sqrt(10) for q in levels for i in levels]

def channel_matrix(nr, nt, channel):
    H = [[0j for _ in range(nt)] for _ in range(nr)]
    if channel == "AWGN":
        for i in range(min(nr,nt)): H[i][i] = 1+0j
    elif channel == "Rayleigh":
        for r in range(nr):
            for c in range(nt):
                H[r][c] = complex(random.gauss(0,1),random.gauss(0,1))/math.sqrt(2)
    else:
        K=5.0
        los=math.sqrt(K/(K+1))/math.sqrt(2)
        scatter=math.sqrt(1/(K+1))/math.sqrt(2)
        for r in range(nr):
            for c in range(nt):
                H[r][c] = los + scatter*complex(random.gauss(0,1),random.gauss(0,1))
    return H

def adjoint(A):
    return [[A[r][c].conjugate() for r in range(len(A))] for c in range(len(A[0]))]

def matmul(A,B):
    return [[sum(A[r][k]*B[k][c] for k in range(len(B)))
             for c in range(len(B[0]))] for r in range(len(A))]

def solve_linear(A,b):
    n=len(A); M=[A[i][:]+[b[i]] for i in range(n)]
    for col in range(n):
        pivot=max(range(col,n),key=lambda r:abs(M[r][col]))
        if abs(M[pivot][col]) < 1e-13:
            raise ValueError("Singular/ill-conditioned matrix; try another channel realization.")
        M[col],M[pivot]=M[pivot],M[col]
        p=M[col][col]
        for j in range(col,n+1): M[col][j] /= p
        for r in range(n):
            if r==col: continue
            f=M[r][col]
            for j in range(col,n+1): M[r][j] -= f*M[col][j]
    return [M[i][n] for i in range(n)]

def linear_weights(H, noise_variance, Es, mmse):
    nr=len(H); nt=len(H[0]); Hh=adjoint(H); G=matmul(Hh,H)
    reg=noise_variance/Es if mmse else 0.0
    for i in range(nt): G[i][i] += reg
    W=[[0j for _ in range(nr)] for _ in range(nt)]
    for c in range(nr):
        sol=solve_linear(G,[Hh[i][c] for i in range(nt)])
        for i in range(nt): W[i][c]=sol[i]
    return W

def nearest(z,S):
    return min(S,key=lambda s:abs(z-s)**2)

def slice_vector(x,S):
    return [nearest(v,S) for v in x]

def linear_detect(y,H,S,nv,Es,mmse):
    W=linear_weights(H,nv,Es,mmse)
    z=[sum(W[i][r]*y[r] for r in range(len(y))) for i in range(len(W))]
    return slice_vector(z,S)

def sic_detect(y,H,S,nv,Es,mmse,adaptive=True):
    nr=len(H); remaining=list(range(len(H[0]))); residual=y[:]
    xhat=[0j]*len(remaining); order=[]
    while remaining:
        Hr=[[H[r][c] for c in remaining] for r in range(nr)]
        W=linear_weights(Hr,nv,Es,mmse)
        z=[sum(W[i][r]*residual[r] for r in range(nr)) for i in range(len(remaining))]
        if adaptive:
            gains=[sum(abs(W[i][r])**2 for r in range(nr)) for i in range(len(remaining))]
            pick=min(range(len(remaining)),key=lambda i:gains[i])
        else:
            norms=[sum(abs(Hr[r][i])**2 for r in range(nr)) for i in range(len(remaining))]
            pick=max(range(len(remaining)),key=lambda i:norms[i])
        stream=remaining[pick]; xhat[stream]=nearest(z[pick],S)
        for r in range(nr): residual[r]-=H[r][stream]*xhat[stream]
        order.append(stream+1); remaining.pop(pick)
    return xhat,order

def ml_detect(y,H,S):
    nt=len(H[0]); candidates=len(S)**nt
    if candidates > MAX_ML_CANDIDATES:
        raise ValueError(f"ML search requires {candidates} candidates; limit is {MAX_ML_CANDIDATES}. Reduce Nt or use QPSK.")
    best=None; best_metric=float("inf")
    for x in product(S,repeat=nt):
        metric=0.0
        for r in range(len(H)):
            err=y[r]-sum(H[r][c]*x[c] for c in range(nt))
            metric += abs(err)**2
        if metric < best_metric:
            best_metric=metric; best=list(x)
    return best

def gs_complex(B):
    m=len(B); n=len(B[0])
    star=[[0j]*n for _ in range(m)]
    mu=[[0j]*n for _ in range(n)]
    norms=[0.0]*n
    for i in range(n):
        v=[B[r][i] for r in range(m)]
        for j in range(i):
            if norms[j] > 1e-14:
                mu[i][j]=sum(B[r][i].conjugate()*star[r][j] for r in range(m))/norms[j]
                for r in range(m): v[r]-=mu[i][j]*star[r][j]
        for r in range(m): star[r][i]=v[r]
        norms[i]=sum(abs(v[r])**2 for r in range(m))
    return mu,norms

def nint(x):
    return math.floor(x+0.5) if x>=0 else math.ceil(x-0.5)

def complex_lll(H,delta):
    m=len(H); n=len(H[0]); B=[row[:] for row in H]
    T=[[complex(1 if r==c else 0,0) for c in range(n)] for r in range(n)]
    if n<2: return B,T,0
    k=1; iterations=0; max_iter=5000
    while k<n and iterations<max_iter:
        iterations+=1; mu,norms=gs_complex(B)
        for j in range(k-1,-1,-1):
            q=complex(nint(mu[k][j].real),nint(mu[k][j].imag))
            if abs(q)>0:
                for r in range(m): B[r][k]-=q*B[r][j]
                for r in range(n): T[r][k]-=q*T[r][j]
                mu,norms=gs_complex(B)
        if norms[k-1]<1e-14: raise ValueError("Channel is not full column rank for LLL.")
        if norms[k] >= (delta-abs(mu[k][k-1])**2)*norms[k-1]:
            k+=1
        else:
            for r in range(m): B[r][k],B[r][k-1]=B[r][k-1],B[r][k]
            for r in range(n): T[r][k],T[r][k-1]=T[r][k-1],T[r][k]
            k=max(k-1,1)
    if iterations>=max_iter: print("Warning: LLL reached iteration limit.")
    return B,T,iterations

def lr_detect(y,H,S,nv,Es,delta):
    Hred,T,it=complex_lll(H,delta)
    W=linear_weights(Hred,nv,Es,True)
    z=[sum(W[i][r]*y[r] for r in range(len(y))) for i in range(len(W))]
    xraw=[sum(T[i][j]*z[j] for j in range(len(z))) for i in range(len(T))]
    return slice_vector(xraw,S),Hred,T,z,it

def mat_cond_2(H):
    # Condition number estimate from Jacobi eigenvalues of real representation of H^H H.
    G=matmul(adjoint(H),H); n=len(G); N=2*n
    A=[[0.0]*N for _ in range(N)]
    for i in range(n):
        for j in range(n):
            a=G[i][j].real; b=G[i][j].imag
            A[i][j]=a; A[i][j+n]=-b; A[i+n][j]=b; A[i+n][j+n]=a
    for _ in range(100*N*N):
        p,q=0,1; largest=0.0
        for i in range(N):
            for j in range(i+1,N):
                if abs(A[i][j])>largest: largest=abs(A[i][j]); p,q=i,j
        if largest<1e-12: break
        app,aqq,apq=A[p][p],A[q][q],A[p][q]
        tau=(aqq-app)/(2*apq)
        t=(1 if tau>=0 else -1)/(abs(tau)+math.sqrt(1+tau*tau))
        c=1/math.sqrt(1+t*t); s=t*c
        for k in range(N):
            if k!=p and k!=q:
                akp,akq=A[k][p],A[k][q]
                A[k][p]=A[p][k]=c*akp-s*akq
                A[k][q]=A[q][k]=s*akp+c*akq
        A[p][p]=app-t*apq; A[q][q]=aqq+t*apq; A[p][q]=A[q][p]=0.0
    eig=[A[i][i] for i in range(N) if A[i][i]>1e-12]
    if len(eig)<N: return float("inf")
    return math.sqrt(max(eig)/min(eig))

def bits(index,k):
    return [(index >> shift)&1 for shift in range(k-1,-1,-1)]

def main():
    print("="*62)
    print(f" {DETECTOR} MIMO DETECTION SIMULATION")
    print("="*62)
    try:
        nt=int(input("Transmit antennas Nt (1-4): "))
        nr=int(input("Receive antennas Nr (1-4): "))
        if not (1<=nt<=4 and 1<=nr<=4): raise ValueError("Nt and Nr must be 1..4.")
        if DETECTOR in ("ZF","ZF-SIC","LR") and nr<nt:
            raise ValueError(f"{DETECTOR} requires Nr >= Nt for this implementation.")
        print("Modulation: 1=QPSK, 2=16-QAM")
        mc=input("Choose modulation: ").strip()
        if mc not in ("1","2"): raise ValueError("Choose 1 or 2.")
        mod="QPSK" if mc=="1" else "16-QAM"; S=constellation(mod); M=len(S); bps=int(math.log2(M))
        print("Channel: 1=Rayleigh, 2=AWGN, 3=Rician")
        cc=input("Choose channel: ").strip()
        channel={"1":"Rayleigh","2":"AWGN","3":"Rician"}.get(cc)
        if channel is None: raise ValueError("Choose channel 1, 2 or 3.")
        snr=float(input("Selected SNR in dB (-5 to 40): "))
        frames=int(input("Frames (1-2000): "))
        symbols=int(input("Symbols per frame (10-1000): "))
        delta=0.75
        if DETECTOR=="LR":
            delta=float(input("LLL delta (0.51 to 0.99): "))
            if not 0.51<=delta<=0.99: raise ValueError("Delta must be 0.51..0.99.")
        if not -5<=snr<=40 or not 1<=frames<=2000 or not 10<=symbols<=1000:
            raise ValueError("Input outside allowed range.")
        random.seed(41)
        snrs=list(range(0,31,5))
        if snr<0 or snr>30: snrs.append(snr)
        Es=sum(abs(s)**2 for s in S)/M
        ber=[]; last=None
        print("\nRunning simulation...")
        for si in snrs:
            nv=10**(-si/10); errors=0; total=0
            for f in range(frames):
                for t in range(symbols):
                    txidx=[random.randrange(M) for _ in range(nt)]
                    x=[S[i] for i in txidx]
                    H=channel_matrix(nr,nt,channel)
                    y=[sum(H[r][c]*x[c] for c in range(nt))+
                       complex(random.gauss(0,math.sqrt(nv/2)),random.gauss(0,math.sqrt(nv/2)))
                       for r in range(nr)]
                    extra={}
                    if DETECTOR=="ML": xhat=ml_detect(y,H,S)
                    elif DETECTOR=="ZF": xhat=linear_detect(y,H,S,nv,Es,False)
                    elif DETECTOR=="MMSE": xhat=linear_detect(y,H,S,nv,Es,True)
                    elif DETECTOR=="ZF-SIC": xhat,order=sic_detect(y,H,S,nv,Es,False,True); extra["order"]=order
                    elif DETECTOR=="MMSE-SIC": xhat,order=sic_detect(y,H,S,nv,Es,True,True); extra["order"]=order
                    else:
                        xhat,Hred,T,z,it=lr_detect(y,H,S,nv,Es,delta)
                        extra.update(Hred=Hred,T=T,z=z,it=it)
                    for k in range(nt):
                        errors += sum(a!=b for a,b in zip(bits(txidx[k],bps),bits(min(range(M),key=lambda i:abs(xhat[k]-S[i])**2),bps)))
                        total += bps
                    if si==snrs[-1] and f==0 and t==0:
                        last=(x[:],xhat[:],H,y[:],extra)
            ber.append(errors/total if total else 0.0)
        print("\nBER versus SNR")
        print(f"{'SNR(dB)':<12}{'BER':<14}")
        for s,b in zip(snrs,ber): print(f"{s:<12}{b:.8e}")
        x,xhat,H,y,extra=last
        def pv(label,v): print(label+": "+"  ".join(f"{a.real:+.4f}{a.imag:+.4f}j" for a in v))
        pv("\nTransmitted x",x); pv("Received y",y); pv("Detected xhat",xhat)
        if "order" in extra: print("Detection order:",extra["order"])
        if DETECTOR=="LR":
            print("LLL iterations:",extra["it"])
            print("cond(H):",mat_cond_2(H),"cond(H*T):",mat_cond_2(extra["Hred"]))
            print("Transformation T:")
            for row in extra["T"]: print("  ",row)
            print("Reduced channel H*T:")
            for row in extra["Hred"]: print("  ",row)
        print("\nDetector:",DETECTOR,"| Modulation:",mod,"| Channel:",channel)
        print(f"Antenna configuration: {nt} x {nr}; selected SNR: {snr} dB")
        print("Simulation completed.")
    except (ValueError, OverflowError) as e:
        print("Input/simulation error:",e)

if __name__=="__main__":
    main()
