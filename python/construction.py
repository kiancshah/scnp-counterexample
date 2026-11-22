"""The explicit bottom u for a top w containing 4231 (the proof of the 'if' direction), tested."""
from core import *
import sys, time

def occurrence(w):
    """A 4231 occurrence i<j<k<l (1-based) with k-j minimal; then no c in (j,k) has w(j)<w(c)<w(k)."""
    n=len(w); best=None
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                for l in range(k+1,n):
                    if w[i]>w[k]>w[j]>w[l]:
                        if best is None or k-j<best[2]-best[1]: best=(i,j,k,l)
    return best

def bottom(w):
    occ=occurrence(w)
    if occ is None: return None
    i,j,k,l=occ; n=len(w)
    A,B,C,D=w[l],w[j],w[k],w[i]
    assert not any(B<w[c]<C for c in range(j+1,k))
    Pset=[p for p in range(n) if A<=w[p]<=D]
    low=[];high=[]
    for p in Pset:
        if p in (i,j,k,l): continue
        if p<j: low.append(p)
        elif p>k: high.append(p)
        elif w[p]<B: low.append(p)
        else:
            assert w[p]>C
            high.append(p)
    u=list(w)
    val=A
    for p in low: u[p]=val; val+=1
    u[i]=val; u[k]=val+1; u[j]=val+2; u[l]=val+3; val+=4
    for p in high: u[p]=val; val+=1
    assert val==D+1
    return tuple(u),(i+1,j+1,k+1,l+1)

def check(w):
    u,(i,j,k,l)=bottom(w)
    assert sorted(u)==list(range(1,len(w)+1))
    def sw(v,a,b):
        x=list(v); x[a-1],x[b-1]=x[b-1],x[a-1]; return tuple(x)
    vL=sw(u,i,j); vR=sw(u,k,l)
    assert (i,j,vL) in covers_up(u) and (k,l,vR) in covers_up(u)
    assert inv(vL)<=inv(w) and inv(vR)<=inv(w)
    assert (j,k) in inv(u) and (j,k) not in inv(w)
    assert bruhat_leq(u,w)
    return u

if __name__=="__main__":
    N=int(sys.argv[1])
    for n in range(4,N+1):
        t=time.time(); cnt=0; bad=0
        for w in permutations(range(1,n+1)):
            if contains(w,(4,2,3,1)) is None: continue
            cnt+=1
            u=check(w)
            if not scnp(u,w): bad+=1
        print(f"n={n}: w containing 4231: {cnt}; explicit u passes all structural checks; non-SCNP confirmed by DP: {bad}  ({time.time()-t:.1f}s)",flush=True)
