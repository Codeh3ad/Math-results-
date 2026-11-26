"""Check the finite incidence and interpolation-position claims without dependencies."""
from itertools import combinations,product
from collections import defaultdict
def check(s,Q,L,t):
    D=next(d for d in range(2,s) if pow(d,(s-1)//2,s)==s-1)
    F=list(product(range(s),repeat=2));one=(1,0)
    def mul(u,v):return ((u[0]*v[0]+D*u[1]*v[1])%s,(u[0]*v[1]+u[1]*v[0])%s)
    def sub(u,v):return ((u[0]-v[0])%s,(u[1]-v[1])%s)
    def inv(u):
        a,b=u;m=pow((a*a-D*b*b)%s,-1,s);return (a*m%s,-b*m%s)
    U=F[1:1+L];I=[inv(u) for u in U]
    assert L<=2*Q and 2*s*t<=L
    # lambda extracts the second coordinate; all nonzero cosets are affine lines avoiding 0.
    colours=[defaultdict(set) for _ in range(s)]
    for x,z in combinations(range(L),2):
        c=(-mul(U[x],U[z])[1])%s
        colours[c][x].add(z);colours[c][z].add(x)
    assert all(len(g[x]&g[z])<=1 for g in colours[1:] for x,z in combinations(range(L),2))
    assert all(len(colours[0][x])<=s-1 for x in range(L))
    directions=[(1,b) for b in range(s)]+[(0,1)]
    def lines(a):
        out=[[] for _ in range(s)]
        for w in F:out[(a[0]*w[0]+a[1]*w[1])%s].append(w)
        return out
    groups=[]
    def group(tag,left,right=None):
        h={}
        for side,pair in enumerate([left]+([right] if right else [])):
            x,C=pair
            for m,w in enumerate(sorted(C)):
                for i in range(t):h[x,w,i]=(side*s+m)*t+i
        assert len(h)<=L and len(set(h.values()))==len(h)
        groups.append((tag,h))
    for d,a in enumerate(directions):
        for x in range(L):
            for C in lines(a):group(('f',d),(x,C))
    joint=defaultdict(set);taken=[set() for _ in U]
    for x,z in combinations(range(L),2):
        collision=x%Q==z%Q
        coeff=inv(sub(I[z],I[x]))
        if collision:
            cells=lines(directions[0])
            for i,C in enumerate(cells):
                for j,E in enumerate(cells):
                    group(('m',(i+j)%s),(x,C),(z,E))
                    joint[x,z].update(product(C,E))
        elif z in colours[0][x]:
            colour=next(v for v in range(2*s) if v not in (taken[x]|taken[z]))
            taken[x].add(colour);taken[z].add(colour)
            cells=[[] for _ in range(s)]
            for w in F:cells[mul(w,coeff)[1]].append(w)
            for C in cells:
                group(('z',colour),(x,C),(z,C));joint[x,z].update(product(C,C))
        for w,v in product(F,repeat=2):
            a=mul(sub(v,w),coeff)
            c=(a[1]+mul(U[x],U[z])[1])%s
            assert (((w,v) in joint[x,z])==(collision or (a[1]==0 and c==0)))
    # Check all overlapping groups. Shared positions differ by a single translation.
    incidence=defaultdict(list);overlaps=0
    for j,(tag,h) in enumerate(groups):
        seen={i for col in h for i in incidence[col]}
        for i in seen:
            prev,old=groups[i]
            common=h.keys()&old.keys()
            assert prev!=tag and len({h[c]-old[c] for c in common})==1
            overlaps+=1
        for c in h:incidence[c].append(j)
    assert max(map(len,incidence.values()))<=3*s+1
    print(s,Q,L,t,len(groups),overlaps,max(map(len,incidence.values())))
if __name__=='__main__':
    print('s Q L t groups overlapping_pairs max_incidence')
    for v in [(3,5,8,1),(5,13,24,2)]:check(*v)
