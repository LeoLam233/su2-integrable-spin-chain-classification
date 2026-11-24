#!/usr/bin/env python3
"""Exact supplementary checks for the reconstructed proof.
The uniform argument is in classification.tex; finite spin checks are regressions,
not evidence by extrapolation. No polynomial elimination is used in the proof.
Requires SymPy (tested version printed below).
"""
import sympy as s
from sympy.physics.wigner import wigner_6j
from math import factorial as fac
from itertools import product
from functools import lru_cache
Q=s.Rational
n=s.Symbol('n',integer=True,positive=True)

def z(x): return s.cancel(x)==0

def allz(M): return all(s.simplify(x)==0 for x in M)

# Exact symbolic evaluation of just the diagonal Racah summands in Appendix A.
# fr(b,a,c) = (b+a)!/(b+c)!, a ratio of finite products.
@lru_cache(None)
def fr(b,a,c):
    return s.prod(b+t for t in range(c+1,a+1)) if a>=c else 1/s.prod(b+t for t in range(a+1,c+1))

def diag_terms(k,i):
    vals=[]
    for r in range(max(0,2*i-k),i+1):
        pref=Q((-1)**r*fac(i)*fac(k-i),fac(i-r)**2*fac(r)*fac(k+r-2*i))
        vals.append(s.factor(pref*fr(n,-i,i-k-r)**2*fr(n,-k+i,-2*i+r)*fr(2*n,-k-i,1-i)*fr(3*n,1-k-r,1-k-i)*(2*n-2*i+1)))
    return vals

expect={1:(n-3)*(11*n-20)/(8*(2*n-5)*(2*n-3)),3:-(5*n*n-32*n+42)/(4*(2*n-7)*(2*n-3)),5:-n*(n-1)/(8*(2*n-7)*(2*n-5))}
for i,r in expect.items():
    assert z(sum(diag_terms(5,i))-r)
    print(f'block 5, channel {i}, {len(diag_terms(5,i))} diagonal Racah summands: PASS')
assert z(sum(diag_terms(6,3)).subs(n,6)+Q(1,2))
p1=(3*n-8)*(3*n-5)/(4*(2*n-5)*(2*n-3))
p3=n*(3*n-8)/(6*(2*n-7)*(2*n-3))
p5=5*(n-4)*(3*n-7)/(12*(2*n-7)*(2*n-5))
p=[p1,p3,p5]
assert all(z(pj-(1+2*expect[j])/3) for j,pj in zip([1,3,5],p))
assert z(sum(p)-1)
d=n/(5*n-6)
assert z(3*(p1+d*p3+p5)-2*(1+d)-(n-6)*(3*n-7)/((2*n-7)*(5*n-6)))
assert z(3*(p1+d*p3+d*p5)-2*(1+d)+3*n*(3*n-7)/(2*(2*n-5)*(5*n-6)))
print('sign-line weights and two mean obstructions (uniform symbolic identities): PASS')
# Reflection identity specialization and the lower J=1 obstruction.
pv,X=s.symbols('p X');mu=Q(1,4)+pv*X;v=Q(1,4)+pv*X**2-mu**2
assert z(1-3*mu+3*mu**2-v-Q(1,4)-pv*X*((4*pv-1)*X-1))
a=[1,0,Q(2,7)];p4=[Q(2,5),Q(1,4),Q(7,20)];mu=sum(x*y for x,y in zip(a,p4));v=sum(y*(x-mu)**2 for x,y in zip(a,p4))
assert mu==Q(1,2) and v==Q(5,28) and -2*(Q(1,4)-v)==-Q(1,7)
print('reflection bifurcation and n=4 lower-block obstruction: PASS')
# Formal orbit expansion. Each term is contracted with the Gram matrix; no
# independence is needed. Signed permutations implement site transpositions.
p,x,y,alpha,beta,c=s.symbols('p x y alpha beta c')
G3=lambda a:(1-a)*s.eye(3)+a*s.ones(3)
C=s.Matrix([[0,c,-c],[-c,0,c],[c,-c,0]])
G=G3(alpha).row_join(C).col_join((-C).row_join(G3(beta)))
P=s.diag(p*s.Matrix([[1,0,0],[0,0,1],[0,1,0]]),-p*s.Matrix([[1,0,0],[0,0,1],[0,1,0]]))
R=s.diag(p*s.Matrix([[0,0,1],[0,1,0],[1,0,0]]),-p*s.Matrix([[0,0,1],[0,1,0],[1,0,0]]))
U=s.diag(x,0,0,y,0,0);V=s.diag(0,x,0,0,y,0)
def mono(letters):
    out=s.zeros(6)
    for bits in product([0,1],repeat=3):
        if not any(bits):continue
        left=s.eye(6);m=None
        for label,bit in zip(letters,bits):
            group,defect=(P,U) if label=='A' else (R,V)
            if bit:m=left*defect if m is None else m*G*defect
            elif m is None:left=left*group
            else:m=m*group.T
        out+=m
    return out
D=2*mono('ABA')-2*mono('BAB')-mono('BAA')-mono('AAB')+mono('ABB')+mono('BBA')
D=D.applyfunc(lambda q:s.factor(s.rem(s.Poly(q,p),s.Poly(p**2-1,p)).as_expr()))
au=2*x*(alpha**2*x**2+alpha*p*x+c**2*x*y+1);bu=x*(p*x*(2*alpha-1)-2);eu=-2*c*p*x*y
assert allz(D[:3,:3]-s.Matrix([[au,0,bu],[0,-au,-bu],[bu,-bu,0]]))
assert allz(D[:3,5]-s.Matrix([eu,eu,0]))
T=s.zeros(5,6);T[0,0]=T[1,1]=1;T[0,2]=T[1,2]=-1;T[2,3]=T[3,4]=T[4,5]=1
assert z((T*D.subs(alpha,-Q(1,2))*T.T)[0,4]-eu)
# Check normalization from the first exceptional defect.
t=-n/(4*(2*n-3));defect=2/(1-2*t);eta=-2*defect*(t*t*defect**2-t*defect+1)
assert z(1-defect/2-d) and z(-Q(1,2)-eta/8+d**3/2)
# Final singlet scalar: independent full formal orbit expansion with dependent orbits.
sub={alpha:Q(1,2),beta:Q(1,7),c:1/s.sqrt(7),p:-1,x:Q(3,2)}
g=G.subs(sub); D6=D.subs(sub);Teta=s.zeros(6);Teta[1,1]=1;Teta[0,0]=-1
f=s.Matrix([1,1,1,0,0,0]); e0=s.eye(6)[:,0]; W=s.eye(6)[:,3:5]; zz=e0-s.eye(6)[:,1]-(1/s.sqrt(7))/(1+Q(1,7))*(s.eye(6)[:,3]+s.eye(6)[:,4])
assert allz(W.T*g*zz) and allz(W.T*g*f)
assert s.simplify((f.T*g*(D6+Q(63,16)*Teta)*g*zz)[0]+Q(27,28)*(3*y+7))==0
print('orbit identities, degenerate reduction, and final singlet scalar: PASS')
# Independent exact Wigner implementation. Finite cross-checks only.
@lru_cache(None)
def Fphys(N,J):
    S=Q(N,2); J=Q(J); js=list(range(min(N,int(J+S)),int(abs(J-S))-1,-1))
    return js,s.Matrix([[(-1)**int(3*S+J)*s.sqrt((2*j+1)*(2*l+1))*wigner_6j(S,S,j,S,J,l) for l in js] for j in js])
js,F=Fphys(4,1);assert js==[3,2,1] and list(F.diagonal())==[Q(1,5),Q(1,2),Q(3,10)]
js,F=Fphys(6,3);assert F[3,3]==-Q(1,2) and F[6,6]==Q(1,7) and F[3,6]**2==Q(1,7)
for N in (5,6,7,10,17):
    js,F=Fphys(N,Q(3*N,2)-5)
    for i in (1,3,5):assert s.simplify(F[i,i]-expect[i].subs(n,N))==0
    O=F.extract([1,3,5],[1,3,5]);Proj=(s.eye(3)+2*O)/3
    assert allz(Proj*Proj-Proj) and s.trace(Proj)==1
print('independent Wigner diagonal and sign-projection cross-check: PASS')

def rc(N,J,a,e):
    js,F=Fphys(N,J);assert allz(F*F-s.eye(len(js)))
    A=s.diag(*[a[N-j] for j in js]);E=s.diag(*[e[N-j] for j in js]);B=F*A*F;comm=A*B-B*A
    return allz(comm*(A+B)-(A+B)*comm-F*E*F+E)
a=[0,1,0,Q(1,4),0,1,Q(7,6)];e=[0,-Q(1,2),0,-Q(1,128),0,-Q(1,2),-Q(217,432)]
for J in range(10):assert rc(6,J,a,e),J
print('exceptional density verified on every physical sector J=0,...,9: PASS')
# A few independent checks of the shared arbitrary-spin recoupling formulas.
for N in (3,4,6,9):
    for k in range(1,N):
        m=N-k;_,F=Fphys(N,Q(3*N,2)-k);_,F2=Fphys(N,Q(3*N,2)-k-1)
        t=Q((-1)**k*fac(2*m+1)*fac(N),fac(m)*fac(N+m+1))
        assert F[k,k]==t and F2[k,k]==t*Q(m-k*(N+2*m+1),2*N) and F2[k+1,k+1]==-t*Q(N+m+1,2*(2*m+1))
        norm=[]
        for i in range(k+1):
            j=N-i; a=2*j*(j+1)-N*(N+m+1)
            assert F2[i,k+1]!=0
            if a:norm.append(s.simplify(F2[i,k]/(F2[i,k+1]*a)))
            else:assert F2[i,k]==0
        assert len(set(norm))==1 and norm[0]!=0
print('overlap formulas independently cross-checked on representative blocks: PASS')
print('SymPy',s.__version__,'; ALL CHECKS PASSED')
