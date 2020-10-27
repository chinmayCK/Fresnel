import numpy as np
import math, sys 
from scipy.linalg import null_space
from scipy.optimize import minimize

pi=math.pi
theta=pi/3  # angle of incidence with the normal to the surface
phi=pi/3    # angle made by the in-plane wavevector with x-axis of the geometry

## define the material parameters here 
exx=2.+0.1j
exy=0.j
ep=np.array([[exx,exy,0],[-exy,exx,0],[0,0,exx]])
mu=(1.+1e-6*1j)*np.eye(3)
#xi=np.zeros((3,3))
xi=0.1j*np.eye(3)
zeta=-xi
#zeta=np.zeros((3,3))
MM=np.block([[ep,xi],[zeta,mu]])
print(MM)

# check the passivity constraint
Mc=-1j*MM
Mc=Mc+Mc.conj().T

def isPSD(A, tol=1e-8):
  E = np.linalg.eigvalsh(A)
  return np.all(E > tol)

if isPSD(Mc):
    print('Passive medium')
else:
    print('Not a passive medium')
    sys.exit()

kp=math.sin(theta)
kx=kp*math.cos(phi)
ky=kp*math.sin(phi)

#print(kx,ky,kp)

epx=math.sqrt(MM[0,0].real)
#print(epx)

# search for solution of kz around randomly distributed initial points
Nk=5
km=max(epx,1.5*kp)
kzr=np.random.uniform(-km,km,Nk)
kzi=np.random.uniform(-km,km,Nk)

def matcond(x0,MM,kx,ky):
    kz=x0[0]+1j*x0[1]
    kd=np.array([[0,-kz,ky],[kz,0,-kx],[-ky,kx,0]])
    MMk=np.block([[np.zeros((3,3)),kd],[-kd,np.zeros((3,3))]])
    M=MM+MMk
    C=-1*abs(np.linalg.cond(M))
    return C

#x0=[kzr[2],kzi[5]]
#fun= lambda x: matcond(x,MM,kx,ky)
#res=minimize(fun,x0,method='Nelder-Mead')
#print(res.fun)
#print(res.x[0],res.x[1])

def compare(A,z):
    #compare given complex number z with previously obtained numbers in A
    h=abs(A-z)/max(np.abs(A))
    #if (h[h<1e-6].size>0):
    if (len(h[h<1e-6])>0):
        y=0  # z is a solution previously obtained
    else:
        y=1  # z is not obtained previously
    return y

m1=0
ks=[]        
for j in range(0,Nk):
    for l in range(0,Nk):
        #print(j,l,m1)
        #print(ks)
        x0=[kzr[j],kzi[l]]
        fun= lambda x: matcond(x,MM,kx,ky)
        res=minimize(fun,x0,method='Nelder-Mead')
        z=res.x[0]+1j*res.x[1]
        #print(z)
        if (res.fun<=-1e15) and (abs(z)<100) and (res.x[1]<=0):
            # condition number singularity check
            # reasonable solution check
            # Im(kz)<0 for expo. decaying solutions for z<0
            # following part is to omit repeated solutions
            if (m1>0):
              if compare(ks,z)==1:  # if the solution is not obtained prev.
                ks.append(z)
                m1=m1+1
            else:
                ks.append(z)
                m1=m1+1


kz=sorted(ks, key=lambda x: x.real)   # sorting by real part 
print(kz)            

if len(kz)==0:
  print('No solution of kz found inside the material')
  rss=rps=rsp=rpp=0
  print(rss,rps,rsp,rpp)
  sys.exit()
  
def findfields(kx,ky,kz,MM):
  kd=np.array([[0,-kz,ky],[kz,0,-kx],[-ky,kx,0]])
  MMk=np.block([[np.zeros((3,3)),kd],[-kd,np.zeros((3,3))]])
  M=MM+MMk
  E=null_space(M);
  return E

# find the null-space eignevectors inside the medium
# there are always going to be 2 solutions
E=np.zeros((6,2),dtype=complex)
p=0
for j in range(0,len(kz)):
    f=findfields(kx,ky,kz[j],MM)
    for l in range(0,np.size(f,1)):
      for m in range(0,np.size(f,0)):
        E[m][p]=f[m,l]
      p=p+1  
        #E.append(f[:,l])

print(E[:,0],E[:,1])
if np.all((E[:,0]==0)) or np.all((E[:,1]==0)):
  print('Two eigenstates are not found inside the material')
  rss=rps=rsp=rpp=0
  print(rss,rps,rsp,rpp)
  sys.exit()
  
# sp polarization modes in vacuum
def findspvectors(kp,kz,phi,MM):
  kx=kp*math.cos(phi)
  ky=kp*math.sin(phi)
  kd=np.array([[0,-kz,ky],[kz,0,-kx],[-ky,kx,0]])
  ep=MM[0:3,0:3]
  xi=MM[0:3,3:6]
  zeta=MM[3:6,0:3]
  mu=MM[3:6,3:6]
  # s or TE polarization
  ev=([[math.sin(phi)],[-math.cos(phi)],[0]])
  hv=np.linalg.solve(mu,np.matmul(kd-zeta,ev))
  esm=np.vstack((ev,hv))
  esm=esm/np.linalg.norm(esm)
  # p or TM polarization 
  hv=([[math.sin(phi)],[-math.cos(phi)],[0]])
  ev=np.linalg.solve(ep,np.matmul(-(xi+kd),hv))
  epm=np.vstack((ev,hv))
  epm=epm/np.linalg.norm(epm)
  return esm, epm

MMv=np.eye(6)
kz0=np.sqrt(1-kp**2)
esi, epi=findspvectors(kp,-kz0,phi,MMv)
esr, epr=findspvectors(kp,kz0,phi,MMv)

# reflection of s-polarized light
As=np.array([[esr[0],epr[0],-E[0,0],-E[0,1]],
    [esr[1],epr[1],-E[1,0],-E[1,1]],
    [esr[3],epr[3],-E[3,0],-E[3,1]],
             [esr[4],epr[4],-E[4,0],-E[4,1]]],dtype=complex)
Bs=np.array([esi[0],esi[1],esi[3],esi[4]],dtype=complex)
rs=np.linalg.solve(As,-Bs)
rps=rs[1]
rss=rs[0]
print(rss,rps)

# reflection of p-polarized light 
Bp=np.array([epi[0],epi[1],epi[3],epi[4]],dtype=complex)
rp=np.linalg.solve(As,-Bp)
rsp=rp[0]
rpp=rp[1]
print(rsp,rpp)
