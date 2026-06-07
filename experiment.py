import csv
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix,diags
from scipy.sparse.linalg import factorized

WIDTH=0.5
DEPTH=0.5
HEIGHTS=[0.5,0.75,1.0,1.25,1.5,2.0]
NX=4
NY=8
NZ=4
EYOUNGS=20000
POISSON=0.2
RHO=400
GRAVITY=5
FPS=60
SECONDS=10.0
EXPLICIT_SUBSTEP=5
IMPLICIT_SUBSTEP=1
DAMPING=5
OUT_DIR=Path("data")

MU=EYOUNGS/(2*(1+POISSON))
LAMBDA=EYOUNGS*POISSON/((1+POISSON)*(1-2*POISSON))

def theory_compression(height):
    return RHO*GRAVITY*height*height/(2*EYOUNGS)

def build_box_mesh(height):
    points=[]
    for j in range(NY+1):
        y=height*j/NY
        for i in range(NX+1):
            x=WIDTH*(i/NX-0.5)
            for k in range(NZ+1):
                z=DEPTH*(k/NZ-0.5)
                points.append((x,y,z))
    def idx(i,j,k):
        return (j*(NX+1)+i)*(NZ+1)+k
    tets=[]
    for j in range(NY):
        for i in range(NX):
            for k in range(NZ):
                v000=idx(i,j,k)
                v001=idx(i,j,k+1)
                v010=idx(i,j+1,k)
                v011=idx(i,j+1,k+1)
                v100=idx(i+1,j,k)
                v101=idx(i+1,j,k+1)
                v110=idx(i+1,j+1,k)
                v111=idx(i+1,j+1,k+1)
                tets.extend([
                    (v000,v001,v011,v111),
                    (v000,v010,v011,v111),
                    (v000,v001,v101,v111),
                    (v000,v100,v101,v111),
                    (v000,v010,v110,v111),
                    (v000,v100,v110,v111),
                ])
    return np.array(points,dtype=np.float64),np.array(tets,dtype=np.int32)

def assemble(points,tets):
    rows=[]
    cols=[]
    vals=[]
    mass=np.zeros(len(points),dtype=np.float64)
    I=np.eye(3,dtype=np.float64)
    for tet in tets:
        ids=list(tet)
        p0,p1,p2,p3=[points[i] for i in ids]
        Rs=np.column_stack((p1-p0,p2-p0,p3-p0))
        volume=abs(np.linalg.det(Rs))/6.0
        B=np.linalg.inv(Rs)
        g1=B[0]
        g2=B[1]
        g3=B[2]
        g0=-(g1+g2+g3)
        g=[g0,g1,g2,g3]
        for a in ids:
            mass[a]+=RHO*volume/4.0
        for a in range(4):
            for b in range(4):
                ga=g[a]
                gb=g[b]
                Kab=volume*(MU*np.dot(ga,gb)*I+MU*np.outer(gb,ga)+LAMBDA*np.outer(ga,gb))
                for ii in range(3):
                    for kk in range(3):
                        rows.append(3*ids[a]+ii)
                        cols.append(3*ids[b]+kk)
                        vals.append(Kab[ii,kk])
    dof=3*len(points)
    K=coo_matrix((vals,(rows,cols)),shape=(dof,dof)).tocsr()
    M=diags(np.repeat(mass,3),format="csr")
    return K,M,mass

def free_dofs(points):
    free=[]
    for i,p in enumerate(points):
        if p[1]>1e-12:
            free.extend([3*i,3*i+1,3*i+2])
    return np.array(free,dtype=np.int32)

def top_nodes(points,height):
    return np.where(np.abs(points[:,1]-height)<1e-12)[0]

def write_csv(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",newline="") as f:
        writer=csv.writer(f)
        writer.writerow(["height","explicit_10s","implicit_10s","theory"])
        writer.writerows(rows)

def compression(u,top):
    return -float(np.mean(u[3*top+1]))

def run_explicit(K,M,mass,free,top):
    h=1.0/(FPS*EXPLICIT_SUBSTEP)
    frames=int(SECONDS*FPS)
    Kf=K[free][:,free]
    mf=np.repeat(mass,3)[free]
    fg=np.zeros_like(mf)
    fg[1::3]=-mf[1::3]*0.0
    for i,dof in enumerate(free):
        if dof%3==1:
            fg[i]=-mf[i]*GRAVITY
    u=np.zeros(len(free),dtype=np.float64)
    v=np.zeros(len(free),dtype=np.float64)
    for _ in range(frames):
        for _ in range(EXPLICIT_SUBSTEP):
            f=-Kf@u+fg-DAMPING*mf*v
            v+=h*f/mf
            u+=h*v
    full_u=np.zeros(3*len(mass),dtype=np.float64)
    full_u[free]=u
    return compression(full_u,top)

def run_implicit(K,M,mass,free,top):
    h=1.0/(FPS*IMPLICIT_SUBSTEP)
    frames=int(SECONDS*FPS)
    Kf=K[free][:,free]
    Mf=M[free][:,free]
    mf=np.repeat(mass,3)[free]
    fg=np.zeros_like(mf)
    for i,dof in enumerate(free):
        if dof%3==1:
            fg[i]=-mf[i]*GRAVITY
    A=(Mf+h*h*Kf+h*DAMPING*diags(mf,format="csr")).tocsc()
    solve=factorized(A)
    u=np.zeros(len(free),dtype=np.float64)
    v=np.zeros(len(free),dtype=np.float64)
    for _ in range(frames):
        for _ in range(IMPLICIT_SUBSTEP):
            f=-Kf@u+fg-DAMPING*mf*v
            rhs=h*(f-h*(Kf@v))
            dv=solve(rhs)
            v+=dv
            u+=h*v
    full_u=np.zeros(3*len(mass),dtype=np.float64)
    full_u[free]=u
    return compression(full_u,top)

def main():
    rows=[]
    for height in HEIGHTS:
        points,tets=build_box_mesh(height)
        K,M,mass=assemble(points,tets)
        free=free_dofs(points)
        top=top_nodes(points,height)
        explicit=run_explicit(K,M,mass,free,top)
        implicit=run_implicit(K,M,mass,free,top)
        theory=theory_compression(height)
        rows.append((height,explicit,implicit,theory))
        print(height,explicit,implicit,theory)
    write_csv(OUT_DIR/"experiment.csv",rows)
    print("wrote",OUT_DIR/"experiment.csv")

if __name__=="__main__":
    main()
