import taichi as ti
import numpy as np
from scipy.sparse import coo_matrix,diags
from scipy.sparse.linalg import spsolve
ti.init(arch=ti.cpu)

def load_ball_mesh():
    with open("ball_mesh.txt","r") as f:
        n,m=map(int,f.readline().split())
        p=[list(map(float,f.readline().split())) for _ in range(n)]
        t=[list(map(int,f.readline().split())) for _ in range(m)]
    return np.array(p,dtype=np.float32),np.array(t,dtype=np.int32)

def build_edges(tet):
    e=[]
    for a,b,c,d in tet:
        e.extend([a,b,b,c,c,a,a,d,b,d,c,d])
    return np.array(e,dtype=np.int32)

def build_faces(tet):
    faces={}
    for a,b,c,d in tet:
        for face in ((a,c,b),(a,b,d),(b,c,d),(c,a,d)):
            key=tuple(sorted(face))
            if key in faces:
                faces[key]=None
            else:
                faces[key]=face
    idx=[]
    for face in faces.values():
        if face is not None:
            a,b,c=face
            idx.extend([a,b,c,a,c,b])
    return np.array(idx,dtype=np.int32)

def delta(i,j):
    return 1 if i==j else 0
    
mesh_pos,mesh_tet=load_ball_mesh()
mesh_edge=build_edges(mesh_tet)
mesh_face=build_faces(mesh_tet)
vec3=ti.types.vector(3,dtype=ti.f32)
mat3=ti.types.matrix(3,3,dtype=ti.f32)
EYE3=ti.Matrix([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
RESTITUTION=0.5
NUM_POINTS=mesh_pos.shape[0]
NUM_TETS=mesh_tet.shape[0]
NUM_EDGES=NUM_TETS*6
NUM_FACE_INDICES=mesh_face.shape[0]
EYOUNGS=20000
POISSON=0.2
RHO=400
GRAVITY=5
SUBSTEP=1
MU=EYOUNGS/(2*(1+POISSON))
LAMBDA=EYOUNGS*POISSON/((1+POISSON)*(1-2*POISSON))
GROUND=-3
GROUND_SIZE=3.0
dt=ti.field(dtype=ti.f32,shape=())
origin_pos=ti.Vector.field(3,dtype=ti.f32,shape=(NUM_POINTS))
Rs_inverse=ti.Matrix.field(3,3,dtype=ti.f32,shape=(NUM_TETS))
pos=ti.Vector.field(3,dtype=ti.f32,shape=(NUM_POINTS))
vel=ti.Vector.field(3,dtype=ti.f32,shape=(NUM_POINTS))
tetid=ti.Vector.field(4,dtype=ti.i32,shape=(NUM_TETS))
strain=ti.Matrix.field(3,3,dtype=ti.f32,shape=(NUM_TETS))
stress=ti.Matrix.field(3,3,dtype=ti.f32,shape=(NUM_TETS))
force=ti.Vector.field(3,dtype=ti.f32,shape=(NUM_POINTS))
mass=ti.field(dtype=ti.f32,shape=(NUM_POINTS))
volume=ti.field(dtype=ti.f32,shape=(NUM_TETS))
edgeid=ti.field(dtype=ti.i32,shape=(NUM_EDGES*2))
faceid=ti.field(dtype=ti.i32,shape=(NUM_FACE_INDICES))
ground_pos=ti.Vector.field(3,dtype=ti.f32,shape=(4))
groundid=ti.field(dtype=ti.i32,shape=(12))
K=None
M=None
def build_K():
    rows=[]
    cols=[]
    vals=[]
    B_all=Rs_inverse.to_numpy()
    V_all=volume.to_numpy()
    I=np.eye(3,dtype=np.float32)
    for eid in range(NUM_TETS):
        B=B_all[eid]
        g1=np.array([B[0,0],B[0,1],B[0,2]],dtype=np.float32)
        g2=np.array([B[1,0],B[1,1],B[1,2]],dtype=np.float32)
        g3=np.array([B[2,0],B[2,1],B[2,2]],dtype=np.float32)
        g0=-g1-g2-g3
        g=[g0,g1,g2,g3]
        ids=mesh_tet[eid]
        for a in range(4):
            for b in range(4):
                ga=g[a]
                gb=g[b]
                Kab=V_all[eid]*(MU*np.dot(ga,gb)*I+MU*np.outer(gb,ga)+LAMBDA*np.outer(ga,gb))
                for i in range(3):
                    for k in range(3):
                        rows.append(3*ids[a]+i)
                        cols.append(3*ids[b]+k)
                        vals.append(Kab[i,k])
    return coo_matrix((vals,(rows,cols)),shape=(3*NUM_POINTS,3*NUM_POINTS)).tocsr()
def build_M():
    m=mass.to_numpy()
    diag=np.repeat(m,3)
    return diags(diag,format="csr")

def init_mesh():
    origin_pos.from_numpy(mesh_pos)
    pos.from_numpy(mesh_pos)
    tetid.from_numpy(mesh_tet)
    edgeid.from_numpy(mesh_edge)
    faceid.from_numpy(mesh_face)

@ti.kernel
def init_ground():
    ground_pos[0]=ti.Vector([-GROUND_SIZE,GROUND,-GROUND_SIZE])
    ground_pos[1]=ti.Vector([GROUND_SIZE,GROUND,-GROUND_SIZE])
    ground_pos[2]=ti.Vector([GROUND_SIZE,GROUND,GROUND_SIZE])
    ground_pos[3]=ti.Vector([-GROUND_SIZE,GROUND,GROUND_SIZE])
    groundid[0]=0
    groundid[1]=1
    groundid[2]=2
    groundid[3]=0
    groundid[4]=2
    groundid[5]=3
    groundid[6]=0
    groundid[7]=2
    groundid[8]=1
    groundid[9]=0
    groundid[10]=3
    groundid[11]=2

@ti.kernel
def init():
    dt[None]=1.0/60.0
    for i in range(NUM_POINTS):
        vel[i]=ti.Vector([0.0,0.0,0.0])
        force[i]=ti.Vector([0.0,0.0,0.0])
        mass[i]=0.0
    for i in range(NUM_TETS):
        p0,p1,p2,p3=origin_pos[tetid[i][0]],origin_pos[tetid[i][1]],origin_pos[tetid[i][2]],origin_pos[tetid[i][3]]
        Rs=ti.Matrix.cols([(p1-p0),(p2-p0),(p3-p0)])
        Rs_inverse[i]=Rs.inverse()
        volume[i]=abs(Rs.determinant())/6.0
        tet_mass=volume[i]*RHO
        mass[tetid[i][0]]+=tet_mass/4.0
        mass[tetid[i][1]]+=tet_mass/4.0
        mass[tetid[i][2]]+=tet_mass/4.0
        mass[tetid[i][3]]+=tet_mass/4.0

def solve_force(delta_t:np.float32):
    h=float(delta_t)
    y=pos.to_numpy().reshape(-1)
    x=origin_pos.to_numpy().reshape(-1)
    v=vel.to_numpy().reshape(-1)
    m=mass.to_numpy()
    mass_vec=np.repeat(m,3)
    f=-K@(y-x)
    f[1::3]-=m*GRAVITY
    rhs=h*(f-h*(K@v))
    A=M+h*h*K
    dv=spsolve(A,rhs)
    force_np=(mass_vec*dv/h).reshape(NUM_POINTS,3).astype(np.float32)
    force.from_numpy(force_np)

@ti.kernel
def update_pos(delta_t:ti.f32): # type: ignore
    for i in range(NUM_POINTS):
        if ti.abs(pos[i][1]-GROUND)<1e-5 or pos[i][1]<GROUND:
            vel[i]+=force[i]*delta_t/mass[i]
            vel[i][1]=ti.max(0.0,vel[i][1])
            pos[i]+=vel[i]*delta_t
            pos[i][1]=ti.max(GROUND,pos[i][1])
        else:
            vel[i]+=force[i]*delta_t/mass[i]
            pos[i]+=vel[i]*delta_t
            if ti.abs(pos[i][1]-GROUND)<1e-5 or pos[i][1]<GROUND:
                pos[i][1]=GROUND
                if vel[i][1]<0:
                    vel[i][1]=-RESTITUTION*vel[i][1]
        force[i]=ti.Vector([0.0,0.0,0.0])

def main():
    global K,M
    init_mesh()
    init_ground()
    init()
    K=build_K()
    M=build_M()
    window=ti.ui.Window("implicit tetra",(1024,768),vsync=True)
    canvas=window.get_canvas()
    scene=window.get_scene()
    camera=ti.ui.Camera()
    camera.position(5,0.0,5.0)
    camera.lookat(0.0,0.0,0.0)
    camera.up(0.0,1.0,0.0)
    while window.running:
        camera.track_user_inputs(window,movement_speed=0.03,hold_key=ti.ui.RMB)
        scene.set_camera(camera)
        scene.ambient_light((0.45,0.45,0.45))
        scene.point_light(pos=(2.0,3.0,4.0),color=(1.0,1.0,1.0))
        for _ in range(SUBSTEP):
            solve_force(dt[None]/SUBSTEP)
            update_pos(dt[None]/SUBSTEP)
        scene.mesh(ground_pos,indices=groundid,color=(0.35,0.42,0.38))
        scene.mesh(pos,indices=faceid,color=(0.72,0.35,0.24))
        scene.lines(pos,width=1.0,indices=edgeid,color=(0.2,0.45,0.9))
        scene.particles(pos,radius=0.005,color=(0.9,0.1,0.08))
        canvas.scene(scene)
        window.show()

if __name__=="__main__":
    main()
