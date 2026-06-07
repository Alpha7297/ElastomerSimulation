#include "geometry.hpp"
#include<iostream>
#include<algorithm>
#include<cstdio>
#include<cstdlib>
#include<fstream>
#include<map>
#include<set>
#include<string>

double r_avg=0.05;
int k_poss=30;
using std::vector;

typedef struct RES3D{
    vector<point> pr;
    vector<tetrahedron> tetr;
}RES3D;

typedef struct face{
    int id1,id2,id3;
    bool operator<(const face& other)const{
        if(id1!=other.id1)return id1<other.id1;
        if(id2!=other.id2)return id2<other.id2;
        return id3<other.id3;
    }
}face;

vector<int> gen_active;
vector<point> gen_p;
vector<tetrahedron> gen_tet;

face make_face(int a,int b,int c){
    if(a>b)std::swap(a,b);
    if(a>c)std::swap(a,c);
    if(b>c)std::swap(b,c);
    return {a,b,c};
}

void add_tet(int ia,int ib,int ic,int id){
    point a=gen_p[ia];
    point b=gen_p[ib];
    point c=gen_p[ic];
    point d=gen_p[id];
    tetrahedron t={a,b,c,d,ia,ib,ic,id};
    if(t.volume>EPS){
        gen_tet.push_back(t);
    }
}

int in_sphere(const tetrahedron& t,const point& p){
    if(t.r<EPS)return 0;
    point pp=p;
    point oo=t.o;
    return distance(pp,oo)<t.r-EPS;
}

void addpoint(point neww){
    neww.edge=0;
    gen_p.push_back(neww);
    int new_idx=gen_p.size()-1;
    vector<int> to_delete;
    std::map<face,int> faces;
    for(int i=0;i<gen_tet.size();i++){
        if(in_sphere(gen_tet[i],gen_p[new_idx])){
            to_delete.push_back(i);
            faces[make_face(gen_tet[i].ia,gen_tet[i].ib,gen_tet[i].ic)]++;
            faces[make_face(gen_tet[i].ia,gen_tet[i].ib,gen_tet[i].id)]++;
            faces[make_face(gen_tet[i].ia,gen_tet[i].ic,gen_tet[i].id)]++;
            faces[make_face(gen_tet[i].ib,gen_tet[i].ic,gen_tet[i].id)]++;
        }
    }
    vector<tetrahedron> new_tet;
    int id2=0;
    for(int i=0;i<gen_tet.size();i++){
        if(id2<to_delete.size()&&to_delete[id2]==i){
            id2++;
            continue;
        }
        new_tet.push_back(gen_tet[i]);
    }
    gen_tet.swap(new_tet);
    for(const auto& pair:faces){
        if(pair.second==1){
            const face& f=pair.first;
            add_tet(f.id1,f.id2,f.id3,new_idx);
        }
    }
    gen_active.push_back(new_idx);
}

void del_super_tets(int super_cnt){
    vector<tetrahedron> new_tet;
    for(int i=0;i<gen_tet.size();i++){
        if(gen_tet[i].ia<super_cnt||gen_tet[i].ib<super_cnt||gen_tet[i].ic<super_cnt||gen_tet[i].id<super_cnt){
            continue;
        }
        gen_tet[i].ia-=super_cnt;
        gen_tet[i].ib-=super_cnt;
        gen_tet[i].ic-=super_cnt;
        gen_tet[i].id-=super_cnt;
        new_tet.push_back(gen_tet[i]);
    }
    gen_tet.swap(new_tet);
    gen_p.erase(gen_p.begin(),gen_p.begin()+super_cnt);
}

RES3D tetrahedron_gen(vector<point> samples){
    for(int i=0;i<samples.size();i++){
        samples[i].x+=(rand()%100)*1e-10;
        samples[i].y+=(rand()%100)*1e-10;
        samples[i].z+=(rand()%100)*1e-10;
    }
    gen_p.clear();
    gen_tet.clear();
    gen_active.clear();
    if(samples.empty()){
        return {gen_p,gen_tet};
    }
    double minx=samples[0].x,maxx=samples[0].x;
    double miny=samples[0].y,maxy=samples[0].y;
    double minz=samples[0].z,maxz=samples[0].z;
    for(int i=1;i<samples.size();i++){
        minx=std::min(minx,samples[i].x);
        maxx=std::max(maxx,samples[i].x);
        miny=std::min(miny,samples[i].y);
        maxy=std::max(maxy,samples[i].y);
        minz=std::min(minz,samples[i].z);
        maxz=std::max(maxz,samples[i].z);
    }
    double cx=(minx+maxx)/2.0;
    double cy=(miny+maxy)/2.0;
    double cz=(minz+maxz)/2.0;
    double span=std::max(maxx-minx,std::max(maxy-miny,maxz-minz));
    if(span<EPS)span=1.0;
    double R=span*20.0+1.0;
    gen_p.push_back({cx,cy+3.0*R,cz});
    gen_p.push_back({cx-2.0*R,cy-R,cz-2.0*R});
    gen_p.push_back({cx+2.0*R,cy-R,cz-2.0*R});
    gen_p.push_back({cx,cy-R,cz+2.0*R});
    add_tet(0,1,2,3);
    for(int i=0;i<samples.size();i++){
        addpoint(samples[i]);
        gen_p.back().edge=samples[i].edge;
    }
    del_super_tets(4);
    return {gen_p,gen_tet};
}

vector<point> gen_ball_points(double radius,int target_count){
    vector<point> res;
    int guard=0;
    while(res.size()<target_count&&guard<target_count*1000){
        guard++;
        double x=(rand()%20000/10000.0-1.0)*radius;
        double y=(rand()%20000/10000.0-1.0)*radius;
        double z=(rand()%20000/10000.0-1.0)*radius;
        point p={x,y,z};
        if(norm(p)>radius)continue;
        int ok=1;
        for(auto& q:res){
            if(distance(p,q)<r_avg){
                ok=0;
                break;
            }
        }
        if(ok)res.push_back(p);
    }
    return res;
}

int main(int argc,char** argv){
    double radius=1.0;
    int target_points=1700;
    unsigned seed=1;
    std::string output="ball_mesh.txt";

    for(int i=1;i<argc;i++){
        std::string arg=argv[i];
        if(arg=="--radius"&&i+1<argc){
            radius=std::atof(argv[++i]);
        }
        else if((arg=="--points"||arg=="--samples")&&i+1<argc){
            target_points=std::atoi(argv[++i]);
        }
        else if(arg=="--min-distance"&&i+1<argc){
            r_avg=std::atof(argv[++i]);
        }
        else if(arg=="--seed"&&i+1<argc){
            seed=(unsigned)std::atoi(argv[++i]);
        }
        else if((arg=="--output"||arg=="-o")&&i+1<argc){
            output=argv[++i];
        }
        else if(arg=="--help"||arg=="-h"){
            std::cout<<"Usage: gen.exe [--radius R] [--points N] [--min-distance D] [--seed S] [--output FILE]\n";
            return 0;
        }
        else{
            std::cerr<<"Unknown argument: "<<arg<<"\n";
            return 1;
        }
    }

    std::srand(seed);
    vector<point> samples=gen_ball_points(radius,target_points);
    RES3D mesh=tetrahedron_gen(samples);

    std::ofstream fout(output);
    if(!fout){
        std::cerr<<"Cannot open output file: "<<output<<"\n";
        return 1;
    }
    fout<<mesh.pr.size()<<" "<<mesh.tetr.size()<<"\n";
    for(const auto& p:mesh.pr){
        fout<<p.x<<" "<<p.y<<" "<<p.z<<"\n";
    }
    for(const auto& t:mesh.tetr){
        fout<<t.ia<<" "<<t.ib<<" "<<t.ic<<" "<<t.id<<"\n";
    }

    std::cout<<"Generated "<<mesh.pr.size()<<" vertices and "<<mesh.tetr.size()<<" tetrahedra to "<<output<<"\n";
    return 0;
}
