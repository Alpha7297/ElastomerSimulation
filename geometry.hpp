#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP
#define EPS 1e-10
#define PI 3.141592653589
#include<cmath>
#include<vector>

class point{
public:
    double x,y,z;
    int edge;
    std::vector<int> tri_idx;
    point(double xx=0,double yy=0,double zz=0):x(xx),y(yy),z(zz),edge(0){}
    double dist(point& other)const{
        return sqrt(pow(x-other.x,2)+pow(y-other.y,2)+pow(z-other.z,2));
    }
    point operator+(const point& other)const{
        return {x+other.x,y+other.y,z+other.z};
    }
    point operator-(const point& other)const{
        return {x-other.x,y-other.y,z-other.z};
    }
    point operator/(const double& other)const{
        return {x/other,y/other,z/other};
    }
    point operator*(const double& other)const{
        return {x*other,y*other,z*other};
    }
};

inline double dot(point a,point b){
    return a.x*b.x+a.y*b.y+a.z*b.z;
}

inline point cross(point a,point b){
    return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};
}

inline double norm(point a){
    return sqrt(dot(a,a));
}

inline double det3(point a,point b,point c){
    return dot(a,cross(b,c));
}

inline double distance(point& a,point& b){
    return sqrt(pow(a.x-b.x,2)+pow(a.y-b.y,2)+pow(a.z-b.z,2));
}

class polygon{
public:
    std::vector<point> edges;
    polygon(std::vector<point> edgess){
        edges=edgess;
    }
    int in(const point& other){
        int count=0;
        for(int i=0;i<edges.size();i++){
            point a=edges[i];
            point b=edges[(i+1)%edges.size()];
            if(fabs((b.y-a.y)*(other.x-a.x)-(b.x-a.x)*(other.y-a.y))<EPS){
                if(other.x>=fmin(a.x,b.x)-EPS&&other.x<=fmax(a.x,b.x)+EPS){
                    if(other.y>=fmin(a.y,b.y)-EPS&&other.y<=fmax(a.y,b.y)+EPS){
                        return 1;
                    }
                }
            }
            if(((a.y>other.y)!=(b.y>other.y))){
                double x_inter=a.x+(other.y-a.y)*(b.x-a.x)/(b.y-a.y);
                if(x_inter>other.x+EPS)count++;
            }
        }
        return count%2;
    }
};

class triangle{
public:
    point a,b,c;
    int ia,ib,ic;
    double x1,x2,x3,y1,y2,y3;
    double p1[3],p2[3];
    double p3[3];
    double d;
    point o;
    double r;
    triangle(point aa=point(),point bb=point(),point cc=point(),int iaa=0,int ibb=0,int icc=0):a(aa),b(bb),c(cc),ia(iaa),ib(ibb),ic(icc){
        x1=a.x,x2=b.x,x3=c.x,y1=a.y,y2=b.y,y3=c.y;
        d=x1*(y2-y3)+y1*(x3-x2)+x2*y3-x3*y2;
        if(fabs(d)<EPS){
            o={0,0,0};
            r=0;
            return;
        }
        p1[0]=(y2-y3)/d,p1[1]=(y3-y1)/d,p1[2]=(y1-y2)/d;
        p2[0]=(x3-x2)/d,p2[1]=(x1-x3)/d,p2[2]=(x2-x1)/d;
        p3[0]=(x2*y3-x3*y2)/d,p3[1]=(x3*y1-y3*x1)/d,p3[2]=(x1*y2-x2*y1)/d;
        double temp=(a.x*(b.y-c.y)+b.x*(c.y-a.y)+c.x*(a.y-b.y))*2;
        if(fabs(temp)<EPS){
            o={0,0,0};
            r=0;
            return;
        }
        double a1=a.x*a.x+a.y*a.y;
        double a2=b.x*b.x+b.y*b.y;
        double a3=c.x*c.x+c.y*c.y;
        o.x=(a1*(b.y-c.y)+a2*(c.y-a.y)+a3*(a.y-b.y))/temp;
        o.y=(a1*(c.x-b.x)+a2*(a.x-c.x)+a3*(b.x-a.x))/temp;
        o.z=0;
        r=sqrt((o.x-a.x)*(o.x-a.x)+(o.y-a.y)*(o.y-a.y));
    }
};

class tetrahedron{
public:
    point a,b,c,d;
    int ia,ib,ic,id;
    point o;
    double r;
    double volume;
    tetrahedron(point aa=point(),point bb=point(),point cc=point(),point dd=point(),int iaa=0,int ibb=0,int icc=0,int idd=0):a(aa),b(bb),c(cc),d(dd),ia(iaa),ib(ibb),ic(icc),id(idd){
        point ba=b-a;
        point ca=c-a;
        point da=d-a;
        double orient=det3(ba,ca,da);
        volume=fabs(orient)/6.0;
        if(fabs(orient)<EPS){
            o={0,0,0};
            r=0;
            return;
        }
        point row1=ba*2.0;
        point row2=ca*2.0;
        point row3=da*2.0;
        point rhs={dot(b,b)-dot(a,a),dot(c,c)-dot(a,a),dot(d,d)-dot(a,a)};
        double mdet=det3(row1,row2,row3);
        if(fabs(mdet)<EPS){
            o={0,0,0};
            r=0;
            return;
        }
        double dx=det3({rhs.x,row1.y,row1.z},{rhs.y,row2.y,row2.z},{rhs.z,row3.y,row3.z});
        double dy=det3({row1.x,rhs.x,row1.z},{row2.x,rhs.y,row2.z},{row3.x,rhs.z,row3.z});
        double dz=det3({row1.x,row1.y,rhs.x},{row2.x,row2.y,rhs.y},{row3.x,row3.y,rhs.z});
        o={dx/mdet,dy/mdet,dz/mdet};
        r=norm(o-a);
    }
};

inline double atan4(point a){
    if(fabs(a.x)<EPS){
        if(a.y>0){
            return PI/2.0;
        }
        else{
            return 1.5*PI;
        }
    }
    if(fabs(a.y)<EPS){
        if(a.x>0){
            return 0;
        }
        else{
            return PI;
        }
    }
    if(a.x>0){
        if(a.y>0){
            return atan(a.y/a.x);
        }
        else{
            return 2*PI-atan(-a.y/a.x);
        }
    }
    else{
        return PI+atan(a.y/a.x);
    }
}
#endif
