#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <bits/stdc++.h>
#include <windows.h>
#include <glut.h>
#include "bitmap_image.hpp"

using namespace std;

#define pi (2*acos(0.0))
#define epsilon (1.0e-6)





int flag_space;
//Texture image
bitmap_image Given_Photo,Image2;

//Intersecting t, object, index
double IntersectingT;
int  IntersectingObj, IntersectingIndex;
//IntersectingObj == 1: CHecker, 2: Sphere, 3: Triangle
double nearDis, farDias, fovY_axis, aspect_ro;
int level_rec, num_pxl;
double ckr_width;
double ckr_am, ckr_diff, ckr_ref;
int num_of_objs;

//Sphere er jnno
double cenX,cenY,cenZ;
double radiuss;

//pyramid er jnno
double low_posX,low_posY,low_posZ;
double hgt,wdt;

//sokol object er jnno
double colrX, colrY, colrZ;
double am, diff, spec, refl;
int shini;

//light src
int num_of_light;
double light_srcX, light_srcY, light_srcZ, light_fall_off;

//spot light src
int num_of_sptLight;
double spt_srcX, spt_srcY, spt_srcZ, spt_fall_off;
double lk_X,lk_Y,lk_Z;
double cutt_off;


class Color {
public:
    double r, g, b;
    Color(double r, double g, double b) {
        this->r = r;
        this->g = g;
        this->b = b;
    }
    Color() {
    }
};
Color checker_color;



struct point
{
	double x,y,z;
};



struct point eye_pos, l, r, u;
void drawSphere(double radius,int slices,int stacks);

struct point pnt_addition(struct point X, struct point Y)
{
    struct point result;
    result.x=X.x+Y.x;
    result.y=X.y+Y.y;
    result.z=X.z+Y.z;

    return result;
};

struct point pnt_subtract(struct point X, struct point Y)
{
    struct point result;
    result.x=X.x-Y.x;
    result.y=X.y-Y.y;
    result.z=X.z-Y.z;

    return result;
};

struct point pnt_multi(struct point pnt, double num)
{
    struct point result;
    result.x=num*pnt.x;
    result.y=num*pnt.y;
    result.z=num*pnt.z;

    return result;
};

struct point pnt_div(struct point pnt, double num)
{
    struct point result;
    result.x=pnt.x/num;
    result.y=pnt.y/num;
    result.z=pnt.z/num;

    return result;
};

struct point pnt_normalize(struct point X)
{
    struct point result;

    double lobdhi;
    lobdhi=sqrt(X.x*X.x+X.y*X.y+X.z*X.z);

    result.x=X.x/lobdhi;
    result.y=X.y/lobdhi;
    result.z=X.z/lobdhi;

    return result;
};

double pnt_Dot(struct point X, struct point Y)
{
    return X.x*Y.x+X.y*Y.y+X.z*Y.z;
}

struct point pnt_Cross(struct point X, struct point Y)
{
    struct point result;

    result.x=X.y*Y.z-X.z*Y.y;
    result.y=-(X.x*Y.z-X.z*Y.x);
    result.z=X.x*Y.y-X.y*Y.x;

    return result;
};

struct point pnt_rotationn(struct point refer, struct point vec, double angl)
{
    struct point perpendi=pnt_Cross(refer,vec);
    struct point result;

    result.x=vec.x*cos(angl)+perpendi.x*sin(angl);
    result.y=vec.y*cos(angl)+perpendi.y*sin(angl);
    result.z=vec.z*cos(angl)+perpendi.z*sin(angl);

    return result;
};

class Ray_line
{
public:
   struct point Ro,Rd;
   Ray_line(struct point strt, struct point direction)
   {
       Ro.x=strt.x;
       Ro.y=strt.y;
       Ro.z=strt.z;

       Rd.x=direction.x;
       Rd.y=direction.y;
       Rd.z=direction.z;

       Rd=pnt_normalize(Rd);
   }
};
void Nearest_T_Obj_Index(Ray_line RRLnn);


class Checker_Brd
{
public:
    double sizz;
    double check_am,check_diff,check_reflct;
    Color colr;
    struct point mainPoint;

    Checker_Brd()
    {

    }

    Checker_Brd(double szC, double amC, double diffC, double reflctC)
    {
        sizz=szC;
        check_am=amC;
        check_diff=diffC;
        check_reflct=reflctC;
        mainPoint={-10000.0, -10000.0, 0.0};
    }

    void Checker_Draw()
    {
        glPushMatrix();
        {
            glBegin(GL_QUADS);
            {
                for (int i = 0; i < 20000/int(sizz); i++)
                    for (int j = 0; j < 20000/int(sizz); j++)
                    {
                        if((i+j)%2==0)
                        {
                            glColor3f(1.0,1.0,1.0);
                        }
                        else{
                            glColor3f(0.0,0.0,0.0);
                        }

                        glVertex3f(mainPoint.x+(i*sizz), mainPoint.y+(j*sizz), 0.0);
                        glVertex3f(mainPoint.x+(i*sizz), mainPoint.y+((j+1)*sizz), 0.0);
                        glVertex3f(mainPoint.x+((i+1)*sizz), mainPoint.y+((j+1)*sizz), 0.0);
                        glVertex3f(mainPoint.x+((i+1)*sizz), mainPoint.y+(j*sizz), 0.0);
                    }
            }
            glEnd();
        }
        glPopMatrix();
    }
};

Checker_Brd ckrr1;

class TriAngle
{
public:
    struct point T1,T2,T3;
    Color TriColor;
    double Tri_am,Tri_diff,Tri_spec,Tri_reflec;
    int Tri_Shini;

    TriAngle(struct point p1, struct point p2, struct point p3, Color clr2, double am2, double diff2, double spec2, double reflec2, int shini2)
    {
        T1=p1;
        T2=p2;
        T3=p3;
        TriColor=clr2;
        Tri_am=am2;
        Tri_diff=diff2;
        Tri_spec=spec2;
        Tri_reflec=reflec2;
        Tri_Shini=shini2;
    }
};
vector<TriAngle> vector_of_Triangle;


class Pyramid
{
public:
    struct point low_pnt;
    double py_wid, py_height;
    Color py_colr;
    double py_am, py_diff,py_spec, py_reflec;
    int py_shini;

    Pyramid()
    {

    }
    Pyramid(double x, double y, double z, double wid,double hgt, double colX, double colY, double colZ, double am, double diff, double spec, double reflec, int shini)
    {
        low_pnt={x,y,z};
        py_wid=wid;
        py_height=hgt;
        Color ddd(colX, colY, colZ);
        py_colr=ddd;
        py_am=am;
        py_diff=diff;
        py_spec=spec;
        py_reflec=reflec;
        py_shini=shini;
    }
    void py_Draw()
    {
        struct point pnt2={low_pnt.x+1*py_wid, low_pnt.y, low_pnt.z};
        struct point pnt3={low_pnt.x+1*py_wid, low_pnt.y+1*py_wid, low_pnt.z};
        struct point pnt4={low_pnt.x, low_pnt.y+1*py_wid, low_pnt.z};

        struct point HEIGHT={low_pnt.x+0.5*py_wid, low_pnt.y+0.5*py_wid, low_pnt.z+1*py_height};

        glPushMatrix();
        {
            glColor3f(py_colr.r, py_colr.g, py_colr.b);
            glBegin(GL_TRIANGLES);
            {
                glVertex3f(HEIGHT.x,HEIGHT.y,HEIGHT.z);
                glVertex3f(low_pnt.x,low_pnt.y,low_pnt.z);
                glVertex3f(pnt2.x,pnt2.y,pnt2.z);

                glVertex3f(HEIGHT.x,HEIGHT.y,HEIGHT.z);
                glVertex3f(pnt2.x,pnt2.y,pnt2.z);
                glVertex3f(pnt3.x,pnt3.y,pnt3.z);

                glVertex3f(HEIGHT.x,HEIGHT.y,HEIGHT.z);
                glVertex3f(pnt3.x,pnt3.y,pnt3.z);
                glVertex3f(pnt4.x,pnt4.y,pnt4.z);

                glVertex3f(HEIGHT.x,HEIGHT.y,HEIGHT.z);
                glVertex3f(pnt4.x,pnt4.y,pnt4.z);
                glVertex3f(low_pnt.x,low_pnt.y,low_pnt.z);
            }
            glEnd();
        }
        glPopMatrix();


    }
};

Pyramid pyra1;

class SpHere
{
public:
    struct point centr;
    double radi;
    Color Sp_color;
    double Sp_am,Sp_diff,Sp_spec,Sp_reflec;
    int Sp_shini;
    SpHere(double x, double y, double z, double radus, double r, double g, double b, double am, double diff, double spec, double reflec, double shini)
    {
        centr={x,y,z};
        radi=radus;
        Color temp(r,g,b);
        Sp_color=temp;
        Sp_am=am;
        Sp_diff=diff;
        Sp_spec=spec;
        Sp_reflec=reflec;
        Sp_shini=shini;
    }
    void Sp_Draw()
    {
        glPushMatrix();
        {
            glTranslatef(centr.x,centr.y,centr.z);
            glColor3f(Sp_color.r,Sp_color.g,Sp_color.b);
            drawSphere(radi,50,50);
        }
        glPopMatrix();
    }
};

class light_src
{
public:
    struct point lt_positn;
    double lt_fall_off;
    light_src(double x, double y, double z, double fll_off)
    {
        lt_positn={x,y,z};
        lt_fall_off=fll_off;
    }
    void lt_draw()
    {
        glPushMatrix();
        {
            glTranslatef(lt_positn.x,lt_positn.y,lt_positn.z);
            glColor3f(1.0,1.0,1.0);
            drawSphere(5,50,50);
        }
        glPopMatrix();

    }
};

class spot_light_src
{
public:
    struct point spt_positn;
    double spt_falll;
    struct point look;
    double spt_cutt_off;
    spot_light_src(double x, double y, double z, double fll_off, double lkX, double lkY, double lkZ, double ct_off)
    {
        spt_positn={x,y,z};
        spt_falll=fll_off;
        look={lkX,lkY,lkZ};
        spt_cutt_off=ct_off;
    }
    void spt_draw()
    {
        glPushMatrix();
        {
            glTranslatef(spt_positn.x,spt_positn.y,spt_positn.z);
            glColor3f(1.0,1.0,1.0);
            drawSphere(5,50,50);
        }
        glPopMatrix();

    }

};

vector<Checker_Brd> vector_of_Chkr;
vector<Pyramid> vector_of_pyra;
vector<SpHere> vector_of_Sphere;
vector<light_src> vector_of_light;
vector<spot_light_src> vector_of_spot_light;

struct point Reflection(struct point N, Ray_line L)
{
    N=pnt_normalize(N);
    struct point R=pnt_normalize(pnt_subtract(L.Rd, pnt_multi(N,2.0*pnt_Dot(L.Rd,N))));
    return R;

};

//xy palne intersection, z_val=0
struct point Line_Plane_Intersect(Ray_line R, double z_val)
{
    double t=(z_val-R.Ro.z)/R.Rd.z;
    struct point result;

    result.x=R.Ro.x+t*R.Rd.x;
    result.y=R.Ro.y+t*R.Rd.y;
    result.z=R.Ro.z+t*R.Rd.z;

    return result;
}

double Line_Plane_Intersect_TVAl(Ray_line R, double z_val)
{
    if(R.Rd.z!=0){
        double t=(z_val-R.Ro.z)/R.Rd.z;
        return t-epsilon;
    }
    return -9999999;
}

struct point Plane_Normal()
{
    struct point result={0.0,0.0,0.1};
    return result;
};


//Sphere Intersection
double Ray_Sphere_Intersect_TVAL(Ray_line R, SpHere spr)
{
    struct point CenToRo=pnt_subtract(R.Ro,spr.centr);

    double b=2*pnt_Dot(R.Rd,CenToRo);
    double c=pnt_Dot(CenToRo,CenToRo)-spr.radi*spr.radi;
    double deter=b*b-4*c;
    //if(deter>=0) cout<<"Determinamt"<<deter<<endl;

    if(deter==0)
    {
        cout<<"T "<<-(b/2)-epsilon<<endl;
        return -(b/2)-epsilon;
    }
    else if(deter>0)
    {
        double tt=(-b+sqrt(deter))/2.0;
        double tt1=(-b-sqrt(deter))/2.0;

         return min(tt,tt1);
    }
    return -9999999;
};

struct point Ray_Sphere_Intersect(Ray_line R, SpHere spr)
{
    double t=Ray_Sphere_Intersect_TVAL(R,spr);
    struct point result;

    result.x=R.Ro.x+t*R.Rd.x;
    result.y=R.Ro.y+t*R.Rd.y;
    result.z=R.Ro.z+t*R.Rd.z;

    return result;
};

struct point Sphere_Normal(SpHere spr, struct point intrsct_pnt)
{
    struct point result=pnt_subtract(intrsct_pnt,spr.centr);
    return pnt_normalize(result);
};

//Triangle intersect
double Ray_Triangle_Intersect_TVAL(Ray_line R, TriAngle TA)
{
    double axbx=TA.T1.x-TA.T2.x;
    double ayby=TA.T1.y-TA.T2.y;
    double azbz=TA.T1.z-TA.T2.z;

    double axcx=TA.T1.x-TA.T3.x;
    double aycy=TA.T1.y-TA.T3.y;
    double azcz=TA.T1.z-TA.T3.z;

    double axRox=TA.T1.x-R.Ro.x;
    double ayRoy=TA.T1.y-R.Ro.y;
    double azRoz=TA.T1.z-R.Ro.z;

    double A=(axbx)*(aycy*R.Rd.z-azcz*R.Rd.y)-(axcx)*(ayby*R.Rd.z-azbz*R.Rd.y)+R.Rd.x*(ayby*azcz-aycy*azbz);

    double beta=((axRox)*(aycy*R.Rd.z-azcz*R.Rd.y)-(axcx)*(ayRoy*R.Rd.z-azRoz*R.Rd.y)+(R.Rd.x)*(ayRoy*azcz-aycy*azRoz))/A;
    double gama=((axbx)*(ayRoy*R.Rd.z-azRoz*R.Rd.y)-(axRox)*(ayby*R.Rd.z-azbz*R.Rd.y)+(R.Rd.x)*(ayby*azRoz-azbz*ayRoy))/A;

    if(beta+gama<1 && beta>0 && gama>0)
    {
        double t=((axbx)*(aycy*azRoz-azcz*ayRoy)-(axcx)*(ayby*azRoz-azbz*ayRoy)+(axRox)*(ayby*azcz-azbz*aycy))/A;
        return t;

    }
    return -9999999;
}

struct point Ray_Triangle_Intersect(Ray_line R, TriAngle TA)
{
    double t=Ray_Triangle_Intersect_TVAL(R,TA);

    struct point result={R.Ro.x+t*R.Rd.x, R.Ro.y+t*R.Rd.y, R.Ro.z+t*R.Rd.z};
    return result;
};

struct point Triangle_Normal(TriAngle TA)
{
    struct point vec1=pnt_subtract(TA.T1,TA.T2);
    struct point vec2=pnt_subtract(TA.T3,TA.T2);
    struct point kross=pnt_Cross(vec2,vec1);
    kross=pnt_normalize(kross);
    return kross;
};

bool is_lt_Illuminte(light_src Src, struct point intersecc)
{
    struct point S=Src.lt_positn;

    struct point PS=pnt_subtract(S,intersecc);
    PS=pnt_normalize(PS);//Dirction
    struct point strt=pnt_addition(intersecc,pnt_multi(PS,epsilon));

    Ray_line Rl(strt,PS);

    double kk;
    for(int i=0;i<vector_of_Triangle.size();i++)
    {
        kk=Ray_Triangle_Intersect_TVAL(Rl,vector_of_Triangle[i]);
        if(kk>0) return false;
    }
    for(int i=0;i<vector_of_Sphere.size();i++)
    {
        kk=Ray_Sphere_Intersect_TVAL(Rl,vector_of_Sphere[i]);
        if(kk>0) return false;
    }

    return true;
}

bool is_Spt_Illuminate(spot_light_src Spot_srcc, struct point intersecc)
{
    struct point S=Spot_srcc.spt_positn;
    struct point SP=pnt_subtract(intersecc,S);
    SP=pnt_normalize(SP);

    struct point PS=pnt_subtract(S, intersecc);
    PS=pnt_normalize(PS);

    struct point look_dir=pnt_subtract(Spot_srcc.look,S);
    look_dir=pnt_normalize(look_dir);

    double anglll=acos(pnt_Dot(SP,look_dir));
    anglll=anglll*360.0/(2*pi);
    if(anglll>Spot_srcc.spt_cutt_off) return false;

    struct point strt=pnt_addition(intersecc,pnt_multi(PS,epsilon));

    Ray_line Rl(strt,PS);

    double kk;
    for(int i=0;i<vector_of_Triangle.size();i++)
    {
        kk=Ray_Triangle_Intersect_TVAL(Rl,vector_of_Triangle[i]);
        if(kk>0) return false;
    }
    for(int i=0;i<vector_of_Sphere.size();i++)
    {
        kk=Ray_Sphere_Intersect_TVAL(Rl,vector_of_Sphere[i]);
        if(kk>0) return false;
    }

    return true;


}

Color Am_Diff_Spec_Reflec(Ray_line RL, struct point intersecc, int level_num, int objType, int index)
{
    double lambrt=0, phng=0;
    Color result(0.0,0.0,0.0);

    if(objType==1)
    {
        Checker_Brd ckkr=vector_of_Chkr[index];
        for(int i=0;i<vector_of_light.size();i++)
        {
            //cout<<"Ahee1"<<endl;
            int flagg=0;
            if(is_lt_Illuminte(vector_of_light[i],intersecc)==false) flagg=1;
            //cout<<is_lt_Illuminte(vector_of_light[i],intersecc)<<endl;
            struct point PS=pnt_subtract(vector_of_light[i].lt_positn,intersecc);
            double dstnc=sqrt(pnt_Dot(PS,PS));
            PS=pnt_normalize(PS);

            struct point Normall=Plane_Normal();
            Normall=pnt_normalize(Normall);

            double scl_factor=exp(-dstnc*dstnc*vector_of_light[i].lt_fall_off);
            lambrt=pnt_Dot(PS,Normall)*scl_factor;

            struct point R_bar=Reflection(Normall,RL);
            R_bar=pnt_normalize(R_bar);

            //phng+=pow(pnt_Dot(R_bar,PS),ckkr.)*scl_factor;

//            if(lambrt>1) lambrt=1;
//            if(lambrt<0) lambrt=0;
//
//            if(phng>1) phng=1;
//            if(phng<0) phng=0;

            if(flagg==0){
                result.r+=checker_color.r*(ckkr.check_am+ckkr.check_diff*lambrt);
                result.g+=checker_color.g*(ckkr.check_am+ckkr.check_diff*lambrt);
                result.b+=checker_color.b*(ckkr.check_am+ckkr.check_diff*lambrt);
            }
            else if(flagg==1)
            {
                result.r+=checker_color.r*(ckkr.check_am);
                result.g+=checker_color.g*(ckkr.check_am);
                result.b+=checker_color.b*(ckkr.check_am);
            }

            if(level_num<level_rec)
            {
                struct point sttrt=pnt_addition(intersecc,pnt_multi(R_bar,epsilon));
                Ray_line RRLnn(sttrt,R_bar);

                Nearest_T_Obj_Index(RRLnn);
                if(IntersectingT!=-1)
                {

                    struct point inns=pnt_addition(RRLnn.Ro,pnt_multi(RRLnn.Rd,IntersectingT));
                    struct point intrr;
                    if(IntersectingObj==1)
                    {
                        intrr=Line_Plane_Intersect(RRLnn,0);
                    }
                    else if(IntersectingObj==2)
                    {
                        intrr=Ray_Sphere_Intersect(RRLnn,vector_of_Sphere[IntersectingIndex]);
                    }
                    else
                    {
                        intrr=Ray_Triangle_Intersect(RRLnn,vector_of_Triangle[IntersectingIndex]);
                    }
                    Color tmmpp=Am_Diff_Spec_Reflec(RRLnn,intrr,level_num+1,IntersectingObj,IntersectingIndex);
                    result.r+=tmmpp.r*ckkr.check_reflct;
                    result.g+=tmmpp.g*ckkr.check_reflct;
                    result.b+=tmmpp.b*ckkr.check_reflct;
                }
            }

        }

        for(int i=0;i<vector_of_spot_light.size();i++)
        {
            //cout<<"Ahee2"<<endl;
            int flagg=0;
            if(is_Spt_Illuminate(vector_of_spot_light[i],intersecc)==false) flagg=1;
            //cout<<is_Spt_Illuminate(vector_of_spot_light[i],intersecc)<<endl;
            struct point PS=pnt_subtract(vector_of_spot_light[i].spt_positn,intersecc);
            double dstnc=sqrt(pnt_Dot(PS,PS));
            PS=pnt_normalize(PS);

            struct point Normall=Plane_Normal();
            Normall=pnt_normalize(Normall);

            double scl_factor=exp(-dstnc*dstnc*vector_of_spot_light[i].spt_falll);
            lambrt=pnt_Dot(PS,Normall)*scl_factor;

            struct point R_bar=Reflection(Normall,RL);
            R_bar=pnt_normalize(R_bar);

            //phng+=pow(pnt_Dot(R_bar,PS),spr.Sp_shini)*scl_factor;

            if(lambrt>1) lambrt=1;
            if(lambrt<0) lambrt=0;


            if(flagg==0){
                result.r+=checker_color.r*(ckkr.check_am+ckkr.check_diff*lambrt);
                result.g+=checker_color.g*(ckkr.check_am+ckkr.check_diff*lambrt);
                result.b+=checker_color.b*(ckkr.check_am+ckkr.check_diff*lambrt);
            }
            else if(flagg==1)
            {
                result.r+=checker_color.r*(ckkr.check_am);
                result.g+=checker_color.g*(ckkr.check_am);
                result.b+=checker_color.b*(ckkr.check_am);
            }

            if(level_num<level_rec)
            {
                struct point sttrt=pnt_addition(intersecc,pnt_multi(R_bar,epsilon));
                Ray_line RRLnn(sttrt,R_bar);

                Nearest_T_Obj_Index(RRLnn);
                if(IntersectingT!=-1)
                {
                    struct point inns=pnt_addition(RRLnn.Ro,pnt_multi(RRLnn.Rd,IntersectingT));
                    struct point intrr;
                    if(IntersectingObj==1)
                    {
                        intrr=Line_Plane_Intersect(RRLnn,0);
                    }
                    else if(IntersectingObj==2)
                    {
                        intrr=Ray_Sphere_Intersect(RRLnn,vector_of_Sphere[IntersectingIndex]);
                    }
                    else
                    {
                        intrr=Ray_Triangle_Intersect(RRLnn,vector_of_Triangle[IntersectingIndex]);
                    }
                    Color tmmpp=Am_Diff_Spec_Reflec(RRLnn,intrr,level_num+1,IntersectingObj,IntersectingIndex);
                    //Color tmmpp=Am_Diff_Spec_Reflec(RRLnn,inns,level_num+1,IntersectingObj,IntersectingIndex);
                    result.r+=tmmpp.r*ckkr.check_reflct;
                    result.g+=tmmpp.g*ckkr.check_reflct;
                    result.b+=tmmpp.b*ckkr.check_reflct;
                }
            }

        }

    }
    else if(objType==2)
    {
        SpHere spr=vector_of_Sphere[index];
        for(int i=0;i<vector_of_light.size();i++)
        {
            int flagg=0;
            //cout<<"Ahee1"<<endl;
            if(is_lt_Illuminte(vector_of_light[i],intersecc)==false) flagg=1;
            //cout<<is_lt_Illuminte(vector_of_light[i],intersecc)<<endl;
            //struct point Sourcc=vector_of_light[i].lt_positn;
            struct point PS=pnt_subtract(vector_of_light[i].lt_positn,intersecc);
            double dstnc=sqrt(pnt_Dot(PS,PS));
            PS=pnt_normalize(PS);

            struct point Normall=Sphere_Normal(spr,intersecc);
            Normall=pnt_normalize(Normall);

            double scl_factor=exp(-dstnc*dstnc*vector_of_light[i].lt_fall_off);
            lambrt+=pnt_Dot(PS,Normall)*scl_factor;

            struct point R_bar=Reflection(Normall,RL);
            R_bar=pnt_normalize(R_bar);

            phng+=pow(pnt_Dot(R_bar,PS),spr.Sp_shini)*scl_factor;

//            if(lambrt>1) lambrt=1;
//            if(lambrt<0) lambrt=0;
//
//            if(phng>1) phng=1;
//            if(phng<0) phng=0;

            if(flagg==0){
                result.r+=spr.Sp_color.r*(spr.Sp_am+spr.Sp_diff*lambrt+spr.Sp_spec*phng);
                result.g+=spr.Sp_color.g*(spr.Sp_am+spr.Sp_diff*lambrt+spr.Sp_spec*phng);
                result.b+=spr.Sp_color.b*(spr.Sp_am+spr.Sp_diff*lambrt+spr.Sp_spec*phng);
            }
            else if(flagg==1)
            {
                result.r+=spr.Sp_color.r*(spr.Sp_am);
                result.g+=spr.Sp_color.g*(spr.Sp_am);
                result.b+=spr.Sp_color.b*(spr.Sp_am);
            }

            if(level_num<level_rec)
            {
                struct point sttrt=pnt_addition(intersecc,pnt_multi(R_bar,epsilon));
                Ray_line RRLnn(sttrt,R_bar);

                Nearest_T_Obj_Index(RRLnn);
                if(IntersectingT!=-1)
                {
                    struct point inns=pnt_addition(RRLnn.Ro,pnt_multi(RRLnn.Rd,IntersectingT));
                    struct point intrr;
                    if(IntersectingObj==1)
                    {
                        intrr=Line_Plane_Intersect(RRLnn,0);
                    }
                    else if(IntersectingObj==2)
                    {
                        intrr=Ray_Sphere_Intersect(RRLnn,vector_of_Sphere[IntersectingIndex]);
                    }
                    else
                    {
                        intrr=Ray_Triangle_Intersect(RRLnn,vector_of_Triangle[IntersectingIndex]);
                    }
                    Color tmmpp=Am_Diff_Spec_Reflec(RRLnn,intrr,level_num+1,IntersectingObj,IntersectingIndex);
                    //Color tmmpp=Am_Diff_Spec_Reflec(RRLnn,inns,level_num+1,IntersectingObj,IntersectingIndex);
                    result.r+=tmmpp.r*spr.Sp_reflec;
                    result.g+=tmmpp.g*spr.Sp_reflec;
                    result.b+=tmmpp.b*spr.Sp_reflec;
                }
            }

        }

        for(int i=0;i<vector_of_spot_light.size();i++)
        {
            //cout<<"Ahee2"<<endl;
            int flagg=0;
            if(is_Spt_Illuminate(vector_of_spot_light[i],intersecc)==false) flagg=1;
            //cout<<is_Spt_Illuminate(vector_of_spot_light[i],intersecc)<<endl;
            struct point PS=pnt_subtract(vector_of_spot_light[i].spt_positn,intersecc);
            double dstnc=sqrt(pnt_Dot(PS,PS));
            PS=pnt_normalize(PS);

            struct point Normall=Sphere_Normal(spr,intersecc);
            Normall=pnt_normalize(Normall);

            double scl_factor=exp(-dstnc*dstnc*vector_of_spot_light[i].spt_falll);
            lambrt+=pnt_Dot(PS,Normall)*scl_factor;

            struct point R_bar=Reflection(Normall,RL);
            R_bar=pnt_normalize(R_bar);

            phng+=pow(pnt_Dot(R_bar,PS),spr.Sp_shini)*scl_factor;

//            if(lambrt>1) lambrt=1;
//            if(lambrt<0) lambrt=0;
//
//            if(phng>1) phng=1;
//            if(phng<0) phng=0;

            if(flagg==0){
                result.r+=spr.Sp_color.r*(spr.Sp_am+spr.Sp_diff*lambrt+spr.Sp_spec*phng);
                result.g+=spr.Sp_color.g*(spr.Sp_am+spr.Sp_diff*lambrt+spr.Sp_spec*phng);
                result.b+=spr.Sp_color.b*(spr.Sp_am+spr.Sp_diff*lambrt+spr.Sp_spec*phng);
            }
            else if(flagg==1)
            {
                result.r+=spr.Sp_color.r*(spr.Sp_am);
                result.g+=spr.Sp_color.g*(spr.Sp_am);
                result.b+=spr.Sp_color.b*(spr.Sp_am);
            }

            if(level_num<level_rec)
            {
                struct point sttrt=pnt_addition(intersecc,pnt_multi(R_bar,epsilon));
                Ray_line RRLnn(sttrt,R_bar);

                Nearest_T_Obj_Index(RRLnn);
                if(IntersectingT!=-1)
                {
                    struct point inns=pnt_addition(RRLnn.Ro,pnt_multi(RRLnn.Rd,IntersectingT));
                    struct point intrr;
                    if(IntersectingObj==1)
                    {
                        intrr=Line_Plane_Intersect(RRLnn,0);
                    }
                    else if(IntersectingObj==2)
                    {
                        intrr=Ray_Sphere_Intersect(RRLnn,vector_of_Sphere[IntersectingIndex]);
                    }
                    else
                    {
                        intrr=Ray_Triangle_Intersect(RRLnn,vector_of_Triangle[IntersectingIndex]);
                    }
                    Color tmmpp=Am_Diff_Spec_Reflec(RRLnn,intrr,level_num+1,IntersectingObj,IntersectingIndex);
                    //Color tmmpp=Am_Diff_Spec_Reflec(RRLnn,inns,level_num+1,IntersectingObj,IntersectingIndex);
                    result.r+=tmmpp.r*spr.Sp_reflec;
                    result.g+=tmmpp.g*spr.Sp_reflec;
                    result.b+=tmmpp.b*spr.Sp_reflec;
                }
            }

        }
    }
    else if(objType==3)
    {
        TriAngle tangl=vector_of_Triangle[index];
        for(int i=0;i<vector_of_light.size();i++)
        {
            //cout<<"Ahee1"<<endl;
            int flagg=0;
            if(is_lt_Illuminte(vector_of_light[i],intersecc)==false) flagg=1;
            //cout<<is_lt_Illuminte(vector_of_light[i],intersecc)<<endl;
            //struct point Sourcc=vector_of_light[i].lt_positn;
            struct point PS=pnt_subtract(vector_of_light[i].lt_positn,intersecc);
            double dstnc=sqrt(pnt_Dot(PS,PS));
            PS=pnt_normalize(PS);

            struct point Normall=Triangle_Normal(tangl);
            Normall=pnt_normalize(Normall);

            double scl_factor=exp(-dstnc*dstnc*vector_of_light[i].lt_fall_off);
            lambrt+=pnt_Dot(PS,Normall)*scl_factor;

            struct point R_bar=Reflection(Normall,RL);
            R_bar=pnt_normalize(R_bar);


            phng+=pow(pnt_Dot(R_bar,PS),tangl.Tri_Shini)*scl_factor;
            //phng=pow(pnt_Dot(R,V),tangl.Tri_Shini)*scl_factor;

//            if(lambrt>1) lambrt=1;
//            if(lambrt<0) lambrt=0;
//
//            if(phng>1) phng=1;
//            if(phng<0) phng=0;

            if(flagg==0){
                result.r+=tangl.TriColor.r*(tangl.Tri_am+tangl.Tri_diff*lambrt+tangl.Tri_spec*phng);
                result.g+=tangl.TriColor.g*(tangl.Tri_am+tangl.Tri_diff*lambrt+tangl.Tri_spec*phng);
                result.b+=tangl.TriColor.b*(tangl.Tri_am+tangl.Tri_diff*lambrt+tangl.Tri_spec*phng);
            }
            else if(flagg==1)
            {
                result.r+=tangl.TriColor.r*(tangl.Tri_am);
                result.g+=tangl.TriColor.g*(tangl.Tri_am);
                result.b+=tangl.TriColor.b*(tangl.Tri_am);
            }

            if(level_num<level_rec)
            {
                struct point sttrt=pnt_addition(intersecc,pnt_multi(R_bar,epsilon));
                Ray_line RRLnn(sttrt,R_bar);

                Nearest_T_Obj_Index(RRLnn);
                if(IntersectingT!=-1)
                {
                    struct point inns=pnt_addition(RRLnn.Ro,pnt_multi(RRLnn.Rd,IntersectingT));
                    struct point intrr;
                    if(IntersectingObj==1)
                    {
                        intrr=Line_Plane_Intersect(RRLnn,0);
                    }
                    else if(IntersectingObj==2)
                    {
                        intrr=Ray_Sphere_Intersect(RRLnn,vector_of_Sphere[IntersectingIndex]);
                    }
                    else
                    {
                        intrr=Ray_Triangle_Intersect(RRLnn,vector_of_Triangle[IntersectingIndex]);
                    }
                    Color tmmpp=Am_Diff_Spec_Reflec(RRLnn,intrr,level_num+1,IntersectingObj,IntersectingIndex);
                    //Color tmmpp=Am_Diff_Spec_Reflec(RRLnn,inns,level_num+1,IntersectingObj,IntersectingIndex);
                    result.r+=tmmpp.r*tangl.Tri_reflec;
                    result.g+=tmmpp.g*tangl.Tri_reflec;
                    result.b+=tmmpp.b*tangl.Tri_reflec;
                }
            }

        }

        for(int i=0;i<vector_of_spot_light.size();i++)
        {
            //cout<<"Ahee2"<<endl;
            int flagg=0;
            if(is_Spt_Illuminate(vector_of_spot_light[i],intersecc)==false) flagg=1;
            //cout<<is_Spt_Illuminate(vector_of_spot_light[i],intersecc)<<endl;
            struct point PS=pnt_subtract(vector_of_spot_light[i].spt_positn,intersecc);
            double dstnc=sqrt(pnt_Dot(PS,PS));
            PS=pnt_normalize(PS);

            struct point Normall=Triangle_Normal(tangl);
            Normall=pnt_normalize(Normall);

            double scl_factor=exp(-dstnc*dstnc*vector_of_spot_light[i].spt_falll);
            lambrt+=pnt_Dot(PS,Normall)*scl_factor;

            struct point R_bar=Reflection(Normall,RL);
            R_bar=pnt_normalize(R_bar);

            phng+=pow(pnt_Dot(R_bar,PS),tangl.Tri_Shini)*scl_factor;

//            if(lambrt>1) lambrt=1;
//            if(lambrt<0) lambrt=0;
//
//            if(phng>1) phng=1;
//            if(phng<0) phng=0;

            if(flagg==0){
                result.r+=tangl.TriColor.r*(tangl.Tri_am+tangl.Tri_diff*lambrt+tangl.Tri_spec*phng);
                result.g+=tangl.TriColor.g*(tangl.Tri_am+tangl.Tri_diff*lambrt+tangl.Tri_spec*phng);
                result.b+=tangl.TriColor.b*(tangl.Tri_am+tangl.Tri_diff*lambrt+tangl.Tri_spec*phng);
            }
            else if(flagg==1)
            {
                result.r+=tangl.TriColor.r*(tangl.Tri_am);
                result.g+=tangl.TriColor.g*(tangl.Tri_am);
                result.b+=tangl.TriColor.b*(tangl.Tri_am);
            }

            if(level_num<level_rec)
            {
                struct point sttrt=pnt_addition(intersecc,pnt_multi(R_bar,epsilon));
                Ray_line RRLnn(sttrt,R_bar);

                Nearest_T_Obj_Index(RRLnn);
                if(IntersectingT!=-1)
                {
                    struct point inns=pnt_addition(RRLnn.Ro,pnt_multi(RRLnn.Rd,IntersectingT));
                    struct point intrr;
                    if(IntersectingObj==1)
                    {
                        intrr=Line_Plane_Intersect(RRLnn,0);
                    }
                    else if(IntersectingObj==2)
                    {
                        intrr=Ray_Sphere_Intersect(RRLnn,vector_of_Sphere[IntersectingIndex]);
                    }
                    else
                    {
                        intrr=Ray_Triangle_Intersect(RRLnn,vector_of_Triangle[IntersectingIndex]);
                    }
                    Color tmmpp=Am_Diff_Spec_Reflec(RRLnn,intrr,level_num+1,IntersectingObj,IntersectingIndex);
                    //Color tmmpp=Am_Diff_Spec_Reflec(RRLnn,inns,level_num+1,IntersectingObj,IntersectingIndex);
                    result.r+=tmmpp.r*tangl.Tri_reflec;
                    result.g+=tmmpp.g*tangl.Tri_reflec;
                    result.b+=tmmpp.b*tangl.Tri_reflec;
                }
            }

        }
    }
    if(result.r>1) result.r=1;
    if(result.r<0) result.r=0;

    if(result.g>1) result.g=1;
    if(result.g<0) result.g=0;

    if(result.b>1) result.b=1;
    if(result.b<0) result.b=0;
    return result;
}


void Nearest_T_Obj_Index(Ray_line R)
{
    IntersectingT=-1;
    IntersectingIndex=-1;
    IntersectingObj=0;

    double tmin=1000000;
    int indx=-1;

    for(int i=0;i<vector_of_Chkr.size();i++)
    {
        if(R.Rd.z!=0){
            double tempT=Line_Plane_Intersect_TVAl(R,0);
            struct point temp=Line_Plane_Intersect(R,0);

            int pxl_x = floor((temp.x - (-10000.0)) / ckr_width);
            int pxl_y = floor((temp.y - (-10000.0)) / ckr_width);
            int c = (pxl_x + pxl_y+1) % 2;
            Color ttt(c,c,c);
            if(tempT>0 && tempT<tmin)
            {
                tmin=tempT;
                indx=i;

                IntersectingObj=1;
                IntersectingT=tmin;
                IntersectingIndex=indx;

                if(flag_space==1)
                {
                    unsigned char r3, g3, b3;
                    int i3, j3;
                    i3 =(int)((int)temp.x - (-10000.0))%(Image2.width()*3/4);
                    j3 = (int)((int)temp.y - (-10000.0))%(Image2.height()*3/4);
//                    i3=temp.x-(pxl_x*ckr_width);
//                    j3=temp.y-(pxl_y*ckr_width);

                    Image2.get_pixel(i3, j3, r3, g3, b3);



                    Color tdd((double(r3) / 255.0), (double(g3) / 255.0),(double(b3) / 255.0));
                    //Color t11(1.0,1.0,1.0);
                    checker_color=tdd;
                    //flag_space=0;


                }
                else{
                    checker_color=ttt;
                }
            }

        }
    }

    for(int i=0;i<vector_of_Sphere.size();i++)
    {
        double tempT=Ray_Sphere_Intersect_TVAL(R,vector_of_Sphere[i]);
        //cout<<"Sphere "<<tempT<<endl;
        if(tempT!=-9999999 && tempT<tmin)
        {
            tmin=tempT;
            indx=i;

            IntersectingObj=2;
            IntersectingT=tmin;
            IntersectingIndex=indx;
        }
    }

    for(int i=0;i<vector_of_Triangle.size();i++)
    {
        double tempT=Ray_Triangle_Intersect_TVAL(R,vector_of_Triangle[i]);
        if(tempT!=-9999999 && tempT<tmin)
        {
            tmin=tempT;
            indx=i;

            IntersectingObj=3;
            IntersectingT=tmin;
            IntersectingIndex=indx;
        }

    }
}

void Building_Image()
{
    Color** TextBuffer;
    TextBuffer=new Color*[num_pxl];


    double near_dist=((double)num_pxl)/(2*tan(fovY_axis*2*pi/(360*2)));
    struct point middle_pnt=pnt_addition(eye_pos,pnt_multi(l,near_dist));
    struct point ekdom_left_corner=pnt_addition(pnt_subtract(middle_pnt,pnt_multi(r,(double)num_pxl/2)),pnt_multi(u,(double)num_pxl/2));

    for(int i=0;i<num_pxl;i++)
    {
        TextBuffer[i]=new Color[num_pxl];
        for(int j=0;j<num_pxl;j++)
        {
            struct point pxl=pnt_subtract(pnt_addition(ekdom_left_corner,pnt_multi(r,(double)i)),pnt_multi(u,(double)j));
            Ray_line eye_to_pixel(eye_pos,pnt_subtract(pxl,eye_pos));

            Nearest_T_Obj_Index(eye_to_pixel);
            Color blk(0.0,0.0,0.0);

            //cout<<"T er val"<<IntersectingT<<endl;

            if(IntersectingT!=-1)
            {
                if(IntersectingObj==1)//CHecker board
                {
//                    TextBuffer[i][j].r=1;
//                    TextBuffer[i][j].g=1;
//                    TextBuffer[i][j].b=1;
                    struct point intrsc=Line_Plane_Intersect(eye_to_pixel,0);
                    TextBuffer[i][j]=Am_Diff_Spec_Reflec(eye_to_pixel,intrsc,1,1,IntersectingIndex);

                    if(flag_space==5)
                    {
                        unsigned char r3, g3, b3;
                        int i3, j3;

                        Image2.get_pixel(i%(int)ckr_width,j%(int)ckr_width, r3, g3, b3);

                        double tile_portion = 0.7;
                        double texture_portion = 1.0 - tile_portion;

                        Color t11(TextBuffer[i][j].r*tile_portion + (double(r3) / 255.0)*texture_portion,TextBuffer[i][j].g*tile_portion + (double(g3) / 255.0)*texture_portion,TextBuffer[i][j].b*tile_portion + (double(b3) / 255.0)*texture_portion);
                        //Color t11(1.0,1.0,1.0);
                        TextBuffer[i][j]=t11;
                    }
                    //TextBuffer[i][j]=checker_color;
                }
                else if(IntersectingObj==2)//Sphere
                {
                    struct point intrsc=Ray_Sphere_Intersect(eye_to_pixel,vector_of_Sphere[IntersectingIndex]);
                    TextBuffer[i][j]=Am_Diff_Spec_Reflec(eye_to_pixel,intrsc,1,2,IntersectingIndex);
                    //TextBuffer[i][j]=vector_of_Sphere[IntersectingIndex].Sp_color;
//                    TextBuffer[i][j].r=1;
//                    TextBuffer[i][j].g=1;
//                    TextBuffer[i][j].b=0;
                }
                else if(IntersectingObj==3)//Triangle
                {
                    struct point intrsc=Ray_Triangle_Intersect(eye_to_pixel,vector_of_Triangle[IntersectingIndex]);
                    TextBuffer[i][j]=Am_Diff_Spec_Reflec(eye_to_pixel,intrsc,1,3,IntersectingObj);
                    //TextBuffer[i][j]=vector_of_Triangle[IntersectingIndex].TriColor;
//                    TextBuffer[i][j].r=1;
//                    TextBuffer[i][j].g=1;
//                    TextBuffer[i][j].b=0;
                }
            }
            else{
                TextBuffer[i][j]=blk;
            }
        }

    }
    cout<<"Kharan"<<endl;

    bitmap_image img(num_pxl,num_pxl);
    for(int i=0;i<num_pxl;i++)
    {
        for(int j=0;j<num_pxl;j++)
        {
            img.set_pixel(i,j,(TextBuffer[i][j].r)*255.0,(TextBuffer[i][j].g)*255.0,(TextBuffer[i][j].b)*255.0);
        }
    }
    img.save_image("outputt.bmp");
    flag_space=0;
    cout<<"Hoia Gese"<<endl;
}

void ReSize()
{
    int width=Given_Photo.width();
    int height=Given_Photo.height();
    bitmap_image noya((int)ckr_width,(int)ckr_width);
    bitmap_image noya1((int)ckr_width,(int)ckr_width);
    bitmap_image noya3((int)ckr_width,(int)ckr_width);
    bitmap_image noya2((int)ckr_width,(int)ckr_width);

    Given_Photo.subsample(noya);
    width=noya.width();
    height=noya.height();

    noya.subsample(noya1);
    noya1.subsample(noya2);

    width=noya2.width();
    height=noya2.height();
    cout<<"noya "<<width<<" "<<height<<endl;

    noya2.save_image("NoyaPic.bmp");
}



void ReadFile()
{
    ifstream inputFile;
    inputFile.open("description.txt");

    inputFile>>nearDis>>farDias>>fovY_axis>>aspect_ro;
    inputFile>>level_rec;
    inputFile>>num_pxl;

    inputFile>>ckr_width;
    inputFile>>ckr_am>>ckr_diff>>ckr_ref;
    //Checker board
    Checker_Brd ckrr2(ckr_width,ckr_am,ckr_diff,ckr_ref);
    ckrr1=ckrr2;
    vector_of_Chkr.push_back(ckrr2);

    inputFile>>num_of_objs;
    cout<<num_of_objs<<endl;

    string cmd;
    for(int i=0;i<num_of_objs;i++)
    {
        inputFile>>cmd;
        if(cmd=="sphere")
        {
            inputFile>>cenX>>cenY>>cenZ;
            inputFile>>radiuss;
            inputFile>>colrX>>colrY>>colrZ;
            cout<<colrX<<" "<<colrY<<" "<<colrZ<<endl;
            inputFile>>am>>diff>>spec>>refl;
            inputFile>>shini;
            cout<<shini<<endl;

            SpHere Spr(cenX,cenY,cenZ, radiuss, colrX,colrY,colrZ, am, diff, spec,refl,shini);
            //glColor3f(colrX,colrY,colrZ);
            vector_of_Sphere.push_back(Spr);
        }
        else if(cmd=="pyramid")
        {
            inputFile>>low_posX>>low_posY>>low_posZ;
            inputFile>>wdt>>hgt;
            inputFile>>colrX>>colrY>>colrZ;
            inputFile>>am>>diff>>spec>>refl;
            inputFile>>shini;

            struct point low_pnt={low_posX, low_posY, low_posZ};
            struct point pnt2={low_posX+1*wdt, low_posY, low_posZ};
            struct point pnt3={low_posX+1*wdt, low_posY+1*wdt, low_posZ};
            struct point pnt4={low_posX, low_posY+1*wdt, low_posZ};
            struct point HEIGHT={low_posX+0.5*wdt, low_posY+0.5*wdt, low_posZ+1*hgt};

            Color colr11(colrX,colrY,colrZ);

            TriAngle tr(HEIGHT,low_pnt,pnt2,colr11,am,diff,spec,refl,shini);
            vector_of_Triangle.push_back(tr);

            TriAngle trr(HEIGHT,pnt2,pnt3,colr11,am,diff,spec,refl,shini);
            vector_of_Triangle.push_back(trr);

            TriAngle tr2(HEIGHT,pnt3,pnt4,colr11,am,diff,spec,refl,shini);
            vector_of_Triangle.push_back(tr2);

            TriAngle tr3(HEIGHT,pnt4,low_pnt,colr11,am,diff,spec,refl,shini);
            vector_of_Triangle.push_back(tr3);

            Pyramid pyra2(low_posX,low_posY,low_posZ,wdt,hgt,colrX,colrY,colrZ,am,diff,spec,refl,shini);
            pyra1=pyra2;
            vector_of_pyra.push_back(pyra2);
        }
    }
    inputFile>>num_of_light;
    for(int i=0;i<num_of_light;i++)
    {
        inputFile>>light_srcX>>light_srcY>>light_srcZ>>light_fall_off;
        light_src D(light_srcX,light_srcY,light_srcZ,light_fall_off);
        vector_of_light.push_back(D);
    }

    inputFile>>num_of_sptLight;
    for(int i=0;i<num_of_sptLight;i++)
    {
        inputFile>>spt_srcX>>spt_srcY>>spt_srcZ>>spt_fall_off;
        inputFile>>lk_X>>lk_Y>>lk_Z;
        inputFile>>cutt_off;
        spot_light_src E(spt_srcX,spt_srcY,spt_srcZ,spt_fall_off,lk_X,lk_Y,lk_Z,cutt_off);
        vector_of_spot_light.push_back(E);
        cout<<cutt_off<<endl;
    }
    inputFile.close();

}







void drawSphere(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
			}glEnd();
		}
	}
}



void keyboardListener(unsigned char key, int x,int y){
	switch(key){
	    case '0':
            Building_Image();
            break;

		case '1':
		    r=pnt_rotationn(u,r,-2*2*pi/360.0);
			l=pnt_rotationn(u,l,-2*2*pi/360.0);
			break;
        case '2':
			r=pnt_rotationn(u,r,2*2*pi/360.0);
			l=pnt_rotationn(u,l,2*2*pi/360.0);
			break;

        case '3':
			u=pnt_rotationn(r,u,2*2*pi/360.0);
			l=pnt_rotationn(r,l,2*2*pi/360.0);
			break;
        case '4':
			u=pnt_rotationn(r,u,-2*2*pi/360.0);
			l=pnt_rotationn(r,l,-2*2*pi/360.0);
			break;

        case '5':
			r=pnt_rotationn(l,r,-2*2*pi/360.0);
			u=pnt_rotationn(l,u,-2*2*pi/360.0);
			break;
        case '6':
			r=pnt_rotationn(l,r,2*2*pi/360.0);
			u=pnt_rotationn(l,u,2*2*pi/360.0);
			break;
        case '7':
            for(int i=0;i<vector_of_Triangle.size();i++)
            {
                cout<<"Triangle :"<<vector_of_Triangle[i].T1.x<<endl;
            }
            break;
        case ' ':
            //cout<<"bched"<<endl;
            flag_space=1;

            break;

		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			//cameraHeight -= 3.0;
			eye_pos=pnt_subtract(eye_pos, l);
			break;
		case GLUT_KEY_UP:		// up arrow key
			eye_pos=pnt_addition(eye_pos,l);
			break;

		case GLUT_KEY_RIGHT:
			eye_pos=pnt_addition(eye_pos,r);
			break;
		case GLUT_KEY_LEFT:
			eye_pos=pnt_subtract(eye_pos,r);
			break;

		case GLUT_KEY_PAGE_UP:
		    eye_pos=pnt_addition(eye_pos,u);
			break;
		case GLUT_KEY_PAGE_DOWN:
		    eye_pos=pnt_subtract(eye_pos,u);
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP

			}
			break;

		case GLUT_RIGHT_BUTTON:
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}



void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	//gluLookAt(0,0,200,	0,0,0,	0,1,0);
	gluLookAt(eye_pos.x,eye_pos.y,eye_pos.z, eye_pos.x+l.x,eye_pos.y+l.y,eye_pos.z+l.z, u.x,u.y,u.z);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects


	//ckrr1.Checker_Draw();
    //pyra1.py_Draw();
    for(int i=0;i<vector_of_Chkr.size();i++)
    {
        vector_of_Chkr[i].Checker_Draw();
    }
    for(int i=0;i<vector_of_pyra.size();i++)
    {
        vector_of_pyra[i].py_Draw();
    }
    for(int i=0;i<vector_of_Sphere.size();i++)
    {
        vector_of_Sphere[i].Sp_Draw();
    }
    for(int i=0;i<vector_of_light.size();i++)
    {
        vector_of_light[i].lt_draw();
    }
    for(int i=0;i<vector_of_spot_light.size();i++)
    {
        vector_of_spot_light[i].spt_draw();
    }
    //cout<<"Triangle "<<vector_of_Triangle.size()<<endl;





	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){

	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization


	ReadFile();

	eye_pos = {100, 140, 30};
    l = { -1.0/sqrt(2.0), -1.0/sqrt(2.0), 0.0 };
	r = { -1.0/sqrt(2.0), 1.0/sqrt(2.0), 0.0 };
	u = { 0.0, 0.0, 1.0 };

	Given_Photo=bitmap_image("texture.bmp");
	ReSize();
	Image2=bitmap_image("NoyaPic.bmp");

	//ckrr1.Checker_Draw();

	//clear the screen
	glClearColor(0,0,0,0);


	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	gluPerspective(fovY_axis,	aspect_ro,	nearDis,	farDias);

}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	ifstream fl;
	fl.open("description.txt");
	int screenXY;
	fl>>screenXY>>screenXY>>screenXY>>screenXY>>screenXY>>screenXY;
	fl.close();
	glutInitWindowSize(screenXY, screenXY);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
