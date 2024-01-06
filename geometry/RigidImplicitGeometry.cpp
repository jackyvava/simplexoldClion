//
// Created by zhangjie on 2024/1/6.
//

#include "RigidImplicitGeometry.h"
#include "LevelSet.h"
#include "MarchingCubes.h"
#include "Constants.h"

#include <vector>
#include <iostream>

//template class RigidImplicitGeometry<2>;
//template class RigidImplicitGeometry<3>;
//template class SphereGeometry<2>;
//template class SphereGeometry<3>;
//template class BoxGeometry<2>;
//template class BoxGeometry<3>;

void RigidImplicitGeometry<2>::Precompute_Mesh(double dx){
    Veci<2> cell_counts;
    cell_counts(0) = int(2*(bbox_max(0)-bbox_min(0))/dx);
    cell_counts(1) = int(2*(bbox_max(1)-bbox_min(1))/dx);
    grid.Initialize(cell_counts, dx, -Vec<2>(cell_counts(0)*dx/2, cell_counts(1)*dx/2));

    LevelSet<2> levelset(grid);
    iterate_cell_d(iter, grid, 2)
    {
        Veci<2> cell = iter.Coord();
        int index = grid.Cell_Index(cell);
        Vec<2> pos = grid.Center(cell);
        levelset.phi(cell) = Get_Phi(pos);
    }

    //levelset.Fast_Marching();
    MarchingCubes<2> marchincubes(levelset);
    marchincubes.Marching();
    mesh = *marchincubes.mesh;//deep copy
}

void RigidImplicitGeometry<3>::Precompute_Mesh(double dx){
    Veci<3> cell_counts;
    cell_counts(0) = int(2*(bbox_max(0)-bbox_min(0))/dx);
    cell_counts(1) = int(2*(bbox_max(1)-bbox_min(1))/dx);
    cell_counts(2) = int(2*(bbox_max(2)-bbox_min(2))/dx);
    grid.Initialize(cell_counts, dx, -Vec<3>(cell_counts(0)*dx/2, cell_counts(1)*dx/2, cell_counts(2)*dx/2));

    LevelSet<3> levelset(grid);
    iterate_cell_d(iter, grid, 3)
    {
        Veci<3> cell = iter.Coord();
        int index = grid.Cell_Index(cell);
        Vec<3> pos = grid.Center(cell);
        levelset.phi(cell) = Get_Phi(pos);
    }
    //levelset.Fast_Marching();
    MarchingCubes<3> marchincubes(levelset);
    marchincubes.Marching();
    mesh = *marchincubes.mesh;//deep copy
}

//Sphere
SphereGeometry<2>::SphereGeometry(const double _r){r = _r;bbox_min=Vec<2>(-r,-r);bbox_max=Vec<2>(r,r);}
void SphereGeometry<2>::Initialize(const double _r){r = _r;bbox_min=Vec<2>(-r,-r);bbox_max=Vec<2>(r,r);}
SphereGeometry<2>& SphereGeometry<2>::operator=(const SphereGeometry& copy){r = copy.r;bbox_min=Vec<2>(-r,-r);bbox_max=Vec<2>(r,r);return *this;}
SphereGeometry<2>::SphereGeometry(const SphereGeometry& copy){r = copy.r;bbox_min=Vec<2>(-r,-r);bbox_max=Vec<2>(r,r);}
double SphereGeometry<2>::Get_Phi(const Vec<2>& pos){return pos.norm() - r;}
double SphereGeometry<2>::Get_Mass(const double density){return 3.1415926 * r * r * density;}
double SphereGeometry<2>::Get_Inertia_Body(double density){return 0.5 * Get_Mass(density) * r * r;}
void SphereGeometry<2>::Project_Out(Vec<2>& pos)
{
    if (pos.norm() < 1e-10)
    {
        pos = Vec<2>::Unit(1) * r;
    }
    else
    {
        pos.normalize();
        pos = pos * r;
    }
}
void SphereGeometry<2>::Project_On(Vec<2>& pos){
    pos.normalize();
    pos = pos*r;
}
void SphereGeometry<2>::Project_Out_Directional(Vec<2>& pos, const Vec<2>& direction){
    double dd = direction.dot(direction);
    double pp = pos.dot(pos);
    double pd = pos.dot(direction);
    double delta = (-pd+sqrt(pd*pd+dd*(r*r-pp)))/dd;
    pos += delta*direction;
}
Vec<2> SphereGeometry<2>::Get_Normal(const Vec<2>& pos){
    if(pos.norm()<1e-10) return Vec<2>::Unit(1);
    else return pos/pos.norm();
}

SphereGeometry<3>::SphereGeometry(const double _r){r = _r;bbox_min=Vec<3>(-r,-r,-r);bbox_max=Vec<3>(r,r,r);}
void SphereGeometry<3>::Initialize(const double _r){r = _r;bbox_min=Vec<3>(-r,-r,-r);bbox_max=Vec<3>(r,r,r);}
SphereGeometry<3>& SphereGeometry<3>::operator=(const SphereGeometry& copy){r = copy.r;bbox_min=Vec<3>(-r,-r,-r);bbox_max=Vec<3>(r,r,r);return *this;}
double SphereGeometry<3>::Get_Phi(const Vec<3>& pos){return pos.norm() - r;}
double SphereGeometry<3>::Get_Mass(const double density){return 4.0 / 3 * pi * r * r * r * density;}
Mat<3> SphereGeometry<3>::Get_Inertia_Body(const double density){return Mat<3>::Identity() * 0.4 * Get_Mass(density) * r * r;}
void SphereGeometry<3>::Project_Out(Vec<3>& pos)
{
    if (pos.norm() < 1e-10)
    {
        pos = Vec<3>::Unit(1) * r;
    }
    else
    {
        pos.normalize();
        pos = pos * r;
    }
}
void SphereGeometry<3>::Project_On(Vec<3>& pos){
    pos.normalize();
    pos = pos*r;
}
void SphereGeometry<3>::Project_Out_Directional(Vec<3>& pos, const Vec<3>& direction){
    double dd = direction.dot(direction);
    double pp = pos.dot(pos);
    double pd = pos.dot(direction);
    double delta = (-pd+sqrt(pd*pd+dd*(r*r-pp)))/dd;
    pos += delta*direction;
}
Vec<3> SphereGeometry<3>::Get_Normal(const Vec<3>& pos){
    if(pos.norm()<1e-10)return Vec<3>::Unit(1);
    else return pos/pos.norm();
}


//Box
BoxGeometry<2>::BoxGeometry(const double _a, const double _b){a = _a;b = _b;bbox_min=Vec<2>(-a,-b);bbox_max=Vec<2>(a,b);}
void BoxGeometry<2>::Initialize(const double _a, const double _b){a = _a;b = _b;bbox_min=Vec<2>(-a,-b);bbox_max=Vec<2>(a,b);}
BoxGeometry<2>& BoxGeometry<2>::operator=(const BoxGeometry<2>& copy){a = copy.a;b = copy.b;bbox_min=Vec<2>(-a,-b);bbox_max=Vec<2>(a,b);return *this;}
double BoxGeometry<2>::Get_Mass(const double density){return a * b * density;}
double BoxGeometry<2>::Get_Inertia_Body(const double density){return 1.0 / 12 * Get_Mass(density) * (a * a + b * b);}
double BoxGeometry<2>::Get_Phi(const Vec<2>& pos)
{
    std::vector<double> positives;
    double d[4];
    positives.clear();
    d[0] = -0.5 * a - pos[0];
    d[1] = pos[0] - 0.5 * a;
    d[2] = -0.5 * b - pos[1];
    d[3] = pos[1] - 0.5 * b;
    for (int i = 0; i < 4; i++)
        if (d[i] > 0)
            positives.push_back(d[i]);
    if (positives.size() == 0)
    {
        double max = d[0];
        for (int i = 1; i < 4; i++)
            if (max < d[i])
                max = d[i];
        return max;
    }
    else if (positives.size() == 1)
        return positives[0];
    else if (positives.size() == 2)
        return std::sqrt(positives[0] * positives[0] + positives[1] * positives[1]);
    else
    {
        std::cerr << "wrong box phi" << std::endl;
        return 0;
    }
}
void BoxGeometry<2>::Project_Out(Vec<2>& pos)
{
    int axis = 0;
    double min_dis = std::min(0.5 * a - pos[0], pos[0] + 0.5 * a);
    double y_dis = std::min(0.5 * b - pos[1], pos[1] + 0.5 * b);
    if (y_dis < min_dis)
    {
        axis = 1;
        min_dis = y_dis;
    }
    if (axis == 0)
    {
        if (pos[0] > 0)
            pos[0] = 0.5 * a;
        else
            pos[0] = -0.5 * a;
    }
    else
    {
        if (pos[1] > 0)
            pos[1] = 0.5 * b;
        else
            pos[1] = -0.5 * b;
    }
}
void BoxGeometry<2>::Project_On(Vec<2>& pos){
    int axis = 0;
    if(pos(0)>a/2) pos(0) = a/2;
    else if(pos(0)<-a/2) pos(0) = -a/2;
    if(pos(1)>b/2) pos(1) = b/2;
    else if(pos(1)<-b/2) pos(1) = -b/2;
}
void BoxGeometry<2>::Project_Out_Directional(Vec<2>& pos, const Vec<2>& direction){
    auto sgn = [](double x) -> double {return double((x>0.0)-(x<0.0));};
    if(std::abs(direction(0))<1e-6) {pos(1) = b/2*sgn(direction(1)); return;}
    if(std::abs(direction(1))<1e-6) {pos(0) = a/2*sgn(direction(0)); return;}

    double ty = pos(1) + (sgn(direction(0))*a/2-pos(0))*direction(1)/direction(0);
    if((ty>=-b/2) && (ty<=b/2)) {pos = Vec<2>(sgn(direction(0))*a/2, ty); return;}
    double tx = pos(0) + (sgn(direction(1))*b/2-pos(1))*direction(0)/direction(1);
    if((tx>=-a/2) && (tx<=a/2)) {pos = Vec<2>(tx, sgn(direction(1))*b/2); return;}
    return;
}
Vec<2> BoxGeometry<2>::Get_Normal(const Vec<2>& pos){
    Vec<2> ap = Vec<2>(abs(pos(0)), abs(pos(1)));
    Vec<2> n = Vec<2>::Zero();
    if(ap(0)<a/2){
        if(ap(1)<b/2){
            if((a/2-ap(0))<(b/2-ap(1))) n(1) = 1.0;
            else n(0) = 1.0;
        }
        else n(1) = 1.0;
    }
    else if(ap(1)<b/2) n(0) = 1.0;
    else n = Vec<2>(ap-Vec<2>(a/2,b/2));
    if(pos(0)<0) n(0) = -n(0);
    if(pos(1)<0) n(1) = -n(1);
    return n.normalized();
}

BoxGeometry<3>::BoxGeometry(const double _a, const double _b, const double _c)
{
    a = _a;
    b = _b;
    c = _c;
    bbox_min=Vec<3>(-a,-b,-c);
    bbox_max=Vec<3>(a,b,c);
}
void BoxGeometry<3>::Initialize(const double _a, const double _b, const double _c)
{
    a = _a;
    b = _b;
    c = _c;
    bbox_min=Vec<3>(-a,-b,-c);
    bbox_max=Vec<3>(a,b,c);
}
BoxGeometry<3>& BoxGeometry<3>::operator=(const BoxGeometry<3>& copy)
{
    a = copy.a;
    b = copy.b;
    c = copy.c;
    bbox_min=Vec<3>(-a,-b,-c);
    bbox_max=Vec<3>(a,b,c);
    return *this;
}
double BoxGeometry<3>::Get_Phi(const Vec<3>& pos)
{
    std::vector<double> positives;
    double d[6];
    positives.clear();
    d[0] = -0.5 * a - pos[0];
    d[1] = pos[0] - 0.5 * a;
    d[2] = -0.5 * b - pos[1];
    d[3] = pos[1] - 0.5 * b;
    d[4] = -0.5 * c - pos[2];
    d[5] = pos[2] - 0.5 * c;
    for (int i = 0; i < 6; i++)
        if (d[i] > 0)
            positives.push_back(d[i]);
    if (positives.size() == 0)
    {
        double max = d[0];
        for (int i = 1; i < 6; i++)
            if (max < d[i])
                max = d[i];
        return max;
    }
    else if (positives.size() == 1)
        return positives[0];
    else if (positives.size() == 2)
        return std::sqrt(positives[0] * positives[0] + positives[1] * positives[1]);
    else if (positives.size() == 3)
        return std::sqrt(positives[0] * positives[0] + positives[1] * positives[1] + positives[2] * positives[2]);
    else
    {
        std::cerr << "wrong box phi" << std::endl;
        return 0;
    }
}
double BoxGeometry<3>::Get_Mass(const double density){return a * b * c * density;}
Mat<3> BoxGeometry<3>::Get_Inertia_Body(const double density)
{
    Mat<3> j = Mat<3>::Zero();
    double m = Get_Mass(density);
    j(0, 0) = 1.0 / 12 * m * (b * b + c * c);
    j(1, 1) = 1.0 / 12 * m * (a * a + c * c);
    j(2, 2) = 1.0 / 12 * m * (a * a + b * b);
    return j;
}
void BoxGeometry<3>::Project_Out(Vec<3>& pos)
{
    int axis = 0;
    double min_dis = std::min(0.5 * a - pos[0], pos[0] + 0.5 * a);
    double y_dis = std::min(0.5 * b - pos[1], pos[1] + 0.5 * b);
    if (y_dis < min_dis)
    {
        axis = 1;
        min_dis = y_dis;
    }
    double z_dis = std::min(0.5 * c - pos[2], pos[2] + 0.5 * c);
    if (z_dis < min_dis)
    {
        axis = 2;
        min_dis = z_dis;
    }
    if (axis == 0)
    {
        if (pos[0] > 0)
            pos[0] = 0.5 * a;
        else
            pos[0] = -0.5 * a;
    }
    else if (axis == 1)
    {
        if (pos[1] > 0)
            pos[1] = 0.5 * b;
        else
            pos[1] = -0.5 * b;
    }
    else
    {
        if (pos[2] > 0)
            pos[2] = 0.5 * c;
        else
            pos[2] = -0.5 * c;
    }
}
void BoxGeometry<3>::Project_On(Vec<3>& pos){
    if(pos(0)>a/2) pos(0)=a/2;
    else if(pos(0)<-a/2) pos(0)=-a/2;
    if(pos(1)>b/2) pos(1)=b/2;
    else if(pos(1)<-b/2) pos(1)=-b/2;
    if(pos(2)>c/2) pos(2)=c/2;
    else if(pos(2)<-c/2) pos(2)=-c/2;
}
void BoxGeometry<3>::Project_Out_Directional(Vec<3>& pos, const Vec<3>& direction){
    auto sgn = [](double x) -> double {return double((x>0.0)-(x<0.0));};
    double sa = a/2*sgn(direction(0));
    double sb = b/2*sgn(direction(1));
    double sc = c/2*sgn(direction(2));
    double tx, ty, tz;
    if(std::abs(direction(0))<1e-6){
        if(std::abs(direction(1))<1e-6) {pos(2) = sc; return;}
        if(std::abs(direction(2))<1e-6) {pos(1) = sb; return;}
        ty = pos(1) + (sc-pos(2))*direction(1)/direction(2);
        if((ty>=-b/2) && (ty<=b/2)) {pos = Vec<3>(pos(0), ty, sc); return;}
        tz = pos(2) + (sb-pos(1))*direction(2)/direction(1);
        if((tz>=-c/2) && (tz<=c/2)) {pos = Vec<3>(pos(0), sb, tz); return;}
    }
    else if(std::abs(direction(1))<1e-6){
        if(std::abs(direction(2))<1e-6) {pos(0) = sa; return;}
        tx = pos(0) + (sc-pos(2))*direction(0)/direction(2);
        if((tx>=-a/2) && (tx<=a/2)) {pos = Vec<3>(tx, pos(1), sc); return;}
        tz = pos(2) + (sa-pos(0))*direction(2)/direction(0);
        if((tz>=-c/2) && (tz<=c/2)) {pos = Vec<3>(sa, pos(1), tz); return;}
    }
    else if(std::abs(direction(2))<1e-6){
        tx = pos(0) + (sb-pos(1))*direction(0)/direction(1);
        if((tx>=-a/2) && (tx<=a/2)) {pos = Vec<3>(tx, sb, pos(2)); return;}
        ty = pos(1) + (sa-pos(0))*direction(1)/direction(0);
        if((ty>=-b/2) && (ty<=b/2)) {pos = Vec<3>(sa, ty, pos(2)); return;}
    }

    ty = pos(1) + (sa-pos(0))*direction(1)/direction(0);
    tz = pos(2) + (sa-pos(0))*direction(2)/direction(0);
    if((std::abs(ty)<=b/2) && (std::abs(tz)<=c/2)) {pos = Vec<3>(sa, ty, tz); return;}
    tx = pos(0) + (sb-pos(1))*direction(0)/direction(1);
    tz = pos(2) + (sb-pos(1))*direction(2)/direction(1);
    if((std::abs(tx)<=a/2) && (std::abs(tz)<=c/2)) {pos = Vec<3>(tx, sb, tz); return;}
    tx = pos(0) + (sc-pos(2))*direction(0)/direction(2);
    ty = pos(1) + (sc-pos(2))*direction(1)/direction(2);
    if((std::abs(tx)<=a/2) && (std::abs(ty)<=c/2)) {pos = Vec<3>(tx, ty, sc); return;}
    return;
}
Vec<3> BoxGeometry<3>::Get_Normal(const Vec<3>& pos){
    Vec<3> normal; for(int i=0;i<3;i++)
        normal[i]=(Get_Phi(pos+1e-2*Vec<3>::Unit(i)*grid.dx)-Get_Phi(pos-1e-2*Vec<3>::Unit(i)*grid.dx))/(2*1e-2*grid.dx);
    return normal.normalized();
}


//Cylinder
CylinderGeometry<2>::CylinderGeometry(double _r, double _h){r=_r;h=_h;bbox_min=Vec<2>(-r,-h/2);bbox_max=Vec<2>(r,h/2);}
void CylinderGeometry<2>::Initialize(double _r, double _h){r=_r;h=_h;bbox_min=Vec<2>(-r,-h/2);bbox_max=Vec<2>(r,h/2);}
CylinderGeometry<2>& CylinderGeometry<2>::operator=(const CylinderGeometry<2>& copy){r=copy.r;h=copy.h;bbox_min=Vec<2>(-r, -h/2);bbox_max=Vec<2>(r,h/2);return *this;}
double CylinderGeometry<2>::Get_Mass(double density){return 2*r*h*density;}
double CylinderGeometry<2>::Get_Inertia_Body(double density){return 1.0/12*Get_Mass(density)*(4*r*r + h*h);}
double CylinderGeometry<2>::Get_Phi(const Vec<2>& pos)
{
    std::vector<double> positives;
    double d[4];
    positives.clear();
    d[0] = - r - pos[0];
    d[1] = pos[0] - r;
    d[2] = -0.5 * h - pos[1];
    d[3] = pos[1] - 0.5 * h;
    for (int i = 0; i < 4; i++)
        if (d[i] > 0)
            positives.push_back(d[i]);
    if (positives.size() == 0)
    {
        double max = d[0];
        for (int i = 1; i < 4; i++)
            if (max < d[i])
                max = d[i];
        return max;
    }
    else if (positives.size() == 1)
        return positives[0];
    else if (positives.size() == 2)
        return std::sqrt(positives[0] * positives[0] + positives[1] * positives[1]);
    else
    {
        std::cerr << "wrong cylinder 2D phi" << std::endl;
        return 0;
    }
}
void CylinderGeometry<2>::Project_Out(Vec<2>& pos){
    int axis = 0;
    double min_dis = std::min(r - pos[0], pos[0] + r);
    double y_dis = std::min(0.5 * h - pos[1], pos[1] + 0.5 * h);
    if (y_dis < min_dis)
    {
        axis = 1;
        min_dis = y_dis;
    }
    if (axis == 0)
    {
        if (pos[0] > 0)
            pos[0] = r;
        else
            pos[0] = -r;
    }
    else
    {
        if (pos[1] > 0)
            pos[1] = 0.5*h;
        else
            pos[1] = -0.5*h;
    }
}
void CylinderGeometry<2>::Project_On(Vec<2>& pos){
    if(pos(0)<-r) pos(0) = -r;
    else if(pos(0)>r) pos(0) = r;
    if(pos(1)<-h/2) pos(1) = -h/2;
    else if(pos(1)>h/2) pos(1) = h/2;
}
void CylinderGeometry<2>::Project_Out_Directional(Vec<2>& pos, const Vec<2>& direction){
    auto sgn = [](double x) -> double {return double((x>0.0)-(x<0.0));};
    if(std::abs(direction(0))<1e-6) {pos(1) = h/2*sgn(direction(1)); return;}
    if(std::abs(direction(1))<1e-6) {pos(0) = r*sgn(direction(0)); return;}

    double ty = pos(1) + (sgn(direction(0))*r-pos(0))*direction(1)/direction(0);
    if((ty>=-h/2) && (ty<=h/2)) {pos = Vec<2>(sgn(direction(0))*r, ty); return;}
    double tx = pos(0) + (sgn(direction(1))*h/2-pos(1))*direction(0)/direction(1);
    if((tx>=-r) && (tx<=r)) {pos = Vec<2>(tx, sgn(direction(1))*h/2); return;}
    return;
}
Vec<2> CylinderGeometry<2>::Get_Normal(const Vec<2>& pos){
    Vec<2> ap = Vec<2>(abs(pos(0)), abs(pos(1)));
    Vec<2> n = Vec<2>::Zero();
    if(ap(0)<r){
        if(ap(1)<h/2){
            if((r-ap(0))<(h/2-ap(1))) n(1) = 1.0;
            else n(0) = 1.0;
        }
        else n(1) = 1.0;
    }
    else if(ap(1)<h/2) n(0) = 1.0;
    else n = Vec<2>(ap-Vec<2>(r,h/2));
    if(pos(0)<0) n(0) = -n(0);
    if(pos(1)<0) n(1) = -n(1);
    return n.normalized();
}

CylinderGeometry<3>::CylinderGeometry(double _r, double _h){r=_r;h=_h;bbox_min=Vec<3>(-r,-r,-h/2);bbox_max=Vec<3>(r,r,h/2);}
void CylinderGeometry<3>::Initialize(double _r, double _h){r=_r;h=_h;bbox_min=Vec<3>(-r,-r,-h/2);bbox_max=Vec<3>(r,r,h/2);}
CylinderGeometry<3>& CylinderGeometry<3>::operator=(const CylinderGeometry<3>& copy){r=copy.r;h=copy.h;bbox_min=Vec<3>(-r,-r,-h/2);bbox_max=Vec<3>(r,r,h/2);return *this;}
double CylinderGeometry<3>::Get_Mass(double density){return pi*r*r*h*density;}
Mat<3> CylinderGeometry<3>::Get_Inertia_Body(double density){
    double mass = Get_Mass(density);
    Mat<3> tI = Mat<3>::Zero();
    tI(0,0) = 1.0/12.0*mass*h*h + 1.0/4.0*mass*r*r;
    tI(1,1) = tI(0,0);
    tI(2,2) = 1.0/2.0*mass*r*r;
    return tI;
}
double CylinderGeometry<3>::Get_Phi(const Vec<3>& pos){
    std::vector<double> positives;
    double d[4];
    double xhat = std::sqrt(pos(0)*pos(0)+pos(1)*pos(1));
    positives.clear();
    d[0] = - r - xhat;
    d[1] = xhat - r;
    d[2] = -0.5 * h - pos(2);
    d[3] = pos(2) - 0.5 * h;
    for (int i = 0; i < 4; i++)
        if (d[i] > 0)
            positives.push_back(d[i]);
    if (positives.size() == 0)
    {
        double max = d[0];
        for (int i = 1; i < 4; i++)
            if (max < d[i])
                max = d[i];
        return max;
    }
    else if (positives.size() == 1)
        return positives[0];
    else if (positives.size() == 2)
        return std::sqrt(positives[0] * positives[0] + positives[1] * positives[1]);
    else
    {
        std::cerr << "wrong cylinder 3D phi" << std::endl;
        return 0;
    }
}
void CylinderGeometry<3>::Project_Out(Vec<3>& pos){
    int axis = 0;
    double xhat = std::sqrt(pos(0)*pos(0)+pos(1)*pos(1));
    double x_dis = std::min(r-xhat, xhat+r);
    double z_dis = std::min(0.5*h-pos(2),pos(2)+0.5*h);
    if(z_dis < x_dis){axis = 1;}
    if(axis==0){
        pos(0) = pos(0)/xhat*r;
        pos(1) = pos(1)/xhat*r;
    }
    else{
        if(pos(2)>0)pos(2)=0.5*h;
        else pos(2)=-0.5*h;
    }
}
void CylinderGeometry<3>::Project_On(Vec<3>& pos){
    double pp = pos(0)*pos(0) + pos(1)*pos(1);
    double k = r/sqrt(pp); pos(0) = k*pos(0); pos(1) = k*pos(1);
    if(pos(2)>h/2) pos(2) = h/2;
    else if(pos(2)<-h/2) pos(2) = -h/2;
}
void CylinderGeometry<3>::Project_Out_Directional(Vec<3>& pos, const Vec<3>& direction){
    auto sgn = [](double x) -> double {return double((x>0.0)-(x<0.0));};
    double dd = direction(0)*direction(0)+direction(1)*direction(1);
    if(dd<1e-8) {pos(2) = sgn(direction(2)*h/2); return;}
    double pp = pos(0)*pos(0)+pos(1)*pos(1);
    double pd = pos(0)*direction(0)+pos(1)*direction(1);
    double delta = (-pd+sqrt(pd*pd+dd*(r*r-pp)))/dd;
    double tz = pos(2) + direction(2)/sqrt(dd)*delta;
    if((tz>=-h/2) && (tz<=h/2)) {pos += delta*direction; return;}
    delta = (sgn(direction(2))*h/2-pos(2))/direction(2);
    pos += delta*direction; return;
}
Vec<3> CylinderGeometry<3>::Get_Normal(const Vec<3>& pos){
    Vec<3> normal; for(int i=0;i<3;i++)
        normal[i]=(Get_Phi(pos+1e-2*Vec<3>::Unit(i)*grid.dx)-Get_Phi(pos-1e-2*Vec<3>::Unit(i)*grid.dx))/(2*1e-2*grid.dx);
    return normal.normalized();

}
