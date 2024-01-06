//
// Created by zhangjie on 2024/1/6.
//

#pragma once

#include "RigidGeometry.h"
#include "Grid.h"
#include "LevelSet.h"
#include "tiny_obj_loader/TinyObjLoader.h"
#include "GH.h"

#include<iostream>
#include<limits>
#include<Eigen/Geometry>

//template class RigidGeometry<2>;
//template class RigidGeometry<3>;
//template class RigidExplicitGeometry<2>;
//template class RigidExplicitGeometry<3>;

void RigidGeometry<2>::Get_Mesh(const double orientation, const Vec<2>& position, SurfaceMesh<2>& tmesh){
    tmesh.vertices = std::make_shared<Array<Vec<2> > >();
    tmesh.elements = mesh.elements;
    for(Vec<2>& pos : *(mesh.vertices)){
        tmesh.vertices->push_back(Eigen::Rotation2Dd(orientation)*pos + position);
    }
}

void RigidGeometry<2>::Get_Mesh(const Mat<2>& orientation, const Vec<2>& position, SurfaceMesh<2>& tmesh){
    tmesh.vertices = std::make_shared<Array<Vec<2> > >();
    tmesh.elements = mesh.elements;
    for(Vec<2>& pos : *(mesh.vertices)){
        tmesh.vertices->push_back(orientation*pos + position);
    }
}

void RigidGeometry<3>::Get_Mesh(const Mat<3>& orientation, const Vec<3>& position, SurfaceMesh<3>& tmesh){
    tmesh.vertices = std::make_shared<Array<Vec<3> > >();
    tmesh.elements = mesh.elements;
    for(Vec<3>& pos : *(mesh.vertices)){
        tmesh.vertices->push_back(orientation*pos + position);
    }
}

bool RigidGeometry<2>::Is_Inside(const Vec<2>& pos){
    return (Get_Phi(pos)<0);
}

bool RigidGeometry<3>::Is_Inside(const Vec<3>& pos){
    return (Get_Phi(pos)<0);
}

//Explicit geometry

RigidExplicitGeometry<2>::RigidExplicitGeometry(const std::string& file, const Vec<2>& scale){
    //it seems tiny obj loader doesn't support 2D obj
    std::cerr << "2D obj loader not implemented" << std::endl;
}

RigidExplicitGeometry<3>::RigidExplicitGeometry(const std::string& file, const Vec<3>& scale){
    std::shared_ptr<SurfaceMesh<3> > mesh_ptr;
    std::vector<std::shared_ptr<SurfaceMesh<3> > > vm = {mesh_ptr};
    Obj::Read_From_Obj_File<SurfaceMesh<3> >(file, vm);
    mesh = *vm[0];
    std::cout << "elements size: " << mesh.Elements().size() << std::endl;
    std::cout << "vertices size: " << mesh.Vertices().size() << std::endl;
    //bbox
    Vec<3> com = Vec<3>::Zero();
    for(int i=0;i<mesh.Vertices().size();++i){
        com += mesh.Vertices()[i];
    }
    com = com/double(mesh.Vertices().size());
    bbox_min = Vec<3>::Zero();
    bbox_max = Vec<3>::Zero();
    for(int i=0;i<mesh.Vertices().size();++i){
        Vec<3> pv = mesh.Vertices()[i];
        pv = Vec<3>(pv(0)*scale(0),pv(1)*scale(1),pv(2)*scale(2));//any better way of doing this?
        mesh.Vertices()[i] = pv;
        for(int j=0;j<3;++j){
            bbox_min(j) = std::min(bbox_min(j),pv(j));
            bbox_max(j) = std::max(bbox_max(j),pv(j));
        }
    }
}

void RigidExplicitGeometry<2>::Precompute_SDF_I(const Veci<2>& cell_counts, const int narrow_band){
    //determine sdf's grid size
    double dx = 0.0;
    for(int i=0;i<2;++i){
        dx = std::max(dx, 2*(bbox_max(i)-bbox_min(i))/cell_counts(i));
    }
    for(int i=0;i<2;++i){
        grid.domain_min(i) = -(dx*cell_counts(i))/2;
        grid.domain_max(i) = (dx*cell_counts(i))/2;
    }
    grid.Initialize(cell_counts, dx, -Vec<2>(cell_counts(0)*dx/2, cell_counts(1)*dx/2));//let everything match

    phi.counts = grid.cell_counts;
    intp.reset(new Interpolation<2>(grid));
    GH::Mesh_To_SDF<2>(mesh, grid, phi.array, narrow_band);

    double dx2 = grid.dx*grid.dx;
    area = 0;
    inertia = 0;
    iterate_cell_d(iter, grid, 2){
        Veci<2> cell = iter.Coord();
        Vec<2> pos = grid.Center(cell);
        if(phi(cell)<0){
            area += dx2;
            inertia += dx2*(pos.dot(pos));
        }
    }
}

void RigidExplicitGeometry<3>::Precompute_SDF_I(const Veci<3>& cell_counts, const int narrow_band){
    //determine sdf's grid size
    std::cout << "Precompute SDF started" << std::endl;
    std::cout << "bounding box: " << bbox_min.transpose() << " to " << bbox_max.transpose() << std::endl;
    double dx = 0.0;
    for(int i=0;i<3;++i){
        dx = std::max(dx, 2*(bbox_max(i)-bbox_min(i))/cell_counts(i));
    }
    for(int i=0;i<3;++i){
        grid.domain_min(i) = -(dx*cell_counts(i))/2;
        grid.domain_max(i) = (dx*cell_counts(i))/2;
    }
    grid.Initialize(cell_counts, dx, -Vec<3>(cell_counts(0)*dx/2, cell_counts(1)*dx/2, cell_counts(2)*dx/2));//let everything match

    std::cout << "SDF grid domain: " << grid.domain_min.transpose() << " to " << grid.domain_max.transpose() << std::endl;
    std::cout << "SDF cell counts: " << grid.cell_counts.transpose() << std::endl;
    phi.counts = grid.cell_counts;
    intp.reset(new Interpolation<3>(grid));
    GH::Mesh_To_SDF<3>(mesh, grid, phi.array, narrow_band);

    double dx3 = grid.dx*grid.dx*grid.dx;
    volume = 0;
    inertia = Mat<3>::Zero();
    iterate_cell_d(iter, grid, 3){
        Veci<3> cell = iter.Coord();
        Vec<3> pos = grid.Center(cell);
        if(phi(cell)<0){
            volume += dx3;
            inertia(0,0) += dx3*(pos(1)*pos(1)+pos(2)*pos(2));
            inertia(1,1) += dx3*(pos(0)*pos(0)+pos(2)*pos(2));
            inertia(2,2) += dx3*(pos(0)*pos(0)+pos(1)*pos(1));
        }
    }
    std::cout << "volume: " << volume << std::endl;
}

double RigidExplicitGeometry<2>::Get_Phi(const Vec<2>& pos){
    //call Precompute_SDF_I before call this function
    if(grid.Inside(pos))
        return intp->Interpolate_Centers(phi, pos);
    else return pos.norm();
}

double RigidExplicitGeometry<3>::Get_Phi(const Vec<3>& pos){
    //call Precompute_SDF_I before call this function
    if(grid.Inside(pos))
        return intp->Interpolate_Centers(phi, pos);
    else return pos.norm();
}

double RigidExplicitGeometry<2>::Get_Mass(const double density){
    //call Precompute_SDF_I before call this function
    return area*density;
}

double RigidExplicitGeometry<3>::Get_Mass(const double density){
    //call Precompute_SDF_I before call this function
    return volume*density;
}

double RigidExplicitGeometry<2>::Get_Inertia_Body(const double density){
    //call Precompute_SDF_I before call this function
    return inertia*density;
}

Mat<3> RigidExplicitGeometry<3>::Get_Inertia_Body(const double density){
    //call Precompute_SDF_I before call this function
    return inertia*density;
}

bool RigidExplicitGeometry<2>::Is_Inside(const Vec<2>& pos){
    //call Precompute_SDF_I before call this function
    return (Get_Phi(pos) < 0);
}

bool RigidExplicitGeometry<3>::Is_Inside(const Vec<3>& pos){
    //call Precompute_SDF_I before call this function
    return (Get_Phi(pos) < 0);
}

Vec<2> RigidExplicitGeometry<2>::Get_Normal(const Vec<2>& pos){
    Vec<2> normal; for(int i=0;i<2;i++)
        normal[i]=(Get_Phi(pos+Vec<2>::Unit(i)*grid.dx)-Get_Phi(pos-Vec<2>::Unit(i)*grid.dx))/(2*grid.dx);
    return normal.normalized();
}

Vec<3> RigidExplicitGeometry<3>::Get_Normal(const Vec<3>& pos){
    Vec<3> normal; for(int i=0;i<3;i++)
        normal[i]=(Get_Phi(pos+Vec<3>::Unit(i)*grid.dx)-Get_Phi(pos-Vec<3>::Unit(i)*grid.dx))/(2*grid.dx);
    return normal.normalized();
}

void RigidExplicitGeometry<2>::Project_Out(Vec<2>& pos){
    //call Precompute_SDF_I before call this function
    //std::cerr << "project out not implemented" << std::endl;
    Vec<2> normal; for(int i=0;i<2;i++)
        normal[i]=(Get_Phi(pos+Vec<2>::Unit(i)*grid.dx)-Get_Phi(pos-Vec<2>::Unit(i)*grid.dx))/(2*grid.dx);
    normal.normalize();
    pos -= normal*Get_Phi(pos);
}

void RigidExplicitGeometry<3>::Project_Out(Vec<3>& pos){
    //call Precompute_SDF_I before call this function
    //std::cerr << "project out not implemented" << std::endl;
    Vec<3> normal; for(int i=0;i<3;i++)
        normal[i]=(Get_Phi(pos+Vec<3>::Unit(i)*grid.dx)-Get_Phi(pos-Vec<3>::Unit(i)*grid.dx))/(2*grid.dx);
    normal.normalize();
    pos -= normal*Get_Phi(pos);
}

void RigidExplicitGeometry<2>::Project_On(Vec<2>& pos){
    Vec<2> normal; for(int i=0;i<2;i++)
        normal[i]=(Get_Phi(pos+Vec<2>::Unit(i)*grid.dx)-Get_Phi(pos-Vec<2>::Unit(i)*grid.dx))/(2*grid.dx);
    normal.normalize();
    pos -= normal*Get_Phi(pos);
}

void RigidExplicitGeometry<3>::Project_On(Vec<3>& pos){
    Vec<3> normal; for(int i=0;i<3;i++)
        normal[i]=(Get_Phi(pos+Vec<3>::Unit(i)*grid.dx)-Get_Phi(pos-Vec<3>::Unit(i)*grid.dx))/(2*grid.dx);
    normal.normalize();
    pos -= normal*Get_Phi(pos);
}

void RigidExplicitGeometry<2>::Project_Out_Directional(Vec<2>& pos, const Vec<2>& direction){
    std::cerr << "Project_Out_Directional() not implemented yet" << std::endl;
}

void RigidExplicitGeometry<3>::Project_Out_Directional(Vec<3>& pos, const Vec<3>& direction){
    double phi = Get_Phi(pos);
    int iter = 0;
    while(phi<1e-2*grid.dx){
        pos -= phi*direction;
        phi = Get_Phi(pos);
        if(phi>=0) return;
        ++iter; if(iter>1000) {std::cerr << "Project_Out_Directional Error" << std::endl; return;}
    }
}

//thin shell

RigidShellGeometry<3>::RigidShellGeometry(const std::string& file, const Vec<3>& scale, int _axis){
    std::shared_ptr<SurfaceMesh<3> > mesh_ptr;
    std::vector<std::shared_ptr<SurfaceMesh<3> > > vm = {mesh_ptr};
    Obj::Read_From_Obj_File<SurfaceMesh<3> >(file, vm);
    mesh = *vm[0];
    std::cout << "elements size: " << mesh.Elements().size() << std::endl;
    std::cout << "vertices size: " << mesh.Vertices().size() << std::endl;
    //bbox
    bbox_min = Vec<3>::Zero();
    bbox_max = Vec<3>::Zero();
    for(int i=0;i<mesh.Vertices().size();++i){
        Vec<3> pv = mesh.Vertices()[i];
        pv = Vec<3>(pv(0)*scale(0),pv(1)*scale(1),pv(2)*scale(2));//any better way of doing this?
        mesh.Vertices()[i] = pv;
        for(int j=0;j<3;++j){
            bbox_min(j) = std::min(bbox_min(j),pv(j));
            bbox_max(j) = std::max(bbox_max(j),pv(j));
        }
    }
    sgn = (_axis>0)? 1 : -1;
    axis = std::abs(_axis);
}
void RigidShellGeometry<3>::Precompute_SDF_I(const Veci<3>& cell_counts, const int narrow_band){
    //determine sdf's grid size
    std::cout << "Precompute SDF started" << std::endl;
    std::cout << "bounding box: " << bbox_min.transpose() << " to " << bbox_max.transpose() << std::endl;
    double dx = 0.0;
    for(int i=0;i<3;++i){
        dx = std::max(dx, 2*(bbox_max(i)-bbox_min(i))/cell_counts(i));
    }
    for(int i=0;i<3;++i){
        grid.domain_min(i) = -(dx*cell_counts(i))/2;
        grid.domain_max(i) = (dx*cell_counts(i))/2;
    }
    grid.Initialize(cell_counts, dx, -Vec<3>(cell_counts(0)*dx/2, cell_counts(1)*dx/2, cell_counts(2)*dx/2));//let everything match

    std::cout << "SDF grid domain: " << grid.domain_min.transpose() << " to " << grid.domain_max.transpose() << std::endl;
    std::cout << "SDF cell counts: " << grid.cell_counts.transpose() << std::endl;
    phi.counts = grid.cell_counts;
    sphi.counts = grid.cell_counts;
    intp.reset(new Interpolation<3>(grid));
    GH::Mesh_To_SDF_Directional<3>(mesh, grid, phi.array, axis, narrow_band);
    for(int i=0;i<phi.array.size();++i) phi.array[i] = -sgn*phi.array[i];
    GH::Mesh_To_SDF_Thin_Shell<3>(mesh, grid, sphi.array, axis, sgn, narrow_band);

    area = 0;
    inertia = Mat<3>::Zero();
    for(int e=0;e<mesh.Elements().size();++e){
        Veci<3> ele = mesh.Elements()[e];
        ArrayF<Vec<3>,3> p;
        for(int q=0;q<3;++q) p[q] = mesh.Vertices()[ele[q]];
        double tarea = MeshFunc::Triangle_Area(p[0], p[1], p[2]);
        area += tarea;
        Vec<3> pos = (p[0]+p[1]+p[2])/3.0;
        inertia(0,0) += tarea*(pos(1)*pos(1)+pos(2)*pos(2));
        inertia(1,1) += tarea*(pos(0)*pos(0)+pos(2)*pos(2));
        inertia(2,2) += tarea*(pos(0)*pos(0)+pos(1)*pos(1));
    }
    std::cout << "area: " << area << std::endl;
}
double RigidShellGeometry<3>::Get_Phi(const Vec<3>& pos){
    if(grid.Inside(pos))
        return intp->Interpolate_Centers(phi, pos);
    /*else{
        Vec<3> tpos = pos;
        for(int i=0;i<3;++i){
            if(tpos(i) < grid.domain_min(i)) tpos(i) = grid.domain_min(i);
            else if(tpos(i) > grid.domain_max(i)) tpos(i) = grid.domain_max(i);
        }
        if(intp->Interpolate_Centers(phi, tpos)>0) return pos.norm();
        else return -pos.norm();
    }*/
    return pos.norm();
}
double RigidShellGeometry<3>::Get_SPhi(const Vec<3>& pos){
    if(grid.Inside(pos))
        return intp->Interpolate_Centers(sphi, pos);
    return pos.norm();
}
double RigidShellGeometry<3>::Get_Mass(const double density){
    return density*area;
}
Mat<3> RigidShellGeometry<3>::Get_Inertia_Body(const double density){
    return inertia*density;
}
bool RigidShellGeometry<3>::Is_Inside(const Vec<3>& pos){
    return (Get_Phi(pos)<0);
}
Vec<3> RigidShellGeometry<3>::Get_Normal(const Vec<3>& pos){
    Vec<3> normal; for(int i=0;i<3;i++)
        normal[i]=(Get_Phi(pos+Vec<3>::Unit(i)*grid.dx)-Get_Phi(pos-Vec<3>::Unit(i)*grid.dx))/(2*grid.dx);
    return normal.normalized();
}
void RigidShellGeometry<3>::Project_Out(Vec<3>& pos){
    double phi = Get_Phi(pos);
    if(phi>0) return;
    pos -= phi*Get_Normal(pos);
}
void RigidShellGeometry<3>::Project_On(Vec<3>& pos){
    pos -= Get_Normal(pos)*Get_Phi(pos);
}
void RigidShellGeometry<3>::Project_Out_Directional(Vec<3>& pos, const Vec<3>& direction){
    std::cerr << "no directional" << std::endl;
}
