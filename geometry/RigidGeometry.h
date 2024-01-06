#pragma once

#include "Mesh.h"
#include "MarchingCubes.h"
#include "Grid.h"
#include "Type.h"
#include "Particles.h"
#include "Interpolation.h"

#include <Eigen/Sparse>
#include<vector>
#include<string>
#include<iostream>

/**
 * @brief rigid geometry class, both implicit and explicit geometry should actually store a surface mesh
 * @tparam d
 */
template <int d>class RigidGeometry{};

template<> class RigidGeometry<2>
{
public:
    Grid<2> grid;///< implicit geometry uses it for generating mesh, explicit geometry uses it for computing sdf

    SurfaceMesh<2> mesh;
    Vec<2> bbox_min;
    Vec<2> bbox_max;

    RigidGeometry(){};
    void Get_Mesh(const double orientation, const Vec<2>& position, SurfaceMesh<2>& tmesh);
    void Get_Mesh(const Mat<2>& rotation, const Vec<2>& position, SurfaceMesh<2>& tmesh);
    virtual void Precompute_Mesh(double dx){};
    virtual void Precompute_SDF_I(const Veci<2>& cell_counts, const int narrow_band){};
    virtual double Get_Phi(const Vec<2>& pos) = 0;
    virtual double Get_SPhi(const Vec<2>& pos){return 0;};
    virtual double Get_Mass(const double density) = 0;
    virtual double Get_Inertia_Body(const double density) = 0;
    bool Is_Inside(const Vec<2>& pos);
    virtual void Project_Out(Vec<2>& pos) = 0;
    virtual void Project_On(Vec<2>& pos) = 0;
    virtual void Project_Out_Directional(Vec<2>& pos, const Vec<2>& direction) = 0;
    virtual Vec<2> Get_Normal(const Vec<2>& pos) = 0;
};

template <> class RigidGeometry<3>
{
public:
    Grid<3> grid;///< implicit geometry uses it for generating mesh, explicit geometry uses it for computing sdf

    SurfaceMesh<3> mesh;
    Vec<3> bbox_min;
    Vec<3> bbox_max;

    RigidGeometry(){};
    void Get_Mesh(const Mat<3>& orientation, const Vec<3>& position, SurfaceMesh<3>& tmesh);
    virtual void Precompute_Mesh(double dx){};
    virtual void Precompute_SDF_I(const Veci<3>& cell_counts, const int narrow_band){};
    virtual double Get_Phi(const Vec<3>& pos) = 0;
    virtual double Get_SPhi(const Vec<3>& pos){return 0;};
    virtual double Get_Mass(const double density) = 0;
    virtual Mat<3> Get_Inertia_Body(const double density) = 0;
    bool Is_Inside(const Vec<3>& pos);
    virtual void Project_Out(Vec<3>& pos) = 0;
    virtual void Project_On(Vec<3>& pos) = 0;
    virtual void Project_Out_Directional(Vec<3>& pos, const Vec<3>& direction) = 0;
    virtual Vec<3> Get_Normal(const Vec<3>& pos) = 0;
};


/**
 * @brief explicit geometry class, note that sdf should be precomputed while system initialization
 * @tparam d
 */
template <int d> class RigidExplicitGeometry :public RigidGeometry<d>{};

template <> class RigidExplicitGeometry<2> : public RigidGeometry<2>{
public:
    Field<double, 2> phi;
    std::unique_ptr<Interpolation<2> > intp=nullptr;
    double area;
    double inertia;

    RigidExplicitGeometry(){};
    RigidExplicitGeometry(const std::string& file, const Vec<2>& scale);
    virtual void Precompute_SDF_I(const Veci<2>& cell_counts, const int narrow_band);
    double Get_Phi(const Vec<2>& pos);
    double Get_Mass(const double density);
    double Get_Inertia_Body(const double density);
    bool Is_Inside(const Vec<2>& pos);
    void Project_Out(Vec<2>& pos);
    void Project_On(Vec<2>& pos);
    void Project_Out_Directional(Vec<2>& pos, const Vec<2>& direction);
    virtual Vec<2> Get_Normal(const Vec<2>& pos);
};

template <> class RigidExplicitGeometry<3> : public RigidGeometry<3>{
public:
    Field<double, 3> phi;
    std::unique_ptr<Interpolation<3> > intp=nullptr;
    double volume;
    Mat<3> inertia;

    RigidExplicitGeometry(){};
    RigidExplicitGeometry(const std::string& file, const Vec<3>& scale);
    virtual void Precompute_SDF_I(const Veci<3>& cell_counts, const int narrow_band);
    double Get_Phi(const Vec<3>& pos);
    double Get_Mass(const double density);
    Mat<3> Get_Inertia_Body(const double density);
    bool Is_Inside(const Vec<3>& pos);
    void Project_Out(Vec<3>& pos);
    void Project_On(Vec<3>& pos);
    void Project_Out_Directional(Vec<3>& pos, const Vec<3>& direction);
    virtual Vec<3> Get_Normal(const Vec<3>& pos);
};

template <int d> class RigidShellGeometry :public RigidGeometry<d>{};

template <> class RigidShellGeometry<2> : public RigidGeometry<2>{
public:
    Field<double, 2> phi;
    std::unique_ptr<Interpolation<2> > intp=nullptr;
    double inertia;
    int axis;

    RigidShellGeometry(){};
    RigidShellGeometry(const std::string& file, const Vec<2>& scale, int _axis){};
    virtual void Precompute_SDF_I(const Veci<3>& cell_counts, const int narrow_band){};
    double Get_Phi(const Vec<2>& pos){return 0;};
    double Get_Mass(const double density){return 0;};
    double Get_Inertia_Body(const double density){return 0;};
    bool Is_Inside(const Vec<2>& pos){return false;};
    void Project_Out(Vec<2>& pos){};
    void Project_On(Vec<2>& pos){};
    void Project_Out_Directional(Vec<2>& pos, const Vec<2>& direction){};
    virtual Vec<2> Get_Normal(const Vec<2>& pos){return Vec<2>::Zero();};
};
template <> class RigidShellGeometry<3> : public RigidGeometry<3>{
public:
    Field<double, 3> phi;
    Field<double, 3> sphi;
    std::unique_ptr<Interpolation<3> > intp=nullptr;
    double area;
    Mat<3> inertia;
    int axis;
    int sgn;

    RigidShellGeometry(){};
    RigidShellGeometry(const std::string& file, const Vec<3>& scale, int _axis);
    virtual void Precompute_SDF_I(const Veci<3>& cell_counts, const int narrow_band);
    double Get_Phi(const Vec<3>& pos);
    double Get_SPhi(const Vec<3>& pos);
    double Get_Mass(const double density);
    Mat<3> Get_Inertia_Body(const double density);
    bool Is_Inside(const Vec<3>& pos);
    void Project_Out(Vec<3>& pos);
    void Project_On(Vec<3>& pos);
    void Project_Out_Directional(Vec<3>& pos, const Vec<3>& direction);
    virtual Vec<3> Get_Normal(const Vec<3>& pos);
};
