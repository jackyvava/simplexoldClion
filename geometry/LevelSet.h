//
// Created by zhangjie on 2024/1/6.
//

#ifndef SIMPLEXOLDCLION_LEVELSET_H
#define SIMPLEXOLDCLION_LEVELSET_H


#include "Grid.h"
#include "Field.h"
#include "FaceField.h"
#include "Interpolation.h"
#include "GeometryPrimitives.h"

template<int d> class LevelSet
{Typedef_VectorDii(d);
public:
    Grid<d> grid;
    Field<real,d> phi;
protected:
    std::unique_ptr<Interpolation<d> > intp=nullptr;
public:
    LevelSet(){}
    LevelSet(const Grid<d>& _grid);
    void Initialize(const Grid<d>& _grid);

    virtual void Set_By_Geom(ImplicitGeometry<d>& geom);
    virtual real Phi(const VectorD& pos) const;
    virtual VectorD Normal(const VectorD& pos) const;	////TOFIX: fix finite difference on the boundary
    virtual VectorD Gradient(const VectorD& pos) const;	////TOFIX: fix finite difference on the boundary
    virtual VectorD Closest_Point(const VectorD& pos,const real epsilon=(real)0) const;
    virtual VectorD Closest_Point_With_Iterations(const VectorD& pos,const int max_iter=5) const;
    virtual void Update_Normals(FaceField<real,d>& normals) const;
    virtual void Update_Normals(Field<VectorD,d>& normals) const;
    virtual real Curvature(const VectorD& pos) const;

    ////Helper functions
    static real Sign(const real phi){return phi<=(real)0?(real)-1:(real)1;}
    static bool Interface(const real phi_1,const real phi_2){return Sign(phi_1)!=Sign(phi_2);}
    static real Theta(const real phi_1,const real phi_2){return phi_1/(phi_1-phi_2);}

    //////////////////////////////////////////////////////////////////////////
    ////Fast marching
    void Fast_Marching(const real band_width=(real)-1);
protected:
    ////Fast marching auxiliary data structures
    Field<real,d> tent;
    Array<char> done;
    Array<int> intf_cells;
    using PRI=std::pair<real,int>;
    std::priority_queue<PRI,Array<PRI>,std::greater<PRI> > heap;
    ////Fast marching helper functions
    void Initialize_Fast_Marching();
    void Precondition(const real band_width);
    real Solve_Eikonal(const VectorDi &cell);
    bool Solve_Quadratic(const real p1,const real p2,const real dx,real& rst);
    bool Solve_Quadratic(const real p1,const real p2,const real p3,const real dx,real& rst);
};

#endif
