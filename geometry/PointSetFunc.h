//
// Created by zhangjie on 2024/1/6.
//

#ifndef SIMPLEXOLDCLION_POINTSETFUNC_H
#define SIMPLEXOLDCLION_POINTSETFUNC_H

#include "Mesh.h"
#include "MeshFunc.h"
#include "MeshAdvFunc.h"
#include "AuxFunc.h"
#include "Constants.h"
#include "GeometryParticles.h"

namespace PointSetFunc
{
    using namespace AuxFunc;

    //////////////////////////////////////////////////////////////////////////
    ////2D points
    real Initialize_Circle_Points(const Vector2& c, const real R, const int p_num, GeometryParticles<2>& particles);
    real Initialize_Oval_Points(const Vector2& c, const real R, const real a, const real b, const int p_num, GeometryParticles<2>& particles);
    int Initialize_Circle_Points(const Vector2& c, const real R, const real dx, GeometryParticles<2>& particles);
    void Initialize_Rectangle_Points(const int nx, const int ny, const real dx, const Vector2& start, GeometryParticles<2>& particles);
    void Initialize_Round_Corner_Rectangle_Points(const int nx, const int ny, const real dx, const real r, const Vector2& start, GeometryParticles<2>& particles);
    void Initialize_Segment_Points(const Vector2& v1, const Vector2& v2, int N, GeometryParticles<2>& particles);
    void Initialize_Curve_Points(const Vector2& c,const real R,const real theta,const int N,GeometryParticles<2>& particles);

    //////////////////////////////////////////////////////////////////////////
    ////3D points
    real Initialize_Sphere_Points(const Vector3& c, const real R, const int sub, GeometryParticles<3>& particles);
    real Add_Sphere_Points(const Vector3& c, const real R, const int sub, GeometryParticles<3>& particles);
    real Initialize_Circle_Points(const Vector3& c, const real R, const Vector3& normal, const real dx, GeometryParticles<3>& particles);
    std::vector<int> Initialize_Circle_Points_Grid(const Vector3& c, const real R, const real dx, GeometryParticles<3>& particles);
    std::vector<int> Initialize_Circle_Points_Grid_Random(const Vector3& c, const real R, const real dx, GeometryParticles<3>& particles);
    std::vector<int> Initialize_Circle_Points_Random(const Vector3& c, const real R, const real dx, GeometryParticles<3>& particles);
    std::vector<int> Initialize_Catenoid_Points(const Vector3& c, const real R, const real dx, GeometryParticles<3>& particles);
    std::vector<int> Initialize_Catenoid_Points_Random(const Vector3& c, const real R, const real dx, GeometryParticles<3>& particles);
    void Initialize_Rectangle_Points(const int nx, const int ny, const real dx, const Vector3& start, const Vector3& normal, GeometryParticles<3>& particles);
    void Initialize_Lattice_Points(const Vector3& domain_min, const Vector2i& counts, Vector3 k1, Vector3 k2, const real dx, GeometryParticles<3>& particles);//nx*ny cells, not nodes
    std::vector<int> Initialize_Lattice_Points2(const Vector3& domain_min, const Vector2i& counts, Vector3 k1, Vector3 k2, const real dx, GeometryParticles<3>& particles);
    real Initialize_Ellipsoid_Points(const Vector3& C, const real R, const int sub, const real a, const real b, const real c, GeometryParticles<3>& particles);
    std::vector<int> Initialize_Catenoid_Points(const Vector3& center, const real R, const int nr, const real height, const int nh, GeometryParticles<3>& particles);

    //////////////////////////////////////////////////////////////////////////
    ////File IO
    template<int d> void Write_Local_Frames_To_File(const std::string file_name, const GeometryParticles<d>& points, const real scale = (real).02)
    {
        SegmentMesh<3> s3;
        for (int i = 0; i < points.Size(); i++) {
            if (points.I(i) == -1)continue;
            for (int j = 0; j < d; j++) {
                Vector<real, d> e = points.E(i).col(j);
                s3.Vertices().push_back(V<3>(points.X(i)));
                Vector<real, d> s_end = points.X(i) + e * scale;
                s3.Vertices().push_back(V<3>(s_end));
                int s = (int)s3.Vertices().size();
                s3.Elements().push_back(Vector2i(s - 2, s - 1));
            }
        }
        s3.Write_To_File_3d(file_name);
    }

    template<int d> void Write_Tracker_Circles_To_File(const std::string file_name, const GeometryParticles<d>& points)
    {
        int pn = points.Size(); Array<Vector<real, d> > normals(pn);
        for (int i = 0; i < pn; i++) {
            normals[i] = points.Normal(i);
            //normals[i] = Fit_Vector<real, d>(0, 0, 1).normalized();
        }
        Write_Vectors_To_File_3d_Fast<d, real>(points.XRef(), normals, file_name);
    }
};
#endif
