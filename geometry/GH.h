#pragma once

#include "Type.h"
#include "Mesh.h"
#include "Grid.h"
#include "LevelSet.h"
#include "MarchingCubes.h"
#include "Codimension.h"

#include <algorithm>
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Geometry>

namespace GH{
    /**
     * @brief the distance from a point to a line segment
     * @tparam d 
     * @param x0 point
     * @param x1 segment vertice 1
     * @param x2 sgement vertice 2
     * @return double distance
     */
    template<int d>
    double Point_Segment_Distance(const Vec<d>& x0, const Vec<d>& x1, const Vec<d>& x2);

    /**
     * @brief the distance from a point to a triangle, only in 3D
     * @param x0 point
     * @param x1 triangle vertice 1
     * @param x2 triangle vertice 2
     * @param x3 triangle vertice 3
     * @return double 
     */
    double Point_Triangle_Distance(const Vec<3>& x0, const Vec<3>& x1, const Vec<3>& x2, const Vec<3>& x3);

    /**
     * @brief robust test of (x0,y0) in the triangle (x1,y1)-(x2,y2)-(x3,y3)
     *  if true is returned, the barycentric coordinates are set in a,b,c.
     * @param x0 point
     * @param x1 triangle vertice 1
     * @param x2 triangle vertice 2
     * @param x3 triangle vertice 3
     * @param a barycentric coordinate to x1
     * @param b barycentric coordinate to x2
     * @param c barycentric coordinate to x3
     * @return true iff in triangle
     * @return false iff out of triangle
     */
    bool Point_In_Triangle_2D(const Vec<2>& x0, const Vec<2>& x1, const Vec<2>& x2, const Vec<2>& x3, double& a, double& b, double& c);

    /**
     * @brief robust overlap test of segment (x1, x2) and box
     * 
     * @param box_center 
     * @param box_half_size 
     * @param x1 segment vertice 1
     * @param x2 segment vertice 2
     * @return true overlaped
     * @return false not overlaped
     */
    bool Segment_Box_Overlap_2D(const Vec<2>& box_center, const Vec<2>& box_half_size, const Vec<2>& x1, const Vec<2>& x2);

    /**
     * @brief robust overlap test of plane and box, modified from: http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/tribox3.txt
     * @param normal the normal vector of the plane
     * @param vert one vertex on the plane
     * @param box_half_size half of the size of the box
     * @return true overlaped
     * @return false not overlaped
     */
    bool Plane_Box_Overlap_3D(const Vec<3>& normal, const Vec<3>& vert, const Vec<3>& box_half_size);

    /**
     * @brief robust overlap test of triangle (x1, x2, x3) and box, modified from: http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/tribox3.txt
     * @param box_center center position of box
     * @param box_half_size half of box's size in 3 axis
     * @param x1 triangle vertice 1
     * @param x2 triangle vertice 2
     * @param x3 triangle vertice 3
     * @return true overlaped
     * @return false not overlaped
     */
    bool Triangle_Box_Overlap_3D(const Vec<3>& box_center, const Vec<3>& box_half_size, const Vec<3>& x1, const Vec<3>& x2, const Vec<3>& x3);

    /**
     * @brief convert mesh to sdf, modified from: https://github.com/christopherbatty/SDFGen
     * @tparam d 
     * @param mesh 
     * @param grid 
     * @param phi the signed distance field
     * @param narrow_band sdf in this narrow band will be computed precisely, other place using sweep
     */
    template<int d>
    void Mesh_To_SDF(const SurfaceMesh<d>& mesh, const Grid<d>& grid, std::vector<double>& phi, int narrow_band=2);

    template<int d>
    void Mesh_To_SDF(const Codimension::CodimMesh<d>& mesh, const Codimension::CodimMesher<d>& mesher, const Grid<d>& grid, std::vector<double>& phi, int narrow_band=2);

    template<int d>
    void Mesh_To_SDF_Narrow_Band(const Codimension::CodimMesh<d>& mesh, const Codimension::CodimMesher<d>& mesher, const Grid<d>& grid, std::vector<double>& phi, int narrow_band=2);

    /**
     * @brief similar to Mesh_To_SDF, designed for fluid surface mesh
     * @tparam d 
     * @param mesh 
     * @param grid 
     * @param phi the signed distance field
     * @param axis The axis for gravity, used for ray's direction when determine sign, 0 for x axis, 1 for y axis, 2 for z axis
     * @param narrow_band sdf in this narrow band will be computed precisely, other place using sweep
     */
    template<int d>
    void Mesh_To_SDF_Directional(const SurfaceMesh<d>& mesh, const Grid<d>& grid, std::vector<double>& phi, int axis, int narrow_band=2);

    /**
     * @brief similar to Mesh_To_SDF, designed for fluid codimension mesh
     * @tparam d 
     * @param mesh 
     * @param grid 
     * @param phi the signed distance field
     * @param axis The axis for gravity, used for ray's direction when determine sign, 0 for x axis, 1 for y axis, 2 for z axis
     * @param narrow_band sdf in this narrow band will be computed precisely, other place using sweep
     */
    template<int d>
    void Mesh_To_SDF_Directional(const Codimension::CodimMesh<d>& mesh, const Codimension::CodimMesher<d>& mesher, const Grid<d>& grid, std::vector<double>& phi, int axis, int narrow_band=2);

    template<int d>
    void Mesh_To_SDF_Directional_Narrow_Band(const Codimension::CodimMesh<d>& mesh, const Codimension::CodimMesher<d>& mesher, const Grid<d>& grid, std::vector<double>& phi, int axis, int narrow_band=2);

    template<int d>
    void Mesh_To_SDF_Thin_Shell(const SurfaceMesh<d>& mesh, const Grid<d>& grid, std::vector<double>& phi, int axis, int sign, int narrow_band=2);
};