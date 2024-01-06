//
// Created by zhangjie on 2024/1/6.
//

#ifndef SIMPLEXOLDCLION_NEIGHBORKDTREE_H
#define SIMPLEXOLDCLION_NEIGHBORKDTREE_H
#include "NeighborSearcher.h"
#include "nanoflann/nanoflann.hpp"

template<int d>
class PointSetAdapter {
    //Declare_Eigen_Types(double, d);
    Typedef_VectorDii(d);
public:
    std::shared_ptr<Array<VectorD> > points_ptr;
    void Initialize(Array<VectorD>& arr);
    // Interface required by nanoflann
    size_t kdtree_get_point_count() const;
    real kdtree_get_pt(const size_t idx, const size_t dim) const;
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
};

template<int d>
class NeighborKDTree: public NeighborSearcher<d> {
    //Declare_Eigen_Types(double, d);
    using Base = NeighborSearcher<d>;
    Typedef_VectorDii(d);
public:
    const int max_leaf = 10;/* max leaf */
    PointSetAdapter<d> points;
    using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<real, PointSetAdapter<d> >,//NOTE: It's actually squared L2 norm
    PointSetAdapter<d>,
    d /* dim */
    >;
    my_kd_tree_t index;
public:
    NeighborKDTree() :index(d, points, nanoflann::KDTreeSingleIndexAdaptorParams(max_leaf)) { index.buildIndex(); }
    virtual void Build_Data(Array<VectorD>& arr);
    virtual size_t Find_Neighbors(const VectorD& pos, const real& radius, Array<int>& results, bool append = false)const;
    virtual int Find_Nearest_Nb(const VectorD& pos)const;
    virtual int Find_K_Nearest_Nb(const VectorD& pos, int k, Array<int>& results)const;
};

#endif


