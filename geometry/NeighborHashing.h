//
// Created by zhangjie on 2024/1/6.
//

#ifndef SIMPLEXOLDCLION_NEIGHBORHASHING_H
#define SIMPLEXOLDCLION_NEIGHBORHASHING_H
#include "NeighborSearcher.h"
#include "SpatialHashing.h"


template<int d, int bs>
class NeighborHashing : public NeighborSearcher<d> {
    Typedef_VectorDii(d);
    using ArrayNBi = ArrayF<int, bs * 4>;
public:
    std::shared_ptr<Array<VectorD> >  points_ptr;
    std::shared_ptr<SpatialHashing<d, bs> > spatial_hashing;
public:
    NeighborHashing(real dx) {
        VectorDi cell_counts = VectorDi::Ones() * 32;	////this number can set to be arbitrary, as the spatial hashing table is infinitely large
        Grid<d> grid;
        grid.Initialize(cell_counts, dx);
        spatial_hashing = std::make_shared<SpatialHashing<d, bs> >(grid);
    }
    virtual void Build_Data(Array<VectorD>& arr) {
        points_ptr = std::make_shared<Array<VectorD> >(arr);
        spatial_hashing->Update_Voxels(arr);
    }
    virtual size_t Find_Neighbors(const VectorD& pos, const real& radius, Array<int>& results, bool append = false)const {
        if (!append) results.clear();
        ArrayNBi temp_res;
        if (spatial_hashing->Find_Nbs(pos, *points_ptr, radius, temp_res)) {
            size_t num = temp_res[0];
            for (size_t i = 0; i < num; i++) {
                results.push_back(temp_res[i + 1]);
            }
            return num;
        }
        else {
            std::cerr << "Error: NeighborHashing::Find_Neighbors fail for pos " << pos.transpose() << std::endl;
            return 0;
        }
    }
    virtual int Find_Nearest_Nb(const VectorD& pos)const {
        return spatial_hashing->Find_Nearest_Nb(pos, *points_ptr);
    }
    virtual int Find_K_Nearest_Nb(const VectorD& pos, int k, Array<int>& results)const {
        std::cout << "NeighborHashing::Find_K_Nearest_Nb error: function not implemented\n";
        exit(0);
        return -1;
    }
};

#endif

