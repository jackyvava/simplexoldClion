//
// Created by zhangjie on 2024/1/6.
//

#include "NeighborSearcher.h"

template<int d>
void NeighborSearcher<d>::Update_Points(Array<VectorD>& points)
{
    search_results.clear();
    this->Build_Data(points);
}

template<int d>
void NeighborSearcher<d>::Update_Points(Array<VectorD>& points, FilterFunc& filter_func)
{
    Array<VectorD> temp_array;
    temp_array.clear(); temp_array.reserve(points.size());
    for (size_t i = 0; i < points.size(); i++) {
        if (filter_func((int)i)) temp_array.push_back(points[i]);
    }
    this->Update_Points(temp_array);
}

template<int d>
size_t NeighborSearcher<d>::Find_Neighbors(const VectorD& pos, const real& radius, FilterFunc& filter_func, Array<int>& results, bool append) const
{
    if (!append) results.clear();
    Array<int> temp_results;
    this->Find_Neighbors(pos, radius, temp_results, false);//append=false, temp_results is cleared
    size_t num = 0;
    for (size_t i = 0; i < temp_results.size(); i++) {
        if (filter_func(temp_results[i])) {
            num++;
            results.push_back(temp_results[i]);
        }
    }
    return num;
}

template<int d>
Array<int> NeighborSearcher<d>::Find_Neighbors(const VectorD& pos, const real& radius) const
{
    Array<int> temp_results;
    this->Find_Neighbors(pos, radius, temp_results, false);//append=false, temp_results is cleared
    return temp_results;
}

template<int d>
Array<int> NeighborSearcher<d>::Find_Neighbors(const VectorD& pos, const real& radius, FilterFunc& filter_func) const
{
    Array<int> temp_results;
    this->Find_Neighbors(pos, radius, filter_func, temp_results, false);//append=false, temp_results are cleared
    return temp_results;
}

template<int d>
ArraySlice<int> NeighborSearcher<d>::Record_Neighbors(const VectorD& pos, const real& radius)
{
    size_t beg = search_results.size();
    this->Find_Neighbors(pos, radius, search_results, true);//append=true
    size_t end = search_results.size();
    return ArraySlice<int>(beg, end, search_results);
}

template<int d>
ArraySlice<int> NeighborSearcher<d>::Record_Neighbors(const VectorD& pos, const real& radius, FilterFunc& filter_func)
{
    size_t beg = search_results.size();
    this->Find_Neighbors(pos, radius, filter_func, search_results, true);//append=true
    size_t end = search_results.size();
    return ArraySlice<int>(beg, end, search_results);
}




template class NeighborSearcher<2>;
template class NeighborSearcher<3>;