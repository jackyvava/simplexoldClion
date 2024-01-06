//
// Created by zhangjie on 2024/1/6.
//

#include "Particles.h"
#include "File.h"
#include "AuxFunc.h"

//////////////////////////////////////////////////////////////////////////
////points

template<int d,typename T> void Points<d,T>::Write_Binary(std::ostream &output) const
{
    int n=Size();File::Write_Binary(output,n);
    if(n>0){File::Write_Binary_Array(output,&(*X())[0],n);}
}

template<int d,typename T> void Points<d,T>::Read_Binary(std::istream &input)
{
    int n=0;File::Read_Binary(input,n);
    Resize(n);
    if(n>0){File::Read_Binary_Array(input,&(*X())[0],n);}
}

template<int d,typename T> void Points<d,T>::Write_To_File_3d(const std::string& file_name) const
{
    if constexpr (d==3){
        File::Write_Binary_To_File(file_name,*this);}
    else{
        Points<3,T> p3;p3.Resize(Size());
        AuxFunc::Dim_Conversion_Array<T,d,3>(*X(),*p3.X(),(T)0);
        File::Write_Binary_To_File(file_name,p3);}
}

template<int d, typename T> void Points<d, T>::Print_Attributes(void) {
    std::cout << "Total " << att_map.size() << " attributes:\n";
    for (auto iter = att_map.begin(); iter != att_map.end(); iter++) {
        std::cout << (iter->first) << " ";
    }
    std::cout << "\n";
}

template<int d, typename T> void Points<d, T>::Resize(const int size) {
    for (auto iter = att_map.begin(); iter != att_map.end(); iter++) {
        (iter->second).Resize(size);
    }
}

template<int d, typename T> void Points<d, T>::Reserve(const int size) {
    for (auto iter = att_map.begin(); iter != att_map.end(); iter++) {
        (iter->second).Reserve(size);
    }
}

template<int d, typename T> int Points<d, T>::Add_Element() {
    int idx = -1;
    for (auto iter = att_map.begin(); iter != att_map.end(); iter++) {
        int now_idx = (iter->second).Add_Element();
        if (idx != -1 && idx != now_idx) {
            std::cerr << "Points::Add_Element error: unexpected size of " << iter->first << "\n";
            exit(0);
        }
        idx = now_idx;
    }
    return idx;
}

template<int d, typename T> int Points<d, T>::Add_Elements(const int n) {
    int idx = -1;
    for (auto iter = att_map.begin(); iter != att_map.end(); iter++) {
        int now_idx = (iter->second).Add_Elements(n);
        if (idx != -1 && idx != now_idx) {
            std::cerr << "Points::Add_Element error: unexpected size of " << iter->first << "\n";
            exit(0);
        }
        idx = now_idx;
    }
    return idx;
}

template<int d, typename T> void Points<d, T>::Copy_Element_From(const int idx, const Points<d, T>& src, const int src_idx) {
    for (auto iter = att_map.begin(); iter != att_map.end(); iter++) {
        (iter->second).Copy_Element(idx, src.att_map.at(iter->first), src_idx);
    }
}

template<int d, typename T> int Points<d, T>::Delete_Elements(const Array<bool>& is_deleted) {
    int deleted_size = 0;
    for (auto iter = att_map.begin(); iter != att_map.end(); iter++) {
        deleted_size = (iter->second).Delete_Elements(is_deleted);
    }
    return deleted_size;
}

template<int d, typename T>
bool Points<d, T>::Save_Snapshot(const std::string& save_folder, const std::string& suffix)
{
    for (auto iter = att_map.begin(); iter != att_map.end(); iter++) {
        std::string file_name = save_folder + "/" + iter->first + suffix;
        if (!(iter->second).Write_Binary(file_name)) {
            std::cerr << "Points::Save_Snapshot error: cannot write file " << file_name << "\n";
            exit(0);
            return false;
        }
    }
    return true;
}

template<int d, typename T>
bool Points<d, T>::Load_Snapshot(const std::string& save_folder, const std::string& suffix)
{
    for (auto iter = att_map.begin(); iter != att_map.end(); iter++) {
        std::string file_name = save_folder + "/" + iter->first + suffix;
        if (!(iter->second).Read_Binary(file_name)) {
            std::cerr << "Points::Load_Snapshot error: cannot read file " << file_name << "\n";
            exit(0);
            return false;
        }
    }
    return true;
}



template<int d,typename T> void Points<d,T>::Copy_Attribute_Self(const int to,const int from)
{X(to)=X(from);}

template<int d,typename T> void Points<d,T>::Copy_Attribute_From(const int to,Points<d,T>& src,const int from)
{X(to)=src.X(from);}

template<int d, typename T> int Points<d, T>::Add_Point(VectorD pos) {
    X()->push_back(pos);
    int n = Size();
    this->Resize(n);//Critical: use this to access child class
    return n - 1;
}

template<int d, typename T> int Points<d, T>::Add_Points(const Array<VectorD>& vertices) {
    int n=(int)vertices.size();
    int idx = this->Add_Elements(n);
    for (int i = 0; i < n; i++) {
        this->X(i + idx) = vertices[i];}
    return idx;
}

template class Points<2,double>;
template class Points<3,double>;
template class Points<2,float>;
template class Points<3,float>;

//////////////////////////////////////////////////////////////////////////
////tracker points

template<int d,typename T> void TrackerPoints<d,T>::Write_Binary(std::ostream& output) const
{
    int n=Size();File::Write_Binary(output,n);
    if(n>0){
        File::Write_Binary_Array(output,&(*X())[0],n);
        File::Write_Binary_Array(output,&(*V())[0],n);
        File::Write_Binary_Array(output,&(*I())[0],n);}
}

template<int d,typename T> void TrackerPoints<d,T>::Read_Binary(std::istream& input)
{
    int n=0;File::Read_Binary(input,n);
    Resize(n);
    if(n>0){
        File::Read_Binary_Array(input,&(*X())[0],n);
        File::Read_Binary_Array(input,&(*V())[0],n);
        File::Read_Binary_Array(input,&(*I())[0],n);}
}

template<int d,typename T> void TrackerPoints<d,T>::Write_To_File_3d(const std::string& file_name) const
{
    if constexpr (d==3){
        File::Write_Binary_To_File(file_name,*this);}
    else{
        TrackerPoints<3,T> p3;p3.Resize(Size());
        AuxFunc::Dim_Conversion_Array<T,d,3>(*X(),*p3.X(),(T)0);
        AuxFunc::Dim_Conversion_Array<T,d,3>(*V(),*p3.V(),(T)0);
        *p3.I()=*I();
        File::Write_Binary_To_File(file_name,p3);}
}

template<int d,typename T>
void Points<d,T>::Write_To_File_3d_Fast(const std::string& file_name) const
{
    ////write x only
    int n=Size();

    float* xf=new float[n*4];
    memset((void*)xf,0,sizeof(float)*n*4);
#pragma omp parallel for
    for(int i=0;i<n;i++){
        if constexpr (d==2){
            xf[i*4]=(float)(*x)[i][0];
            xf[i*4+1]=(float)(*x)[i][1];}
        else if constexpr (d==3){
            xf[i*4]=(float)(*x)[i][0];
            xf[i*4+1]=(float)(*x)[i][1];
            xf[i*4+2]=(float)(*x)[i][2];}}
    std::ofstream output(file_name,std::ios::binary);if(!output)return;
    File::Write_Binary(output,n*4);
    File::Write_Binary_Array(output,xf,n*4);
    delete [] xf;
}

template class TrackerPoints<2,double>;
template class TrackerPoints<3,double>;
template class TrackerPoints<2,float>;
template class TrackerPoints<3,float>;

//////////////////////////////////////////////////////////////////////////
////particles

template<int d,typename T> void Particles<d,T>::Write_Binary(std::ostream &output) const
{
    int n=Size();File::Write_Binary(output,n);
    if(n>0){
        File::Write_Binary_Array(output,&(*X())[0],n);
        File::Write_Binary_Array(output,&(*V())[0],n);
        File::Write_Binary_Array(output,&(*F())[0],n);
        File::Write_Binary_Array(output,&(*M())[0],n);
        File::Write_Binary_Array(output,&(*C())[0],n);
        File::Write_Binary_Array(output,&(*I())[0],n);}
}

template<int d,typename T> void Particles<d,T>::Read_Binary(std::istream &input)
{
    int n=0;File::Read_Binary(input,n);
    Resize(n);
    if(n>0){
        File::Read_Binary_Array(input,&(*X())[0],n);
        File::Read_Binary_Array(input,&(*V())[0],n);
        File::Read_Binary_Array(input,&(*F())[0],n);
        File::Read_Binary_Array(input,&(*M())[0],n);
        File::Read_Binary_Array(input,&(*C())[0],n);
        File::Read_Binary_Array(input,&(*I())[0],n);}
}

template<int d,typename T> void Particles<d,T>::Write_To_File_3d(const std::string& file_name) const
{
    if constexpr (d==3){
        File::Write_Binary_To_File(file_name,*this);}
    else{
        Particles<3,T> p3;p3.Resize(Size());
        AuxFunc::Dim_Conversion_Array<T,d,3>(*X(),*p3.X(),(T)0);
        AuxFunc::Dim_Conversion_Array<T,d,3>(*V(),*p3.V(),(T)0);
        AuxFunc::Dim_Conversion_Array<T,d,3>(*F(),*p3.F(),(T)0);
        *p3.M()=*M();
        *p3.C()=*C();
        *p3.I()=*I();
        File::Write_Binary_To_File(file_name,p3);}
}

template<int d,typename T> void Particles<d,T>::Copy_Attribute_Self(const int to,const int from)
{X(to)=X(from);V(to)=V(from);F(to)=F(from);M(to)=M(from);C(to)=C(from);I(to)=I(from);}

template<int d,typename T> void Particles<d,T>::Copy_Attribute_From(const int to,Particles<d,T>& src,const int from)
{X(to)=src.X(from);V(to)=src.V(from);F(to)=src.F(from);M(to)=src.M(from);C(to)=src.C(from);I(to)=src.I(from);}

template<int d,typename T> void Particles<d,T>::Remove(const std::set<int>& indices_to_remove)
{
    int n=Size()-(int)indices_to_remove.size();
    int i=Size()-1;for(int j=Size()-1;j>=0;j--){
        if(indices_to_remove.find(j)==indices_to_remove.end())continue;
        if(i!=j)Copy_Attribute_Self(j,i);i--;}
    Resize(n);
}

template class Particles<2,double>;
template class Particles<3,double>;
template class Particles<2,float>;
template class Particles<3,float>;

template<int d,class T> void Write_To_File_3d_Fast(const Array<Vector<T,d> >& X,const std::string& file_name)
{
    int n=(int)X.size();
    float* xf=new float[n*4];
    memset((void*)xf,0,sizeof(float)*n*4);
#pragma omp parallel for
    for(int i=0;i<n;i++){
        if constexpr (d==2){
            xf[i*4]=(float)X[i][0];
            xf[i*4+1]=(float)X[i][1];}
        else if constexpr (d==3){
            xf[i*4]=(float)X[i][0];
            xf[i*4+1]=(float)X[i][1];
            xf[i*4+2]=(float)X[i][2];}}
    std::ofstream output(file_name,std::ios::binary);if(!output)return;
    File::Write_Binary(output,n*4);
    File::Write_Binary_Array(output,xf,n*4);
    delete [] xf;
}

template void Write_To_File_3d_Fast<2,real>(const Array<Vector2>&,const std::string&);
template void Write_To_File_3d_Fast<3,real>(const Array<Vector3>&,const std::string&);
template void Write_To_File_3d_Fast<2,float>(const Array<Vector2f>&,const std::string&);
template void Write_To_File_3d_Fast<3,float>(const Array<Vector3f>&,const std::string&);

template<int d,class T> void Write_Segments_To_File_3d_Fast(const Array<Vector<T,d> >& X,const Array<Vector<T,d> >& V,const std::string& file_name,
                                                            const bool use_v_as_displacement/*=true*/)
{
    int n=(int)X.size();
    float* xf=new float[n*8];
    memset((void*)xf,0,sizeof(float)*n*8);
    if(use_v_as_displacement){
        ////use v as the first endpoint
#pragma omp parallel for
        for(int i=0;i<n;i++){
            if constexpr (d==2){
                xf[i*8]=(float)X[i][0];
                xf[i*8+1]=(float)X[i][1];

                xf[i*8+4]=(float)(X[i][0]+V[i][0]);
                xf[i*8+5]=(float)(X[i][1]+V[i][1]);}

            else if constexpr (d==3){
                xf[i*8]=(float)X[i][0];
                xf[i*8+1]=(float)X[i][1];
                xf[i*8+2]=(float)X[i][2];

                xf[i*8+4]=(float)(X[i][0]+V[i][0]);
                xf[i*8+5]=(float)(X[i][1]+V[i][1]);
                xf[i*8+6]=(float)(X[i][2]+V[i][2]);}}}
    else{	////use v as the second endpoint
#pragma omp parallel for
        for(int i=0;i<n;i++){
            if constexpr (d==2){
                xf[i*8]=(float)X[i][0];
                xf[i*8+1]=(float)X[i][1];

                xf[i*8+4]=(float)(V[i][0]);
                xf[i*8+5]=(float)(V[i][1]);}

            else if constexpr (d==3){
                xf[i*8]=(float)X[i][0];
                xf[i*8+1]=(float)X[i][1];
                xf[i*8+2]=(float)X[i][2];

                xf[i*8+4]=(float)(V[i][0]);
                xf[i*8+5]=(float)(V[i][1]);
                xf[i*8+6]=(float)(V[i][2]);}}}

    std::ofstream output(file_name,std::ios::binary);if(!output)return;
    File::Write_Binary(output,n*8);
    File::Write_Binary_Array(output,xf,n*8);
    delete [] xf;
}

template void Write_Segments_To_File_3d_Fast<2,real>(const Array<Vector2>&,const Array<Vector2>&,const std::string&,const bool);
template void Write_Segments_To_File_3d_Fast<3,real>(const Array<Vector3>&,const Array<Vector3>&,const std::string&,const bool);
template void Write_Segments_To_File_3d_Fast<2,float>(const Array<Vector2f>&,const Array<Vector2f>&,const std::string&,const bool);
template void Write_Segments_To_File_3d_Fast<3,float>(const Array<Vector3f>&,const Array<Vector3f>&,const std::string&,const bool);

template<int d,class T> void Write_Vectors_To_File_3d_Fast(const Array<Vector<T,d> >& X,const Array<Vector<T,d> >& V,const std::string& file_name)
{
    int n=(int)X.size();
    float* xf=new float[n*8];
    memset((void*)xf,0,sizeof(float)*n*8);
#pragma omp parallel for
    for(int i=0;i<n;i++){
        if constexpr (d==2){
            xf[i*8]=(float)X[i][0];
            xf[i*8+1]=(float)X[i][1];

            xf[i*8+4]=(float)(V[i][0]);
            xf[i*8+5]=(float)(V[i][1]);}

        else if constexpr (d==3){
            xf[i*8]=(float)X[i][0];
            xf[i*8+1]=(float)X[i][1];
            xf[i*8+2]=(float)X[i][2];

            xf[i*8+4]=(float)(V[i][0]);
            xf[i*8+5]=(float)(V[i][1]);
            xf[i*8+6]=(float)(V[i][2]);}}

    std::ofstream output(file_name,std::ios::binary);if(!output)return;
    File::Write_Binary(output,n*8);
    File::Write_Binary_Array(output,xf,n*8);
    delete [] xf;
}

template void Write_Vectors_To_File_3d_Fast<2,real>(const Array<Vector2>&,const Array<Vector2>&,const std::string&);
template void Write_Vectors_To_File_3d_Fast<3,real>(const Array<Vector3>&,const Array<Vector3>&,const std::string&);
template void Write_Vectors_To_File_3d_Fast<2,float>(const Array<Vector2f>&,const Array<Vector2f>&,const std::string&);
template void Write_Vectors_To_File_3d_Fast<3,float>(const Array<Vector3f>&,const Array<Vector3f>&,const std::string&);

