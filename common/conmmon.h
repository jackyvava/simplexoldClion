//
// Created by zhangjie on 2023/10/17.
//

#ifndef conmmon_h
#define conmmon_h
#include<Eigen/Dense>
#include<Eigen/Geometry>
#include<vector>
#include<list>
#include<queue>
#include<array>
#include<memory>

//使用新名字来简化Eigen命名
using Quaterniond = Eigen::Quaterniond;//四元数
using Quaternionf = Eigen::Quaternionf;
using Rotation2f = Eigen::Rotation2Df;//2D旋转矩阵
using Rotation2d = Eigen::Rotation2Dd;
using Transform2f = Eigen::Transform<float,2,Eigen::Affine>;
using Transform3f=Eigen::Transform<float,3,Eigen::Affine>;
using Transform2d=Eigen::Transform<double,2,Eigen::Affine>;
using Transform3d=Eigen::Transform<double,3,Eigen::Affine>;

#define Declare_Eigen_Types(my_type,t)  \
using real=my_type;                     \
using Vector1=Eigen::Matrix<real,1,1>;  \
using Vector2=Eigen::Vector2##t;        \
using Vector3=Eigen::Vector3##t;        \
using Vector4=Eigen::Vector4##t;        \
using Vector5=Eigen::Matrix<real,5,1>;; \
using Vector6=Eigen::Matrix<real,6,1>;; \
using Vector7=Eigen::Matrix<real,7,1>;; \
using Vector8=Eigen::Matrix<real,8,1>;; \
using Vector9=Eigen::Matrix<real,9,1>;; \
using VectorX=Eigen::VectorX##t;        \
using Matrix2=Eigen::Matrix2##t;        \
using Matrix3=Eigen::Matrix3##t;        \
using Matrix4=Eigen::Matrix4##t;        \
using MatrixX=Eigen::MatrixX##t;        \
using C=std::complex<real>;             \
using Vector1c=Eigen::Matrix<C,1,1>;    \
using VectorXc=Eigen::VectorXc##t;      \
using Vector2c=Eigen::Vector2c##t;      \
using Vector3c=Eigen::Vector3c##t;      \
using Vector4c=Eigen::Vector4c##t;      \
using Vector5c=Eigen::Matrix<C,5,1>;;	\
using Vector6c=Eigen::Matrix<C,6,1>;;	\
using Quaternion=Eigen::Quaternion##t;  \
using AngleAxis=Eigen::AngleAxis##t;	\
using Rotation2=Eigen::Rotation2D##t;	\
using Transform2=Transform2##t;	\
using Transform3=Transform3##t;

#define Declare_Eigen_Vector_Types(type,t)		\
using Vector1##t=Eigen::Matrix<type,1,1>;       \
using Vector2##t=Eigen::Vector2##t;             \
using Vector3##t=Eigen::Vector3##t;             \
using Vector4##t=Eigen::Vector4##t;             \
using VectorX##t=Eigen::VectorX##t;				\
using Vector5##t=Eigen::Matrix<type,5,1>;;		\
using Vector6##t=Eigen::Matrix<type,6,1>;;		\
using Vector7##t=Eigen::Matrix<type,7,1>;;		\
using Vector8##t=Eigen::Matrix<type,8,1>;;		\
using Vector9##t=Eigen::Matrix<type,9,1>;;

#define Declare_Eigen_Matrix_Types(type,t)		\
using Matrix1##t=Eigen::Matrix<type,1,1>;       \
using Matrix2##t=Eigen::Matrix2##t;             \
using Matrix3##t=Eigen::Matrix3##t;             \
using Matrix4##t=Eigen::Matrix4##t;             \
using MatrixX##t=Eigen::MatrixX##t;

#ifdef USE_FLOAT
Declare_Eigen_Types(float,f)
#else
Declare_Eigen_Types(double,d)
#endif
Declare_Eigen_Vector_Types(int,i)
Declare_Eigen_Vector_Types(float,f)
Declare_Eigen_Vector_Types(double,d)
Declare_Eigen_Matrix_Types(int,i)
Declare_Eigen_Matrix_Types(float,f)
Declare_Eigen_Matrix_Types(double,d)

using TI=int;
using uchar=unsigned char;
using ushort=unsigned short;
template<class T,int d> using Vector=Eigen::Matrix<T,d,1>;
template<class T,int d> using Matrix=Eigen::Matrix<T,d,d>;

template<class T> using Array=std::vector<T>;
template<class T> using ArrayPtr=std::shared_ptr<Array<T> >;
using size_type=Array<int>::size_type;

#endif //conmmon_h
