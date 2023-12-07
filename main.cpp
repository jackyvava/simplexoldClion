#include <Eigen/Core>
#include<Eigen/Dense>
#include<Eigen/Geometry>
#include <iostream>


using Transform2f=Eigen::Transform<float,2,Eigen::Affine>;
using Transform3f=Eigen::Transform<float,3,Eigen::Affine>;
using Transform2d=Eigen::Transform<double,2,Eigen::Affine>;
using Transform3d=Eigen::Transform<double,3,Eigen::Affine>;

int main() {
    // 创建一个二维浮点数仿射变换矩阵
    Transform2f transform2f;
    transform2f.setIdentity();  // 将矩阵初始化为单位矩阵

    // 创建一个三维双精度浮点数仿射变换矩阵
    Transform3d transform3d;
    transform3d.setIdentity();

    // 对二维矩阵进行平移操作
    transform2f.translate(Eigen::Vector2f(3.0, 3.0));

    // 对三维矩阵进行旋转操作
   Eigen::AngleAxisd rotation(M_PI / 4.0, Eigen::Vector3d::UnitZ());
   transform3d.rotate(rotation);

    // 打印二维矩阵和三维矩阵

    std::cout << "2D Transform:\n" << transform2f.matrix() << std::endl;
    std::cout << "3D Transform:\n" << transform3d.matrix() << std::endl;

    return 0;
}
