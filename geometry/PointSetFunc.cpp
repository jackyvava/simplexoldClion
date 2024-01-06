//
// Created by zhangjie on 2024/1/6.
//

#include "PointSetFunc.h"
#include <Eigen/StdVector>
#include <iostream>
using namespace AuxFunc;

real PointSetFunc::Initialize_Circle_Points(const Vector2& c, const real R, const int p_num, GeometryParticles<2>& particles)
{
    particles.Resize(p_num);
    for (int i = 0; i < p_num; i++) {
        particles.M(i) = (real)1;
        particles.F(i) = Vector2::Zero();
        particles.V(i) = Vector2::Zero();

        real theta = two_pi / (real)p_num * (real)i;
        Vector2 normal = V<2>(cos(theta), sin(theta));
        Vector2 pos = R * normal;
        particles.X(i) = pos + c;
        Vector2 tang = -Orthogonal_Vector(normal);
        particles.E(i).col(0) = tang;
        particles.E(i).col(1) = normal;
    }
    real dx = two_pi * R / (real)p_num;
    return dx;
}

int PointSetFunc::Initialize_Circle_Points(const Vector2& c, const real R, const real dx, GeometryParticles<2>& particles)
{
    int p_num = (int)ceil(two_pi * R / dx);
    Initialize_Circle_Points(c, R, p_num, particles);
    return p_num;
}

real PointSetFunc::Initialize_Oval_Points(const Vector2& c, const real R, const real a, const real b, const int p_num, GeometryParticles<2>& particles)
{
    particles.Resize(p_num);

    real angle = two_pi / (real)p_num;
    real total_length = 0;
    for (int i = 0; i < p_num; i++) {
        real theta = angle * ((real)i);
        Vector2 pos(a*cos(theta)*R, b*sin(theta)*R);
        particles.X(i) = pos+c;
        if (i != 0) {
            total_length+=(particles.X(i)-particles.X(i-1)).norm();
        }
    }
    real dx = total_length/p_num;
    return dx;
}

// Setup equidistanced points on a segment
// v1 and v2 are two ends, N is the number of points
void PointSetFunc::Initialize_Segment_Points(const Vector2& v1, const Vector2& v2, int N, GeometryParticles<2>& particles)
{
    particles.Resize(N);
    real h = (v2 - v1).norm() / N;

    // Normalized direction
    Vector2  u = (v2 - v1).normalized();
    Vector2  orthogonal_u = Orthogonal_Vector(u);
    particles.X(0) = v1 + h / 2 * u;

    for (int i = 0; i < N; i++) {
        particles.M(i) = (real)1;
        particles.F(i) = Vector2::Zero();
        particles.V(i) = Vector2::Zero();

        if (i != 0) {
            particles.X(i) = particles.X(0) + i * h * u;
        }

        particles.E(i).col(0) = u;
        particles.E(i).col(1) = orthogonal_u;
    }
}

void PointSetFunc::Initialize_Curve_Points(const Vector2& center, const real R, const real theta, const int N, GeometryParticles<2>& particles)
{
    particles.Resize(N+1);
    real angle=pi-(real)2*theta;real step=angle/(real)N;
    for(int i=0;i<N+1;i++){
        Vector2 normal=Vector2(cos(theta+step*(real)i),sin(theta+step*(real)i));
        particles.X(i)=(center+R*normal);
        particles.M(i)=(real)1;
        particles.F(i)=Vector2::Zero();
        particles.V(i)=Vector2::Zero();
        particles.E(i).col(0)=-Orthogonal_Vector(normal);
        particles.E(i).col(1)=normal;
    }
}

void PointSetFunc::Initialize_Rectangle_Points(const int nx, const int ny, const real dx, const Vector2& start, GeometryParticles<2>& particles)
{
    particles.Resize((nx + ny) * 2);

    for (int i = 0; i < nx; i++) {
        particles.X(i) = Vector2(start[0]+i*dx,start[1]);
    }

    for (int i = 0; i < ny; i++) {
        particles.X(i+nx) = Vector2(start[0] + nx * dx, start[1] + i * dx);
    }

    for (int i = 0; i < nx; i++) {
        particles.X(i + nx+ny) = Vector2(start[0] + nx * dx-i*dx, start[1] + ny * dx);
    }

    for (int i = 0; i < ny; i++) {
        particles.X(i + 2*nx + ny) = Vector2(start[0], start[1] + ny * dx-i*dx);
    }
}

void PointSetFunc::Initialize_Round_Corner_Rectangle_Points(const int nx, const int ny, const real dx, const real r, const Vector2& start, GeometryParticles<2>& particles)
{
    real left = start[0] - r;
    real bottom = start[1] - r;
    real right = start[0] + nx * dx + r;
    real top = start[1] + ny * dx + r;
    Vector2 rb_corner(start[0] + nx * dx,start[1]);
    Vector2 rt_corner(start[0] + nx * dx, start[1] + ny * dx);
    Vector2 lt_corner(start[0], start[1] + ny * dx);
    Vector2 lb_corner = start;

    int corner_num = r / dx * pi / 2;
    particles.Resize((nx + ny) * 2+corner_num*4);
    real theta = 2 * pi / corner_num / 4;

    int count = 0;
    for (int i = 0; i < nx; i++) {
        particles.X(i) = Vector2(start[0] + i * dx, bottom);
    }
    count += nx;

    for (int i = 0; i < corner_num; i++) {
        real angle = (real)1.5 *pi + i* theta;
        particles.X(count+i) = Vector2(cos(angle)*r,sin(angle)*r) + rb_corner;
    }
    count += corner_num;

    for (int i = 0; i < ny; i++) {
        particles.X(count+i) = Vector2(right, start[1] + i * dx);
    }
    count += ny;

    for (int i = 0; i < corner_num; i++) {
        real angle = i* theta;
        particles.X(count+i) = Vector2(cos(angle)*r, sin(angle)*r) + rt_corner;
    }
    count += corner_num;

    for (int i = 0; i < nx; i++) {
        particles.X(count+i) = Vector2(start[0] + nx * dx - i * dx, top);
    }
    count += nx;

    for (int i = 0; i < corner_num; i++) {
        real angle = (real)0.5 *pi + i * theta;
        particles.X(count+i) = Vector2(cos(angle)*r, sin(angle)*r) + lt_corner;
    }
    count += corner_num;

    for (int i = 0; i < ny; i++) {
        particles.X(count+i) = Vector2(left, start[1] + ny * dx - i * dx);
    }
    count += ny;

    for (int i = 0; i < corner_num; i++) {
        real angle = pi + i* theta;
        particles.X(i+count) = Vector2(cos(angle)*r, sin(angle)*r) + lb_corner;
    }
}

/*-----------------------------------------3D Functions-----------------------------------------------------------------*/

real PointSetFunc::Initialize_Sphere_Points(const Vector3& c, const real R, const int sub, GeometryParticles<3>& particles)
{
    TriangleMesh<3> mesh;
    MeshFunc::Initialize_Sphere_Mesh(R, &mesh, sub);
    MeshFunc::Update_Normals(mesh);
    particles.Resize((int)mesh.Vertices().size());
    for (int i = 0; i < particles.Size(); i++) {
        particles.X(i) = mesh.Vertices()[i] + c;
        particles.V(i) = Vector3::Zero();
        particles.F(i) = Vector3::Zero();
        particles.M(i) = (real)1;
        Vector3 normal = (*mesh.Normals())[i];
        Vector3 t1 = -Orthogonal_Vector(normal);
        Vector3 t2 = t1.cross(normal);
        particles.E(i).col(0) = t1;
        particles.E(i).col(1) = t2;
        particles.E(i).col(2) = normal;
    }
    real dx = MeshFunc::Average_Edge_Length<3>(mesh.Vertices(), mesh.Elements());
    return dx;
}

real PointSetFunc::Add_Sphere_Points(const Vector3& c, const real R, const int sub, GeometryParticles<3>& particles)
{
    TriangleMesh<3> mesh;
    MeshFunc::Initialize_Sphere_Mesh(R, &mesh, sub);
    MeshFunc::Update_Normals(mesh);
    int osize=particles.Size();
    particles.Resize(osize+(int)mesh.Vertices().size());
    for (int i = osize; i < particles.Size(); i++) {
        particles.X(i) = mesh.Vertices()[i-osize] + c;
        particles.V(i) = Vector3::Zero();
        particles.F(i) = Vector3::Zero();
        particles.M(i) = (real)1;
        Vector3 normal = (*mesh.Normals())[i-osize];
        Vector3 t1 = -Orthogonal_Vector(normal);
        Vector3 t2 = t1.cross(normal);
        particles.E(i).col(0) = t1;
        particles.E(i).col(1) = t2;
        particles.E(i).col(2) = normal;
    }
    real dx = MeshFunc::Average_Edge_Length<3>(mesh.Vertices(), mesh.Elements());
    return dx;
}

real PointSetFunc::Initialize_Circle_Points(const Vector3& c, const real R, const Vector3& normal, const real dx, GeometryParticles<3>& particles)
{
    int p_num = (int)ceil(two_pi * R / dx);
    TriangleMesh<3> mesh;
    MeshFunc::Initialize_Circle_Mesh(c, R, normal, &mesh, p_num);
    MeshFunc::Update_Normals(mesh);
    particles.Resize((int)mesh.Vertices().size());
    for (int i = 0; i < particles.Size(); i++) {
        particles.X(i) = mesh.Vertices()[i];
        particles.V(i) = Vector3::Zero();
        particles.F(i) = Vector3::Zero();
        particles.M(i) = (real)1;
        Vector3 normal = (*mesh.Normals())[i];
        Vector3 t1 = -Orthogonal_Vector(normal);
        Vector3 t2 = t1.cross(normal);
        particles.E(i).col(0) = t1;
        particles.E(i).col(1) = t2;
        particles.E(i).col(2) = normal;
    }
    real avg_dx = MeshFunc::Average_Edge_Length<3>(mesh.Vertices(), mesh.Elements());
    return avg_dx;
}

std::vector<int> PointSetFunc::Initialize_Circle_Points_Grid(const Vector3& c, const real R, const real dx, GeometryParticles<3>& particles)
{
    std::vector<Vector3d, Eigen::aligned_allocator<Vector3d>> pts;
    std::vector<int> is_boundary;
    int p_num = 2 * R / dx;
    for (int i = 0; i <= p_num; i++) {
        for (int j = 0; j <= p_num; j++) {
            real x_coord = 2 * R * (real)i / p_num - R;
            real z_coord = 2 * R * (real)j / p_num - R;
            Vector3d curr_x;
            curr_x << x_coord, 0, z_coord;
            if ((curr_x - Vector3d::Zero()).norm() < R-0.5*dx) {
                //if (1 == 1) {
                pts.push_back(curr_x);
                is_boundary.push_back(0);
            }
        }
    }

    int boundary_multiplier = 1;
    int num_boundary = boundary_multiplier * 2 * pi * R / dx;
    for (int i = 0; i < num_boundary; i++) {
        real x_coord = R * cos(((real)i/num_boundary)*(2*pi));
        real z_coord = R * sin(((real)i/num_boundary)*(2 * pi));
        Vector3d curr_x;
        curr_x << x_coord, 0, z_coord;
        pts.push_back(curr_x);
        is_boundary.push_back(1);
    }

    particles.Resize((int)pts.size());
    for (int i = 0; i < particles.Size(); i++) {
        particles.X(i) = pts[i] + c;
        particles.V(i) = Vector3::Zero();
        particles.F(i) = Vector3::Zero();
        particles.M(i) = (real)1;
        particles.E(i).col(0) = -Vector3d::Unit(2);
        particles.E(i).col(1) = -Vector3d::Unit(0);
        particles.E(i).col(2) = -Vector3d::Unit(1);
    }
    return is_boundary;
}


std::vector<int> PointSetFunc::Initialize_Circle_Points_Random(const Vector3& c, const real R, const real dx, GeometryParticles<3>& particles)
{
    std::vector<Vector3d, Eigen::aligned_allocator<Vector3d>> pts;
    std::vector<int> is_boundary;
    int num_points = 0;
    int p_num = 2 * R / dx; // num points along diameter
    for (int i = 0; i <= p_num; i++) {
        for (int j = 0; j <= p_num; j++) {
            real x_coord = 2 * R * (real)i / p_num - R;
            real z_coord = 2 * R * (real)j / p_num - R;
            Vector3d curr_x;
            curr_x << x_coord, 0, z_coord;
            if ((curr_x - Vector3d::Zero()).norm() < R - 0.5 * dx) {
                num_points++;
            }
        }
    }
    std::cout << "number of points: " <<num_points << std::endl;

    int num_obtained = 0;
    while (num_obtained < num_points) {
        real x_coord, z_coord;
        x_coord = ((real)rand() / (RAND_MAX+1.)) * (2*R)-R;
        z_coord = ((real)rand() / (RAND_MAX + 1.)) * (2*R)-R;
        //std::cout << "rand r: " << rand_r << std::endl;
        //std::cout << "rand theta: " << rand_theta << std::endl;
        Vector3d curr_x;
        curr_x << x_coord, 0, z_coord;
        if ((curr_x - Vector3d::Zero()).norm() < R - 0.5 * dx) {
            //if (1 == 1) {
            pts.push_back(curr_x);
            is_boundary.push_back(0);
            num_obtained++;
        }
    }

    int boundary_multiplier = 1; // how much denser are boundary points
    int num_boundary = boundary_multiplier * 2 * pi * R / dx;
    for (int i = 0; i < num_boundary; i++) {
        real x_coord = R * cos(((real)i / num_boundary) * (2 * pi));
        real z_coord = R * sin(((real)i / num_boundary) * (2 * pi));
        Vector3d curr_x;
        curr_x << x_coord, 0, z_coord;
        pts.push_back(curr_x);
        is_boundary.push_back(1);
    }

    particles.Resize((int)pts.size());
    for (int i = 0; i < particles.Size(); i++) {
        particles.X(i) = pts[i] + c;
        particles.V(i) = Vector3::Zero();
        particles.F(i) = Vector3::Zero();
        particles.M(i) = (real)1;
        particles.E(i).col(0) = -Vector3d::Unit(2);
        particles.E(i).col(1) = -Vector3d::Unit(0);
        particles.E(i).col(2) = -Vector3d::Unit(1);
    }
    return is_boundary;
}

std::vector<int> PointSetFunc::Initialize_Circle_Points_Grid_Random(const Vector3& c, const real R, const real dx, GeometryParticles<3>& particles)
{
    real randomness = 0.5 * dx;
    std::vector<Vector3d, Eigen::aligned_allocator<Vector3d>> pts;
    std::vector<int> is_boundary;
    int p_num = 2 * R / dx;
    for (int i = 0; i <= p_num; i++) {
        for (int j = 0; j <= p_num; j++) {
            real x_coord = 2 * R * (real)i / p_num - R;
            real z_coord = 2 * R * (real)j / p_num - R;
            x_coord += ((real)rand() / (RAND_MAX + 1.)) * (2 * randomness) - randomness;
            z_coord += ((real)rand() / (RAND_MAX + 1.)) * (2 * randomness) - randomness;
            Vector3d curr_x;
            curr_x << x_coord, 0, z_coord;
            if ((curr_x - Vector3d::Zero()).norm() < R - 0.5 * dx) {
                //if (1 == 1) {
                pts.push_back(curr_x);
                is_boundary.push_back(0);
            }
        }
    }

    int boundary_multiplier = 1;
    int num_boundary = boundary_multiplier * 2 * pi * R / dx;
    for (int i = 0; i < num_boundary; i++) {
        real x_coord = R * cos(((real)i / num_boundary) * (2 * pi));
        real z_coord = R * sin(((real)i / num_boundary) * (2 * pi));
        Vector3d curr_x;
        curr_x << x_coord, 0, z_coord;
        pts.push_back(curr_x);
        is_boundary.push_back(1);
    }

    particles.Resize((int)pts.size());
    for (int i = 0; i < particles.Size(); i++) {
        particles.X(i) = pts[i] + c;
        particles.V(i) = Vector3::Zero();
        particles.F(i) = Vector3::Zero();
        particles.M(i) = (real)1;
        particles.E(i).col(0) = -Vector3d::Unit(2);
        particles.E(i).col(1) = -Vector3d::Unit(0);
        particles.E(i).col(2) = -Vector3d::Unit(1);
    }
    return is_boundary;
}


std::vector<int> PointSetFunc::Initialize_Catenoid_Points(const Vector3& c, const real R, const real dx, GeometryParticles<3>& particles)
{
    real randomness = 0.5 * dx;
    std::vector<Vector3d, Eigen::aligned_allocator<Vector3d>> pts;
    std::vector<int> is_boundary;
    real init_height = R / 1.5;

    int num_boundary = 2 * pi * R / dx;
    for (int i = 0; i < num_boundary; i++) {
        real x_coord = R * cos(((real)i / num_boundary) * (2 * pi));
        real z_coord = R * sin(((real)i / num_boundary) * (2 * pi));
        Vector3d curr_x;
        curr_x << x_coord, 0, z_coord;
        pts.push_back(curr_x);
        is_boundary.push_back(1);
        curr_x << x_coord, init_height, z_coord;
        pts.push_back(curr_x);
        is_boundary.push_back(2);
        for (int j = 1; j <= (int)((init_height - 2 * dx) / dx); j++) {
            curr_x << x_coord, j*dx + 1 * dx, z_coord;
            pts.push_back(curr_x);
            is_boundary.push_back(0);
        }
    }

    particles.Resize((int)pts.size());
    for (int i = 0; i < particles.Size(); i++) {
        particles.X(i) = pts[i] + c;
        particles.V(i) = Vector3::Zero();
        particles.F(i) = Vector3::Zero();
        particles.M(i) = (real)1;

        Vector3 k1, k2;
        k2 << particles.X(i)[0], 0, particles.X(i)[2];
        k2.normalize();
        k1 = Vector3::Unit(1);
        Vector3 k3 = k1.cross(k2);
        k2 = k1.cross(k3);
        particles.E(i).col(0) = k1;
        particles.E(i).col(1) = k3;
        particles.E(i).col(2) = k2;
    }
    return is_boundary;
}

std::vector<int> PointSetFunc::Initialize_Catenoid_Points_Random(const Vector3& c, const real R, const real dx, GeometryParticles<3>& particles)
{
    real randomness = 0.5 * dx;
    std::vector<Vector3d, Eigen::aligned_allocator<Vector3d>> pts;
    std::vector<int> is_boundary;
    real init_height = R / 1.5;

    int num_boundary = 2 * pi * R / dx;
    for (int i = 0; i < num_boundary; i++) {
        real x_coord = R * cos(((real)i / num_boundary) * (2 * pi));
        real z_coord = R * sin(((real)i / num_boundary) * (2 * pi));
        Vector3d curr_x;
        curr_x << x_coord, 0, z_coord;
        pts.push_back(curr_x);
        is_boundary.push_back(1);
        curr_x << x_coord, init_height, z_coord;
        pts.push_back(curr_x);
        is_boundary.push_back(2);
    }

    int num_fluid_particles = 8 * pts.size();

    for (int i = 0; i < num_fluid_particles; i++) {
        real angle = ((real)rand() / (RAND_MAX + 1.)) * 2 * pi;
        real height = ((real)rand() / (RAND_MAX + 1.)) * (init_height - 1. * dx);
        Vector3d curr_x;
        curr_x << R * cos(angle), height + 0.5*dx, R* sin(angle);
        pts.push_back(curr_x);
        is_boundary.push_back(0);
    }

    particles.Resize((int)pts.size());
    for (int i = 0; i < particles.Size(); i++) {
        particles.X(i) = pts[i] + c;
        particles.V(i) = Vector3::Zero();
        particles.F(i) = Vector3::Zero();
        particles.M(i) = (real)1;

        Vector3 k1, k2;
        k2 << particles.X(i)[0], 0, particles.X(i)[2];
        k2.normalize();
        k1 = Vector3::Unit(1);
        Vector3 k3 = k1.cross(k2);
        k2 = k1.cross(k3);
        particles.E(i).col(0) = k1;
        particles.E(i).col(1) = k3;
        particles.E(i).col(2) = k2;
    }
    return is_boundary;
}

void PointSetFunc::Initialize_Rectangle_Points(const int nx, const int ny, const real dx, const Vector3& start, const Vector3& normal, GeometryParticles<3>& particles)
{
    ////kxx: TODO
}

void PointSetFunc::Initialize_Lattice_Points(const Vector3& domain_min, const Vector2i& counts, Vector3 k1, Vector3 k2, const real dx, GeometryParticles<3>& particles) {
    k1.normalize();
    k2.normalize();
    Vector3 k3 = k1.cross(k2);
    int n1 = counts[0], n2 = counts[1];
    particles.Resize((n1 + 1) * (n2 + 1));
    int idx = 0;
    for (int i = 0; i <= n1; i++) {
        for (int j = 0; j <= n2; j++) {
            Vector3 pos = domain_min + k1 * i * dx + k2 * j * dx;
            particles.X(idx) = pos;
            particles.V(idx) = Vector3::Zero();
            particles.F(idx) = Vector3::Zero();
            particles.M(idx) = 1.0;
            particles.E(idx).col(0) = k1;
            particles.E(idx).col(1) = k2;
            particles.E(idx).col(2) = k3;
            idx++;
        }
    }
}

//Yitong's modification
std::vector<int> PointSetFunc::Initialize_Lattice_Points2(const Vector3& domain_min, const Vector2i& counts, Vector3 k1, Vector3 k2, const real dx, GeometryParticles<3>& particles) {
    std::vector<int> is_boundary;
    int num_padding = 5;
    k1.normalize();
    k2.normalize();
    Vector3 k3 = k1.cross(k2);
    int n1 = counts[0], n2 = counts[1];
    int x_skip_start = (int)(n1 / 2);
    int x_skip_end = n1;
    real R = 0.3;
    int num_boundary = 2 * pi * R / dx;

    particles.Resize((n1 + 2 * num_padding + 1) * (n2 + 2 * num_padding + 1) - (n2-1) * (x_skip_end - x_skip_start+1) +num_boundary);
    //std::cout << "particles size" << particles.Size() << std::endl;
    int idx = 0;
    for (int i = -num_padding; i <= n1+ num_padding; i++) {
        for (int j = -num_padding; j <= n2+ num_padding; j++) {
            if (j > 0 && j < n2 && i >= x_skip_start && i <= x_skip_end)continue;
            Vector3 pos = domain_min + k1 * i * dx + k3 * j * dx;
            particles.X(idx) = pos;
            particles.V(idx) = Vector3::Zero();
            particles.F(idx) = Vector3::Zero();
            particles.M(idx) = 1.0;
            particles.E(idx).col(0) = -k3;
            particles.E(idx).col(1) = -k1;
            particles.E(idx).col(2) = -k2;
            if (i <= 0 || i >= n1 || j <= 0 || j >= n2) {
                is_boundary.push_back(1);
            }
            else {
                is_boundary.push_back(0);
            }
            idx++;
        }
    }
    std::cout << "num boundary" << num_boundary << std::endl;
    for (int i = 0; i < num_boundary; i++) {
        //std::cout << "i" << i << std::endl;
        real x_coord = R * cos(((real)i / num_boundary) * (2 * pi));
        real z_coord = R * sin(((real)i / num_boundary) * (2 * pi));
        Vector3d curr_x;
        curr_x << x_coord + 0.5, 0, z_coord;
        //std::cout << "curr x\n" << curr_x << std::endl;
        particles.X(idx) = curr_x;
        particles.M(idx) = 1.0;
        particles.E(idx).col(0) = -k3;
        particles.E(idx).col(1) = -k1;
        particles.E(idx).col(2) = -k2;
        is_boundary.push_back(1);
        idx++;
        //std::cout << "idx" << idx << std::endl;
    }

    return is_boundary;
}

real PointSetFunc::Initialize_Ellipsoid_Points(const Vector3& C, const real R, const int sub, const real a, const real b, const real c,  GeometryParticles<3>& particles)
{
    TriangleMesh<3> mesh;
    MeshFunc::Initialize_Ellipsoid_Mesh(R, &mesh, a,b,c, sub);
    MeshFunc::Update_Normals(mesh);
    particles.Resize((int)mesh.Vertices().size());
    for (int i = 0; i < particles.Size(); i++) {
        particles.X(i) = mesh.Vertices()[i] + C;
        particles.V(i) = Vector3::Zero();
        particles.F(i) = Vector3::Zero();
        particles.M(i) = (real)1;
        Vector3 normal = (*mesh.Normals())[i];
        Vector3 t1 = -Orthogonal_Vector(normal);
        Vector3 t2 = t1.cross(normal);
        particles.E(i).col(0) = t1;
        particles.E(i).col(1) = t2;
        particles.E(i).col(2) = normal;
    }
    real dx = MeshFunc::Average_Edge_Length<3>(mesh.Vertices(), mesh.Elements());
    return dx;
}

std::vector<int> PointSetFunc::Initialize_Catenoid_Points(const Vector3& center, const real R, const int nr, const real height, const int nh, GeometryParticles<3>& particles)
{
    assert(nh > 1);
    int total_size = nr * nh;
    Array<int> bd_flags;
    bd_flags.resize(total_size);
    particles.Resize(total_size);
    real y_min = -height / 2, dh = height / (nh - 1.0);
    int idx = 0;
    for (int i = 0; i < nh; i++) {
        real p_y = y_min + dh * (i + 0.0);
        //real w = 1.0 + (p_y * p_y - height * height * 0.25) / (height * height * 0.25) * 0.1;
        real w = cosh(p_y / height) * height;
        for (int j = 0; j < nr; j++) {
            real alpha = 2.0 * pi / (nr + 0.0) * (j + 0.0);
            real r = R * w;
            Vector3 pos(r * cos(alpha), p_y, r * sin(alpha));
            particles.X(idx) = pos + center;
            particles.V(idx) = Vector3::Zero();
            particles.F(idx) = Vector3::Zero();
            particles.M(idx) = 0.0;
            Vector3 normal = pos; normal[1] = 0; normal.normalize();
            Vector3 t1 = -Orthogonal_Vector(normal);
            Vector3 t2 = t1.cross(normal);
            particles.E(idx).col(0) = t1;
            particles.E(idx).col(1) = t2;
            particles.E(idx).col(2) = normal;
            if (i == 0 || i == nh - 1) bd_flags[idx] = 1;
            else bd_flags[idx] = 0;
            idx++;
        }
    }
    return bd_flags;
}

