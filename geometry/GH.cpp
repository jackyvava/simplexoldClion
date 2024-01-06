#include "GH.h"

namespace GH{
    template<class T>
    inline T Min(T a1, T a2, T a3){
        return std::min(a1, std::min(a2, a3));
    }

    template<class T>
    inline T Max(T a1, T a2, T a3){
        return std::max(a1, std::max(a2, a3));
    }

    template<int d>
    //这里做了修改，去除了<d>
    double Point_Segment_Distance(const Vec<d>& x0, const Vec<d>& x1, const Vec<d>& x2){
        Vec<d> dx = x2 - x1;
        double m2 = dx.dot(dx);
        // find parameter value of closest point on segment
        double s12 = dx.dot(x2-x0)/m2;
        if(s12<0){
            s12=0;
        }else if(s12>1){
            s12=1;
        }
        // and find the distance
        return (x0-(s12*x1+(1-s12)*x2)).norm();
    }

    template double Point_Segment_Distance<2>(const Vec<2>& x0, const Vec<2>& x1, const Vec<2>& x2);
    template double Point_Segment_Distance<3>(const Vec<3>& x0, const Vec<3>& x1, const Vec<3>& x2);

    double Point_Triangle_Distance(const Vec<3>& x0, const Vec<3>& x1, const Vec<3>& x2, const Vec<3>& x3){
        // first find barycentric coordinates of closest point on infinite plane
        Vec<3> x13=x1-x3, x23=x2-x3, x03=x0-x3;
        double m13=x13.dot(x13), m23=x23.dot(x23), d=x13.dot(x23);
        double invdet=1.0/std::max(m13*m23-d*d,1e-30);
        double a=x13.dot(x03), b=x23.dot(x03);
        // the barycentric coordinates themselves
        double w23=invdet*(m23*a-d*b);
        double w31=invdet*(m13*b-d*a);
        double w12=1-w23-w31;
        if(w23>=0 && w31>=0 && w12>=0){ // if we're inside the triangle
            return (x0-(w23*x1+w31*x2+w12*x3)).norm();
        }else{ // we have to clamp to one of the edges
            if(w23>0) // this rules out edge 2-3 for us
                return std::min(Point_Segment_Distance<3>(x0,x1,x2), Point_Segment_Distance<3>(x0,x1,x3));
            else if(w31>0) // this rules out edge 1-3
                return std::min(Point_Segment_Distance<3>(x0,x1,x2), Point_Segment_Distance<3>(x0,x2,x3));
            else // w12 must be >0, ruling out edge 1-2
                return std::min(Point_Segment_Distance<3>(x0,x1,x3), Point_Segment_Distance<3>(x0,x2,x3));
        }
    }

    bool Point_In_Triangle_2D(const Vec<2>& x0, const Vec<2>& x1, const Vec<2>& x2, const Vec<2>& x3,
                              double& a, double& b, double& c){
        Vec<2> x10 = x1 - x0;
        Vec<2> x20 = x2 - x0;
        Vec<2> x30 = x3 - x0;

        //return the side of origin to (x1,x2), and also set twice_signed_area for barycentric cooridnate
        auto orientation = [] (const Vec<2>& tx1, const Vec<2>& tx2, double &twice_signed_area) -> int {
            twice_signed_area = tx1(1)*tx2(0)-tx1(0)*tx2(1);
            if(twice_signed_area > 0) return 1;
            else if(twice_signed_area < 0) return -1;
            else if(tx2(1) > tx1(1)) return 1;
            else if(tx2(1) < tx1(1)) return -1;
            else if(tx1(0) > tx2(0)) return 1;
            else if(tx1(0) < tx2(0)) return -1;
            else return 0;
        };

        int signa = orientation(x20, x30, a);
        if(signa == 0) return false;
        int signb = orientation(x30, x10, b);
        if(signb != signa) return false;
        int signc = orientation(x10, x20, c);
        if(signc != signa) return false;
        double sum = a + b + c;
        a /= sum;
        b /= sum;
        c /= sum;
        return true;
    }

    bool Segment_Box_Overlap_2D(const Vec<2>& box_center, const Vec<2>& box_half_size, const Vec<2>& x1, const Vec<2>& x2){
        int count = 0;
        Vec<2> v1 = x1 - box_center;
        Vec<2> v2 = x2 - box_center;
        auto F = [&](Vec<2> p) -> int {
            double tF = (v2(1)-v1(1))*p(0)+(v1(0)-v2(0))*p(1)+v2(0)*v1(1)-v1(0)*v2(1);
            return int(0 < tF) - int(tF < 0);
        };
        count += F(box_half_size);
        count += F(-box_half_size);
        count += F(Vec<2>(box_half_size(0),-box_half_size(1)));
        count += F(Vec<2>(-box_half_size(0),box_half_size(1)));
        if((count==-4)||(count==4))return false;

        if(v1(0)>box_half_size(0) && v2(0)>box_half_size(0))return false;
        if(v1(0)<-box_half_size(0) && v2(0)<-box_half_size(0))return false;
        if(v1(1)>box_half_size(1) && v2(1)>box_half_size(1))return false;
        if(v1(1)<-box_half_size(1) && v2(1)<-box_half_size(1))return false;

        return true;
    }

    bool Plane_Box_Overlap_3D(const Vec<3>& normal, const Vec<3>& vert, const Vec<3>& box_half_size){
        int q;double v;
        Vec<3> vmin, vmax;
        for(q=0;q<3;++q){
            v = vert(q);
            if(normal(q)>0.0){
                vmin(q) = -box_half_size(q)-v;
                vmax(q) = box_half_size(q)-v;
            }
            else{
                vmin(q) = box_half_size(q)-v;
                vmax(q) = -box_half_size(q)-v;
            }
        }
        if(normal.dot(vmin)>0.0) return false;
        if(normal.dot(vmax)>=0.0) return true;
        return false;
    }

    bool Triangle_Box_Overlap_3D(const Vec<3>& box_center, const Vec<3>& box_half_size, const Vec<3>& x1, const Vec<3>& x2, const Vec<3>& x3){
        double min,max,p0,p1,p2,rad,fex,fey,fez;
        Vec<3> v0,v1,v2,e0,e1,e2,normal;
        //centralized vertices
        v0 = x1-box_center;
        v1 = x2-box_center;
        v2 = x3-box_center;
        //edge
        e0 = v1-v0;
        e1 = v2-v1;
        e2 = v0-v2;

        //test macros
        //X-tests
        auto Axis_Test_X01 = [&](const double a, const double b, const double fa, const double fb) -> bool {
            p0 = a*v0(1)-b*v0(2);
            p2 = a*v2(1)-b*v2(2);
            if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;}
            rad = fa*box_half_size(1)+fb*box_half_size(2);
            if((min>rad) || (max<-rad)) return false;
            else return true;
        };

        auto Axis_Test_X2 = [&](const double a, const double b, const double fa, const double fb) -> bool {
            p0 = a*v0(1)-b*v0(2);
            p1 = a*v1(1)-b*v1(2);
            if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;}
            rad = fa*box_half_size(1)+fb*box_half_size(2);
            if((min>rad) || (max<-rad)) return false;
            else return true;
        };

        //Y-tests
        auto Axis_Test_Y02 = [&](const double a, const double b, const double fa, const double fb) -> bool {
            p0 = -a*v0(0)+b*v0(2);
            p2 = -a*v2(0)+b*v2(2);
            if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;}
            rad = fa*box_half_size(0)+fb*box_half_size(2);
            if((min>rad) || (max<-rad)) return false;
            else return true;
        };

        auto Axis_Test_Y1 = [&](const double a, const double b, const double fa, const double fb) -> bool {
            p0 = -a*v0(0)+b*v0(2);
            p1 = -a*v1(0)+b*v1(2);
            if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;}
            rad = fa*box_half_size(0)+fb*box_half_size(2);
            if((min>rad) || (max<-rad)) return false;
            else return true;
        };

        //Z-tests
        auto Axis_Test_Z12 = [&](const double a, const double b, const double fa, const double fb) -> bool {
            p1 = a*v1(0)-b*v1(1);
            p2 = a*v2(0)-b*v2(1);
            if(p2<p1) {min=p2; max=p1;} else {min=p1; max=p2;}
            rad = fa*box_half_size(0)+fb*box_half_size(1);
            if((min>rad) || (max<-rad)) return false;
            else return true;
        };

        auto Axis_Test_Z0 = [&](const double a, const double b, const double fa, const double fb) -> bool {
            p0 = a*v0(0)-b*v0(1);
            p1 = a*v1(0)-b*v1(1);
            if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;}
            rad = fa*box_half_size(0)+fb*box_half_size(1);
            if((min>rad) || (max<-rad)) return false;
            else return true;
        };

        fex = std::abs(e0(0));
        fey = std::abs(e0(1));
        fez = std::abs(e0(2));
        if(!Axis_Test_X01(e0(2), e0(1), fez, fey))return false;
        if(!Axis_Test_Y02(e0(2), e0(0), fez, fex))return false;
        if(!Axis_Test_Z12(e0(1), e0(0), fey, fex))return false;

        fex = std::abs(e1(0));
        fey = std::abs(e1(1));
        fez = std::abs(e1(2));
        if(!Axis_Test_X01(e1(2), e1(1), fez, fey))return false;
        if(!Axis_Test_Y02(e1(2), e1(0), fez, fex))return false;
        if(!Axis_Test_Z0(e1(1), e1(0), fey, fex))return false;

        fex = std::abs(e2(0));
        fey = std::abs(e2(1));
        fez = std::abs(e2(2));
        if(!Axis_Test_X2(e2(2), e2(1), fez, fey))return false;
        if(!Axis_Test_Y1(e2(2), e2(0), fez, fex))return false;
        if(!Axis_Test_Z12(e2(1), e2(0), fey, fex))return false;

        min = Min(v0(0), v1(0), v2(0)); max = Max(v0(0), v1(0), v2(0));
        if((min>box_half_size(0)) || (max<-box_half_size(0))) return false;

        min = Min(v0(1), v1(1), v2(1)); max = Max(v0(1), v1(1), v2(1));
        if((min>box_half_size(1)) || (max<-box_half_size(1))) return false;

        min = Min(v0(2), v1(2), v2(2)); max = Max(v0(2), v1(2), v2(2));
        if((min>box_half_size(2)) || (max<-box_half_size(2))) return false;

        //plane box intersection
        normal = e0.cross(e1);
        if(!Plane_Box_Overlap_3D(normal, v0, box_half_size)) return false;
        return true;
    }

    template<>
    void Mesh_To_SDF<2>(const SurfaceMesh<2>& mesh, const Grid<2>& grid, std::vector<double>& phi, int narrow_band){
        //initialize phi to infinity
        double init_value = grid.cell_counts.sum()*grid.dx;
        phi.resize(grid.cell_counts.prod());
        for(int i=0;i<phi.size();++i)
            phi[i] = init_value;
        //trace the closest segment of each cell
        std::vector<int> closest_seg;closest_seg.resize(grid.cell_counts.prod(), -1);
        //trace the ray intersection at cell(i,j), the ray start from (0,j) in ith direction
        std::vector<int> intersection_count;intersection_count.resize(grid.cell_counts.prod(), 0);


        const std::vector<Vec<2> >& x = mesh.Vertices();
        for(int s=0;s<mesh.elements.size();++s){
            int p,q;
            p = mesh.elements[s](0);
            q = mesh.elements[s](1);
            Veci<2> fp = grid.Cell_Coord(x[p]);
            Veci<2> fq = grid.Cell_Coord(x[q]);
            int i0 = std::clamp(std::min(fp(0), fq(0))-narrow_band, 0, grid.cell_counts(0)-1);
            int i1 = std::clamp(std::max(fp(0), fq(0))+narrow_band, 0, grid.cell_counts(0)-1);
            int j0 = std::clamp(std::min(fp(1), fq(1))-narrow_band, 0, grid.cell_counts(1)-1);
            int j1 = std::clamp(std::max(fp(1), fq(1))+narrow_band, 0, grid.cell_counts(1)-1);

            for(int j=j0;j<=j1;++j){
                for(int i=i0;i<=i1;++i){
                    Vec<2> gx = grid.Center(Veci<2>(i,j));
                    double dist = Point_Segment_Distance<2>(gx, x[p], x[q]);
                    int idx = grid.Cell_Index(Veci<2>(i,j));
                    if(dist < phi[idx]){
                        phi[idx] = dist;
                        closest_seg[idx] = s;
                    }
                }
            }

            // intersection count
            j0 = std::clamp(std::min(fp(1),fq(1)), 0, grid.cell_counts(1)-1);
            j1 = std::clamp(std::max(fp(1),fq(1)), 0, grid.cell_counts(1)-1);
            for(int j=j0; j<=j1; ++j){
                double a;
                double pj = grid.Center(Veci<2>(0,j))(1);//get j's position
                if((pj>=x[p](1) && pj<x[q](1)) || (pj>=x[q](1) && pj<x[p](1))){
                    a = std::fabs((pj-x[q](1))/(x[p](1)-x[q](1)));
                    a = std::clamp(a, 0.0, 1.0);
                    double pi = a*x[p](0) + (1-a)*x[q](0);//get intersection point's ith coordinate
                    if(pi > grid.domain_min(0) && pi < grid.domain_max(0)){
                        Veci<2> icell = grid.Cell_Coord(Vec<2>(pi, pj));
                        if(pi > grid.Center(icell)(0)){
                            //intersection happened after the center of icell
                            icell(0) += 1;//change the sign of next cell
                        }
                        if(!grid.Valid_Cell(icell))continue;
                        int idx = grid.Cell_Index(icell);
                        ++intersection_count[idx];//ray along x coordinate from (0,j)'s center intersect with this tri at icell
                    }
                }
            }
        }

        //update i0's phi using i1
        auto check_neighbor = [&] (const Veci<2>& i0, const Veci<2>& i1) {
            if(closest_seg[grid.Cell_Index(i1)] >= 0){
                int cseg = closest_seg[grid.Cell_Index(i1)];
                int tp = mesh.elements[cseg](0);
                int tq = mesh.elements[cseg](1);
                Vec<2> gx = grid.Center(i0);
                double td = Point_Segment_Distance<2>(gx, x[tp], x[tq]);
                if(td < phi[grid.Cell_Index(i0)]){
                    phi[grid.Cell_Index(i0)] = td;
                    closest_seg[grid.Cell_Index(i0)] = closest_seg[grid.Cell_Index(i1)];
                }
            }
        };

        //sweep in (di,dj,dk) direction
        auto sweep = [&] (int di, int dj) {
            int ti0, ti1;
            if(di>0) {ti0=1;ti1=grid.cell_counts(0);}
            else {ti0=grid.cell_counts(0)-2;ti1=-1;}
            int tj0, tj1;
            if(dj>0) {tj0=1;tj1=grid.cell_counts(1);}
            else {tj0=grid.cell_counts(1)-2;tj1=-1;}
            for(int tj=tj0; tj!=tj1; tj+=dj){
                for(int ti=ti0; ti!=ti1; ti+=di){
                    check_neighbor(Veci<2>(ti,tj), Veci<2>(ti-di,tj));
                    check_neighbor(Veci<2>(ti,tj), Veci<2>(ti,tj-dj));
                    check_neighbor(Veci<2>(ti,tj), Veci<2>(ti-di,tj-dj));
                }
            }
        };

        //sweep in all direction
        for(int pass=0; pass<2; ++pass){
            sweep(+1,+1);
            sweep(-1,-1);
            sweep(+1,-1);
            sweep(-1,+1);
        }

        //determine sign from intersection
        for(int j=0; j<grid.cell_counts(1); ++j){
            int total_count=0;
            for(int i=0; i<grid.cell_counts(0); ++i){
                int idx = grid.Cell_Index(Veci<2>(i,j));
                total_count += intersection_count[idx];
                if(total_count%2 == 1){
                    phi[idx] = -phi[idx];
                }
            }
        }
    }

    template<>
    void Mesh_To_SDF<3>(const SurfaceMesh<3>& mesh, const Grid<3>& grid, std::vector<double>& phi, int narrow_band){
        //initialize phi to infinity
        double init_value = grid.cell_counts.sum()*grid.dx;
        phi.resize(grid.cell_counts.prod());
        for(int i=0;i<phi.size();++i)
            phi[i] = init_value;
        //trace the closest triangle of each cell
        std::vector<int> closest_tri;closest_tri.resize(grid.cell_counts.prod(), -1);
        //trace the ray intersection at cell(i,j,k), the ray start from (0,j,k) in ith direction
        std::vector<int> intersection_count;intersection_count.resize(grid.cell_counts.prod(), 0);


        const std::vector<Vec<3> >& x = mesh.Vertices();
        for(int t=0;t<mesh.elements.size();++t){
            int p,q,r;
            p = mesh.elements[t](0);
            q = mesh.elements[t](1);
            r = mesh.elements[t](2);
            Veci<3> fp = grid.Cell_Coord(x[p]);
            Veci<3> fq = grid.Cell_Coord(x[q]);
            Veci<3> fr = grid.Cell_Coord(x[r]);
            int i0 = std::clamp(Min(fp(0), fq(0), fr(0))-narrow_band, 0, grid.cell_counts(0)-1);
            int i1 = std::clamp(Max(fp(0), fq(0), fr(0))+narrow_band, 0, grid.cell_counts(0)-1);
            int j0 = std::clamp(Min(fp(1), fq(1), fr(1))-narrow_band, 0, grid.cell_counts(1)-1);
            int j1 = std::clamp(Max(fp(1), fq(1), fr(1))+narrow_band, 0, grid.cell_counts(1)-1);
            int k0 = std::clamp(Min(fp(2), fq(2), fr(2))-narrow_band, 0, grid.cell_counts(2)-1);
            int k1 = std::clamp(Max(fp(2), fq(2), fr(2))+narrow_band, 0, grid.cell_counts(2)-1);

            for(int k=k0;k<=k1;++k){
                for(int j=j0;j<=j1;++j){
                    for(int i=i0;i<=i1;++i){
                        Vec<3> gx = grid.Center(Veci<3>(i,j,k));
                        double dist = Point_Triangle_Distance(gx, x[p], x[q], x[r]);
                        int idx = grid.Cell_Index(Veci<3>(i,j,k));
                        if(dist < phi[idx]){
                            phi[idx] = dist;
                            closest_tri[idx] = t;
                        }
                    }
                }
            }

            // intersection count
            j0 = std::clamp(Min(fp(1),fq(1),fr(1)), 0, grid.cell_counts(1)-1);
            j1 = std::clamp(Max(fp(1),fq(1),fr(1)), 0, grid.cell_counts(1)-1);
            k0 = std::clamp(Min(fp(2),fq(2),fr(2)), 0, grid.cell_counts(2)-1);
            k1 = std::clamp(Max(fp(2),fq(2),fr(2)), 0, grid.cell_counts(2)-1);
            for(int k=k0; k<=k1; ++k){
                for(int j=j0; j<=j1; ++j){
                    double a, b, c;
                    Vec<3> ljk = grid.Center(Veci<3>(0,j,k));//get j,k's position
                    if(Point_In_Triangle_2D(ljk.segment<2>(1), x[p].segment<2>(1), x[q].segment<2>(1), x[r].segment<2>(1), a, b, c)){
                        double li = a*x[p](0) + b*x[q](0) + c*x[r](0);//get intersection point's ith coordinate
                        if(li > grid.domain_min(0) && li < grid.domain_max(0)){
                            Veci<3> icell = grid.Cell_Coord(Vec<3>(li,ljk(1),ljk(2)));
                            int idx = grid.Cell_Index(icell);
                            if(li > grid.Center(icell)(0)){
                                //intersection happened after the center of icell
                                icell(0) += 1;//change the sign of next cell
                            }
                            if(!grid.Valid_Cell(icell))continue;
                            ++intersection_count[idx];//ray along x coordinate from (0,j,k)'s center intersect with this tri at icell
                        }
                    }
                }
            }
        }

        //update i0's phi using i1
        auto check_neighbor = [&] (const Veci<3>& i0, const Veci<3>& i1) {
            if(closest_tri[grid.Cell_Index(i1)] >= 0){
                int ctri = closest_tri[grid.Cell_Index(i1)];
                int tp = mesh.elements[ctri](0);
                int tq = mesh.elements[ctri](1);
                int tr = mesh.elements[ctri](2);
                Vec<3> gx = grid.Center(i0);
                double td = Point_Triangle_Distance(gx, x[tp], x[tq], x[tr]);
                if(td < phi[grid.Cell_Index(i0)]){
                    phi[grid.Cell_Index(i0)] = td;
                    closest_tri[grid.Cell_Index(i0)] = closest_tri[grid.Cell_Index(i1)];
                }
            }
        };

        //sweep in (di,dj,dk) direction
        auto sweep = [&] (int di, int dj, int dk) {
            int ti0, ti1;
            if(di>0) {ti0=1;ti1=grid.cell_counts(0);}
            else {ti0=grid.cell_counts(0)-2;ti1=-1;}
            int tj0, tj1;
            if(dj>0) {tj0=1;tj1=grid.cell_counts(1);}
            else {tj0=grid.cell_counts(1)-2;tj1=-1;}
            int tk0, tk1;
            if(dk>0) {tk0=1;tk1=grid.cell_counts(2);}
            else {tk0=grid.cell_counts(2)-2;tk1=-1;}
            for(int tk=tk0; tk!=tk1; tk+=dk){
                for(int tj=tj0; tj!=tj1; tj+=dj){
                    for(int ti=ti0; ti!=ti1; ti+=di){
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti-di,tj,tk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti,tj-dj,tk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti-di,tj-dj,tk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti,tj,tk-dk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti-di,tj,tk-dk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti,tj-dj,tk-dk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti-di,tj-dj,tk-dk));
                    }
                }
            }
        };

        //sweep in all direction
        for(int pass=0; pass<2; ++pass){
            sweep(+1,+1,+1);
            sweep(-1,-1,-1);
            sweep(+1,+1,-1);
            sweep(-1,-1,+1);
            sweep(+1,-1,+1);
            sweep(-1,+1,-1);
            sweep(+1,-1,-1);
            sweep(-1,+1,+1);
        }

        //determine sign from intersection
        for(int k=0;k<grid.cell_counts(2);++k){
            for(int j=0;j<grid.cell_counts(1);++j){
                int total_count=0;
                for(int i=0; i<grid.cell_counts(0); ++i){
                    int idx = grid.Cell_Index(Veci<3>(i,j,k));
                    total_count += intersection_count[idx];
                    if(total_count%2 == 1){
                        phi[idx] = -phi[idx];
                    }
                }
            }
        }
    }

    template<>
    void Mesh_To_SDF<2>(const Codimension::CodimMesh<2>& mesh, const Codimension::CodimMesher<2>& mesher, const Grid<2>& grid, std::vector<double>& phi, int narrow_band){
        //initialize phi to infinity
        double init_value = grid.cell_counts.sum()*grid.dx;
        phi.resize(grid.cell_counts.prod());
        for(int i=0;i<phi.size();++i)
            phi[i] = init_value;
        //trace the closest segment of each cell
        std::vector<int> closest_seg;closest_seg.resize(grid.cell_counts.prod(), -1);
        //trace the ray intersection at cell(i,j), the ray start from (0,j) in ith direction
        std::vector<int> intersection_count;intersection_count.resize(grid.cell_counts.prod(), 0);


        const std::vector<Vec<2> >& x = mesh.Vertices();
        for(int s=0;s<mesh.elements.size();++s){
            if(!mesher.Valid_Element(s))continue;
            int p,q;
            p = mesh.elements[s](0);
            q = mesh.elements[s](1);
            Veci<2> fp = grid.Cell_Coord(x[p]);
            Veci<2> fq = grid.Cell_Coord(x[q]);
            int i0 = std::clamp(std::min(fp(0), fq(0))-narrow_band, 0, grid.cell_counts(0)-1);
            int i1 = std::clamp(std::max(fp(0), fq(0))+narrow_band, 0, grid.cell_counts(0)-1);
            int j0 = std::clamp(std::min(fp(1), fq(1))-narrow_band, 0, grid.cell_counts(1)-1);
            int j1 = std::clamp(std::max(fp(1), fq(1))+narrow_band, 0, grid.cell_counts(1)-1);

            for(int j=j0;j<=j1;++j){
                for(int i=i0;i<=i1;++i){
                    Vec<2> gx = grid.Center(Veci<2>(i,j));
                    double dist = Point_Segment_Distance<2>(gx, x[p], x[q]);
                    int idx = grid.Cell_Index(Veci<2>(i,j));
                    if(dist < phi[idx]){
                        phi[idx] = dist;
                        closest_seg[idx] = s;
                    }
                }
            }

            // intersection count
            j0 = std::clamp(std::min(fp(1),fq(1)), 0, grid.cell_counts(1)-1);
            j1 = std::clamp(std::max(fp(1),fq(1)), 0, grid.cell_counts(1)-1);
            for(int j=j0; j<=j1; ++j){
                double a;
                double pj = grid.Center(Veci<2>(0,j))(1);//get j's position
                if((pj>=x[p](1) && pj<x[q](1)) || (pj>=x[q](1) && pj<x[p](1))){
                    a = std::fabs((pj-x[q](1))/(x[p](1)-x[q](1)));
                    a = std::clamp(a, 0.0, 1.0);
                    double pi = a*x[p](0) + (1-a)*x[q](0);//get intersection point's ith coordinate
                    if(pi > grid.domain_min(0) && pi < grid.domain_max(0)){
                        Veci<2> icell = grid.Cell_Coord(Vec<2>(pi, pj));
                        if(pi > grid.Center(icell)(0)){
                            //intersection happened after the center of icell
                            icell(0) += 1;//change the sign of next cell
                        }
                        if(!grid.Valid_Cell(icell))continue;
                        int idx = grid.Cell_Index(icell);
                        ++intersection_count[idx];//ray along x coordinate from (0,j)'s center intersect with this tri at icell
                    }
                }
            }
        }

        //update i0's phi using i1
        auto check_neighbor = [&] (const Veci<2>& i0, const Veci<2>& i1) {
            if(closest_seg[grid.Cell_Index(i1)] >= 0){
                int cseg = closest_seg[grid.Cell_Index(i1)];
                int tp = mesh.elements[cseg](0);
                int tq = mesh.elements[cseg](1);
                Vec<2> gx = grid.Center(i0);
                double td = Point_Segment_Distance<2>(gx, x[tp], x[tq]);
                if(td < phi[grid.Cell_Index(i0)]){
                    phi[grid.Cell_Index(i0)] = td;
                    closest_seg[grid.Cell_Index(i0)] = closest_seg[grid.Cell_Index(i1)];
                }
            }
        };

        //sweep in (di,dj,dk) direction
        auto sweep = [&] (int di, int dj) {
            int ti0, ti1;
            if(di>0) {ti0=1;ti1=grid.cell_counts(0);}
            else {ti0=grid.cell_counts(0)-2;ti1=-1;}
            int tj0, tj1;
            if(dj>0) {tj0=1;tj1=grid.cell_counts(1);}
            else {tj0=grid.cell_counts(1)-2;tj1=-1;}
            for(int tj=tj0; tj!=tj1; tj+=dj){
                for(int ti=ti0; ti!=ti1; ti+=di){
                    check_neighbor(Veci<2>(ti,tj), Veci<2>(ti-di,tj));
                    check_neighbor(Veci<2>(ti,tj), Veci<2>(ti,tj-dj));
                    check_neighbor(Veci<2>(ti,tj), Veci<2>(ti-di,tj-dj));
                }
            }
        };

        //sweep in all direction
        for(int pass=0; pass<2; ++pass){
            sweep(+1,+1);
            sweep(-1,-1);
            sweep(+1,-1);
            sweep(-1,+1);
        }

        //determine sign from intersection
        for(int j=0; j<grid.cell_counts(1); ++j){
            int total_count=0;
            for(int i=0; i<grid.cell_counts(0); ++i){
                int idx = grid.Cell_Index(Veci<2>(i,j));
                total_count += intersection_count[idx];
                if(total_count%2 == 1){
                    phi[idx] = -phi[idx];
                }
            }
        }
    }

    template<>
    void Mesh_To_SDF<3>(const Codimension::CodimMesh<3>& mesh, const Codimension::CodimMesher<3>& mesher, const Grid<3>& grid, std::vector<double>& phi, int narrow_band){
        //initialize phi to infinity
        double init_value = grid.cell_counts.sum()*grid.dx;
        phi.resize(grid.cell_counts.prod());
        for(int i=0;i<phi.size();++i)
            phi[i] = init_value;
        //trace the closest triangle of each cell
        std::vector<int> closest_tri;closest_tri.resize(grid.cell_counts.prod(), -1);
        //trace the ray intersection at cell(i,j,k), the ray start from (0,j,k) in ith direction
        std::vector<int> intersection_count;intersection_count.resize(grid.cell_counts.prod(), 0);


        const std::vector<Vec<3> >& x = mesh.Vertices();
        for(int t=0;t<mesh.elements.size();++t){
            if(!mesher.Valid_Element(t))continue;
            int p,q,r;
            p = mesh.elements[t](0);
            q = mesh.elements[t](1);
            r = mesh.elements[t](2);
            Veci<3> fp = grid.Cell_Coord(x[p]);
            Veci<3> fq = grid.Cell_Coord(x[q]);
            Veci<3> fr = grid.Cell_Coord(x[r]);
            int i0 = std::clamp(Min(fp(0), fq(0), fr(0))-narrow_band, 0, grid.cell_counts(0)-1);
            int i1 = std::clamp(Max(fp(0), fq(0), fr(0))+narrow_band, 0, grid.cell_counts(0)-1);
            int j0 = std::clamp(Min(fp(1), fq(1), fr(1))-narrow_band, 0, grid.cell_counts(1)-1);
            int j1 = std::clamp(Max(fp(1), fq(1), fr(1))+narrow_band, 0, grid.cell_counts(1)-1);
            int k0 = std::clamp(Min(fp(2), fq(2), fr(2))-narrow_band, 0, grid.cell_counts(2)-1);
            int k1 = std::clamp(Max(fp(2), fq(2), fr(2))+narrow_band, 0, grid.cell_counts(2)-1);

            for(int k=k0;k<=k1;++k){
                for(int j=j0;j<=j1;++j){
                    for(int i=i0;i<=i1;++i){
                        Vec<3> gx = grid.Center(Veci<3>(i,j,k));
                        double dist = Point_Triangle_Distance(gx, x[p], x[q], x[r]);
                        int idx = grid.Cell_Index(Veci<3>(i,j,k));
                        if(dist < phi[idx]){
                            phi[idx] = dist;
                            closest_tri[idx] = t;
                        }
                    }
                }
            }

            // intersection count
            j0 = std::clamp(Min(fp(1),fq(1),fr(1)), 0, grid.cell_counts(1)-1);
            j1 = std::clamp(Max(fp(1),fq(1),fr(1)), 0, grid.cell_counts(1)-1);
            k0 = std::clamp(Min(fp(2),fq(2),fr(2)), 0, grid.cell_counts(2)-1);
            k1 = std::clamp(Max(fp(2),fq(2),fr(2)), 0, grid.cell_counts(2)-1);
            for(int k=k0; k<=k1; ++k){
                for(int j=j0; j<=j1; ++j){
                    double a, b, c;
                    Vec<3> ljk = grid.Center(Veci<3>(0,j,k));//get j,k's position
                    if(Point_In_Triangle_2D(ljk.segment<2>(1), x[p].segment<2>(1), x[q].segment<2>(1), x[r].segment<2>(1), a, b, c)){
                        double li = a*x[p](0) + b*x[q](0) + c*x[r](0);//get intersection point's ith coordinate
                        if(li > grid.domain_min(0) && li < grid.domain_max(0)){
                            Veci<3> icell = grid.Cell_Coord(Vec<3>(li,ljk(1),ljk(2)));
                            int idx = grid.Cell_Index(icell);
                            if(li > grid.Center(icell)(0)){
                                //intersection happened after the center of icell
                                icell(0) += 1;//change the sign of next cell
                            }
                            if(!grid.Valid_Cell(icell))continue;
                            ++intersection_count[idx];//ray along x coordinate from (0,j,k)'s center intersect with this tri at icell
                        }
                    }
                }
            }
        }

        //update i0's phi using i1
        auto check_neighbor = [&] (const Veci<3>& i0, const Veci<3>& i1) {
            if(closest_tri[grid.Cell_Index(i1)] >= 0){
                int ctri = closest_tri[grid.Cell_Index(i1)];
                int tp = mesh.elements[ctri](0);
                int tq = mesh.elements[ctri](1);
                int tr = mesh.elements[ctri](2);
                Vec<3> gx = grid.Center(i0);
                double td = Point_Triangle_Distance(gx, x[tp], x[tq], x[tr]);
                if(td < phi[grid.Cell_Index(i0)]){
                    phi[grid.Cell_Index(i0)] = td;
                    closest_tri[grid.Cell_Index(i0)] = closest_tri[grid.Cell_Index(i1)];
                }
            }
        };

        //sweep in (di,dj,dk) direction
        auto sweep = [&] (int di, int dj, int dk) {
            int ti0, ti1;
            if(di>0) {ti0=1;ti1=grid.cell_counts(0);}
            else {ti0=grid.cell_counts(0)-2;ti1=-1;}
            int tj0, tj1;
            if(dj>0) {tj0=1;tj1=grid.cell_counts(1);}
            else {tj0=grid.cell_counts(1)-2;tj1=-1;}
            int tk0, tk1;
            if(dk>0) {tk0=1;tk1=grid.cell_counts(2);}
            else {tk0=grid.cell_counts(2)-2;tk1=-1;}
            for(int tk=tk0; tk!=tk1; tk+=dk){
                for(int tj=tj0; tj!=tj1; tj+=dj){
                    for(int ti=ti0; ti!=ti1; ti+=di){
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti-di,tj,tk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti,tj-dj,tk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti-di,tj-dj,tk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti,tj,tk-dk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti-di,tj,tk-dk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti,tj-dj,tk-dk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti-di,tj-dj,tk-dk));
                    }
                }
            }
        };

        //sweep in all direction
        for(int pass=0; pass<2; ++pass){
            sweep(+1,+1,+1);
            sweep(-1,-1,-1);
            sweep(+1,+1,-1);
            sweep(-1,-1,+1);
            sweep(+1,-1,+1);
            sweep(-1,+1,-1);
            sweep(+1,-1,-1);
            sweep(-1,+1,+1);
        }

        //determine sign from intersection
        for(int k=0;k<grid.cell_counts(2);++k){
            for(int j=0;j<grid.cell_counts(1);++j){
                int total_count=0;
                for(int i=0; i<grid.cell_counts(0); ++i){
                    int idx = grid.Cell_Index(Veci<3>(i,j,k));
                    total_count += intersection_count[idx];
                    if(total_count%2 == 1){
                        phi[idx] = -phi[idx];
                    }
                }
            }
        }
    }

    template<>
    void Mesh_To_SDF_Narrow_Band<2>(const Codimension::CodimMesh<2>& mesh, const Codimension::CodimMesher<2>& mesher, const Grid<2>& grid, std::vector<double>& phi, int narrow_band){
        //initialize phi to infinity
        double init_value = grid.cell_counts.sum()*grid.dx;
        phi.resize(grid.cell_counts.prod());
        for(int i=0;i<phi.size();++i)
            phi[i] = init_value;
        //trace the closest segment of each cell
        std::vector<int> closest_seg;closest_seg.resize(grid.cell_counts.prod(), -1);
        //trace the ray intersection at cell(i,j), the ray start from (0,j) in ith direction
        std::vector<int> intersection_count;intersection_count.resize(grid.cell_counts.prod(), 0);


        const std::vector<Vec<2> >& x = mesh.Vertices();
        for(int s=0;s<mesh.elements.size();++s){
            if(!mesher.Valid_Element(s))continue;
            int p,q;
            p = mesh.elements[s](0);
            q = mesh.elements[s](1);
            Veci<2> fp = grid.Cell_Coord(x[p]);
            Veci<2> fq = grid.Cell_Coord(x[q]);
            int i0 = std::clamp(std::min(fp(0), fq(0))-narrow_band, 0, grid.cell_counts(0)-1);
            int i1 = std::clamp(std::max(fp(0), fq(0))+narrow_band, 0, grid.cell_counts(0)-1);
            int j0 = std::clamp(std::min(fp(1), fq(1))-narrow_band, 0, grid.cell_counts(1)-1);
            int j1 = std::clamp(std::max(fp(1), fq(1))+narrow_band, 0, grid.cell_counts(1)-1);

            for(int j=j0;j<=j1;++j){
                for(int i=i0;i<=i1;++i){
                    Vec<2> gx = grid.Center(Veci<2>(i,j));
                    double dist = Point_Segment_Distance<2>(gx, x[p], x[q]);
                    int idx = grid.Cell_Index(Veci<2>(i,j));
                    if(dist < phi[idx]){
                        phi[idx] = dist;
                        closest_seg[idx] = s;
                    }
                }
            }

            // intersection count
            j0 = std::clamp(std::min(fp(1),fq(1)), 0, grid.cell_counts(1)-1);
            j1 = std::clamp(std::max(fp(1),fq(1)), 0, grid.cell_counts(1)-1);
            for(int j=j0; j<=j1; ++j){
                double a;
                double pj = grid.Center(Veci<2>(0,j))(1);//get j's position
                if((pj>=x[p](1) && pj<x[q](1)) || (pj>=x[q](1) && pj<x[p](1))){
                    a = std::fabs((pj-x[q](1))/(x[p](1)-x[q](1)));
                    a = std::clamp(a, 0.0, 1.0);
                    double pi = a*x[p](0) + (1-a)*x[q](0);//get intersection point's ith coordinate
                    if(pi > grid.domain_min(0) && pi < grid.domain_max(0)){
                        Veci<2> icell = grid.Cell_Coord(Vec<2>(pi, pj));
                        if(pi > grid.Center(icell)(0)){
                            //intersection happened after the center of icell
                            icell(0) += 1;//change the sign of next cell
                        }
                        if(!grid.Valid_Cell(icell))continue;
                        int idx = grid.Cell_Index(icell);
                        ++intersection_count[idx];//ray along x coordinate from (0,j)'s center intersect with this tri at icell
                    }
                }
            }
        }

        //determine sign from intersection
        for(int j=0; j<grid.cell_counts(1); ++j){
            int total_count=0;
            for(int i=0; i<grid.cell_counts(0); ++i){
                int idx = grid.Cell_Index(Veci<2>(i,j));
                total_count += intersection_count[idx];
                if(total_count%2 == 1){
                    phi[idx] = -phi[idx];
                }
            }
        }
    }

    template<>
    void Mesh_To_SDF_Narrow_Band<3>(const Codimension::CodimMesh<3>& mesh, const Codimension::CodimMesher<3>& mesher, const Grid<3>& grid, std::vector<double>& phi, int narrow_band){
        //initialize phi to infinity
        double init_value = grid.cell_counts.sum()*grid.dx;
        phi.resize(grid.cell_counts.prod());
        for(int i=0;i<phi.size();++i)
            phi[i] = init_value;
        //trace the closest triangle of each cell
        std::vector<int> closest_tri;closest_tri.resize(grid.cell_counts.prod(), -1);
        //trace the ray intersection at cell(i,j,k), the ray start from (0,j,k) in ith direction
        std::vector<int> intersection_count;intersection_count.resize(grid.cell_counts.prod(), 0);


        const std::vector<Vec<3> >& x = mesh.Vertices();
        for(int t=0;t<mesh.elements.size();++t){
            if(!mesher.Valid_Element(t))continue;
            int p,q,r;
            p = mesh.elements[t](0);
            q = mesh.elements[t](1);
            r = mesh.elements[t](2);
            Veci<3> fp = grid.Cell_Coord(x[p]);
            Veci<3> fq = grid.Cell_Coord(x[q]);
            Veci<3> fr = grid.Cell_Coord(x[r]);
            int i0 = std::clamp(Min(fp(0), fq(0), fr(0))-narrow_band, 0, grid.cell_counts(0)-1);
            int i1 = std::clamp(Max(fp(0), fq(0), fr(0))+narrow_band, 0, grid.cell_counts(0)-1);
            int j0 = std::clamp(Min(fp(1), fq(1), fr(1))-narrow_band, 0, grid.cell_counts(1)-1);
            int j1 = std::clamp(Max(fp(1), fq(1), fr(1))+narrow_band, 0, grid.cell_counts(1)-1);
            int k0 = std::clamp(Min(fp(2), fq(2), fr(2))-narrow_band, 0, grid.cell_counts(2)-1);
            int k1 = std::clamp(Max(fp(2), fq(2), fr(2))+narrow_band, 0, grid.cell_counts(2)-1);

            for(int k=k0;k<=k1;++k){
                for(int j=j0;j<=j1;++j){
                    for(int i=i0;i<=i1;++i){
                        Vec<3> gx = grid.Center(Veci<3>(i,j,k));
                        double dist = Point_Triangle_Distance(gx, x[p], x[q], x[r]);
                        int idx = grid.Cell_Index(Veci<3>(i,j,k));
                        if(dist < phi[idx]){
                            phi[idx] = dist;
                            closest_tri[idx] = t;
                        }
                    }
                }
            }

            // intersection count
            j0 = std::clamp(Min(fp(1),fq(1),fr(1)), 0, grid.cell_counts(1)-1);
            j1 = std::clamp(Max(fp(1),fq(1),fr(1)), 0, grid.cell_counts(1)-1);
            k0 = std::clamp(Min(fp(2),fq(2),fr(2)), 0, grid.cell_counts(2)-1);
            k1 = std::clamp(Max(fp(2),fq(2),fr(2)), 0, grid.cell_counts(2)-1);
            for(int k=k0; k<=k1; ++k){
                for(int j=j0; j<=j1; ++j){
                    double a, b, c;
                    Vec<3> ljk = grid.Center(Veci<3>(0,j,k));//get j,k's position
                    if(Point_In_Triangle_2D(ljk.segment<2>(1), x[p].segment<2>(1), x[q].segment<2>(1), x[r].segment<2>(1), a, b, c)){
                        double li = a*x[p](0) + b*x[q](0) + c*x[r](0);//get intersection point's ith coordinate
                        if(li > grid.domain_min(0) && li < grid.domain_max(0)){
                            Veci<3> icell = grid.Cell_Coord(Vec<3>(li,ljk(1),ljk(2)));
                            int idx = grid.Cell_Index(icell);
                            if(li > grid.Center(icell)(0)){
                                //intersection happened after the center of icell
                                icell(0) += 1;//change the sign of next cell
                            }
                            if(!grid.Valid_Cell(icell))continue;
                            ++intersection_count[idx];//ray along x coordinate from (0,j,k)'s center intersect with this tri at icell
                        }
                    }
                }
            }
        }

        //determine sign from intersection
        for(int k=0;k<grid.cell_counts(2);++k){
            for(int j=0;j<grid.cell_counts(1);++j){
                int total_count=0;
                for(int i=0; i<grid.cell_counts(0); ++i){
                    int idx = grid.Cell_Index(Veci<3>(i,j,k));
                    total_count += intersection_count[idx];
                    if(total_count%2 == 1){
                        phi[idx] = -phi[idx];
                    }
                }
            }
        }
    }

    template<>
    void Mesh_To_SDF_Directional<2>(const SurfaceMesh<2>& mesh, const Grid<2>& grid, std::vector<double>& phi, int axis, int narrow_band){
        if(axis==2){std::cerr << "no z axis in 2D" << std::endl;}
        //initialize phi to infinity
        double init_value = grid.cell_counts.sum()*grid.dx;
        phi.resize(grid.cell_counts.prod());
        for(int i=0;i<phi.size();++i)
            phi[i] = init_value;
        //trace the closest segment of each cell
        std::vector<int> closest_seg;closest_seg.resize(grid.cell_counts.prod(), -1);
        //trace the ray intersection at cell(i,j), the ray start from (0,j) in ith direction
        std::vector<int> intersection_count;intersection_count.resize(grid.cell_counts.prod(), 0);


        const std::vector<Vec<2> >& x = mesh.Vertices();
        for(int s=0;s<mesh.elements.size();++s){
            int p,q;
            p = mesh.elements[s](0);
            q = mesh.elements[s](1);
            Veci<2> fp = grid.Cell_Coord(x[p]);
            Veci<2> fq = grid.Cell_Coord(x[q]);
            int i0 = std::clamp(std::min(fp(0), fq(0))-narrow_band, 0, grid.cell_counts(0)-1);
            int i1 = std::clamp(std::max(fp(0), fq(0))+narrow_band, 0, grid.cell_counts(0)-1);
            int j0 = std::clamp(std::min(fp(1), fq(1))-narrow_band, 0, grid.cell_counts(1)-1);
            int j1 = std::clamp(std::max(fp(1), fq(1))+narrow_band, 0, grid.cell_counts(1)-1);

            for(int j=j0;j<=j1;++j){
                for(int i=i0;i<=i1;++i){
                    Vec<2> gx = grid.Center(Veci<2>(i,j));
                    double dist = Point_Segment_Distance<2>(gx, x[p], x[q]);
                    int idx = grid.Cell_Index(Veci<2>(i,j));
                    if(dist < phi[idx]){
                        phi[idx] = dist;
                        closest_seg[idx] = s;
                    }
                }
            }

            // intersection count, ray send in axis direction
            int tj0 = std::clamp(std::min(fp(1-axis),fq(1-axis)), 0, grid.cell_counts(1-axis)-1);
            int tj1 = std::clamp(std::max(fp(1-axis),fq(1-axis)), 0, grid.cell_counts(1-axis)-1);
            for(int tj=tj0; tj<=tj1; ++tj){
                double a;
                Veci<2> ftj = Veci<2>::Zero();
                ftj(1-axis) = tj;
                double ptj = grid.Center(ftj)(1-axis);//get tj's position
                if((ptj>=x[p](1-axis) && ptj<x[q](1-axis)) || (ptj>=x[q](1-axis) && ptj<x[p](1-axis))){
                    a = std::fabs((ptj-x[q](1-axis))/(x[p](1-axis)-x[q](1-axis)));
                    a = std::clamp(a, 0.0, 1.0);
                    double pti = a*x[p](axis) + (1-a)*x[q](axis);//get intersection point's ith coordinate
                    if(pti > grid.domain_min(axis) && pti < grid.domain_max(axis)){
                        Vec<2> gt = Vec<2>::Zero(); gt(axis) = pti; gt(1-axis)=ptj;
                        Veci<2> icell = grid.Cell_Coord(gt);
                        if(pti > grid.Center(icell)(axis)){
                            //intersection happened after the center of icell
                            icell(axis) += 1;//change the sign of next cell
                        }
                        if(!grid.Valid_Cell(icell))continue;
                        int idx = grid.Cell_Index(icell);
                        ++intersection_count[idx];//ray along x coordinate from (0,j)'s center intersect with this tri at icell
                    }
                }
            }
        }

        //update i0's phi using i1
        auto check_neighbor = [&] (const Veci<2>& i0, const Veci<2>& i1) {
            if(closest_seg[grid.Cell_Index(i1)] >= 0){
                int cseg = closest_seg[grid.Cell_Index(i1)];
                int tp = mesh.elements[cseg](0);
                int tq = mesh.elements[cseg](1);
                Vec<2> gx = grid.Center(i0);
                double td = Point_Segment_Distance<2>(gx, x[tp], x[tq]);
                if(td < phi[grid.Cell_Index(i0)]){
                    phi[grid.Cell_Index(i0)] = td;
                    closest_seg[grid.Cell_Index(i0)] = closest_seg[grid.Cell_Index(i1)];
                }
            }
        };

        //sweep in (di,dj,dk) direction
        auto sweep = [&] (int di, int dj) {
            int ti0, ti1;
            if(di>0) {ti0=1;ti1=grid.cell_counts(0);}
            else {ti0=grid.cell_counts(0)-2;ti1=-1;}
            int tj0, tj1;
            if(dj>0) {tj0=1;tj1=grid.cell_counts(1);}
            else {tj0=grid.cell_counts(1)-2;tj1=-1;}
            for(int tj=tj0; tj!=tj1; tj+=dj){
                for(int ti=ti0; ti!=ti1; ti+=di){
                    check_neighbor(Veci<2>(ti,tj), Veci<2>(ti-di,tj));
                    check_neighbor(Veci<2>(ti,tj), Veci<2>(ti,tj-dj));
                    check_neighbor(Veci<2>(ti,tj), Veci<2>(ti-di,tj-dj));
                }
            }
        };

        //sweep in all direction
        for(int pass=0; pass<2; ++pass){
            sweep(+1,+1);
            sweep(-1,-1);
            sweep(+1,-1);
            sweep(-1,+1);
        }

        //determine sign from intersection
        for(int j=0; j<grid.cell_counts(1-axis); ++j){
            int total_count=0;
            for(int i=0; i<grid.cell_counts(axis); ++i){
                Veci<2> cell = Veci<2>::Zero(); cell(axis)=i; cell(1-axis)=j;
                int idx = grid.Cell_Index(cell);
                total_count += intersection_count[idx];
                //from bottom to top, bottom's sign is negative
                if(total_count%2 == 0){
                    phi[idx] = -phi[idx];
                }
            }
        }
    }

    template<>
    void Mesh_To_SDF_Directional<3>(const SurfaceMesh<3>& mesh, const Grid<3>& grid, std::vector<double>& phi, int axis, int narrow_band){
        //initialize phi to infinity
        double init_value = grid.cell_counts.sum()*grid.dx;
        phi.resize(grid.cell_counts.prod());
        for(int i=0;i<phi.size();++i)
            phi[i] = init_value;
        //trace the closest triangle of each cell
        std::vector<int> closest_tri;closest_tri.resize(grid.cell_counts.prod(), -1);
        //trace the ray intersection at cell(i,j,k), the ray start from (0,j,k) in ith direction
        std::vector<int> intersection_count;intersection_count.resize(grid.cell_counts.prod(), 0);

        int ai = axis;
        int aj = (axis+1)%3;
        int ak = (axis+2)%3;
        const std::vector<Vec<3> >& x = mesh.Vertices();
        for(int t=0;t<mesh.elements.size();++t){
            int p,q,r;
            p = mesh.elements[t](0);
            q = mesh.elements[t](1);
            r = mesh.elements[t](2);
            Veci<3> fp = grid.Cell_Coord(x[p]);
            Veci<3> fq = grid.Cell_Coord(x[q]);
            Veci<3> fr = grid.Cell_Coord(x[r]);
            int i0 = std::clamp(Min(fp(0), fq(0), fr(0))-narrow_band, 0, grid.cell_counts(0)-1);
            int i1 = std::clamp(Max(fp(0), fq(0), fr(0))+narrow_band, 0, grid.cell_counts(0)-1);
            int j0 = std::clamp(Min(fp(1), fq(1), fr(1))-narrow_band, 0, grid.cell_counts(1)-1);
            int j1 = std::clamp(Max(fp(1), fq(1), fr(1))+narrow_band, 0, grid.cell_counts(1)-1);
            int k0 = std::clamp(Min(fp(2), fq(2), fr(2))-narrow_band, 0, grid.cell_counts(2)-1);
            int k1 = std::clamp(Max(fp(2), fq(2), fr(2))+narrow_band, 0, grid.cell_counts(2)-1);

            for(int k=k0;k<=k1;++k){
                for(int j=j0;j<=j1;++j){
                    for(int i=i0;i<=i1;++i){
                        Vec<3> gx = grid.Center(Veci<3>(i,j,k));
                        double dist = Point_Triangle_Distance(gx, x[p], x[q], x[r]);
                        int idx = grid.Cell_Index(Veci<3>(i,j,k));
                        if(dist < phi[idx]){
                            phi[idx] = dist;
                            closest_tri[idx] = t;
                        }
                    }
                }
            }

            // intersection count
            j0 = std::clamp(Min(fp(aj),fq(aj),fr(aj)), 0, grid.cell_counts(aj)-1);
            j1 = std::clamp(Max(fp(aj),fq(aj),fr(aj)), 0, grid.cell_counts(aj)-1);
            k0 = std::clamp(Min(fp(ak),fq(ak),fr(ak)), 0, grid.cell_counts(ak)-1);
            k1 = std::clamp(Max(fp(ak),fq(ak),fr(ak)), 0, grid.cell_counts(ak)-1);
            for(int k=k0; k<=k1; ++k){
                for(int j=j0; j<=j1; ++j){
                    double a, b, c;
                    Veci<3> fjk = Veci<3>::Zero(); fjk(aj) = j; fjk(ak) = k;
                    Vec<3> ljk = grid.Center(fjk);//get j,k's position
                    if(Point_In_Triangle_2D(Vec<2>(ljk(aj), ljk(ak)),
                                            Vec<2>(x[p](aj),x[p](ak)), Vec<2>(x[q](aj),x[q](ak)), Vec<2>(x[r](aj),x[r](ak)), a, b, c)){
                        double li = a*x[p](ai) + b*x[q](ai) + c*x[r](ai);//get intersection point's ith coordinate
                        Vec<3> gt = Vec<3>::Zero(); gt(ai) = li; gt(aj) = ljk(aj); gt(ak) = ljk(ak);
                        Veci<3> icell = grid.Cell_Coord(gt);
                        if(li > grid.Center(icell)(ai)){
                            //intersection happened after the center of icell
                            icell(ai) += 1;//change the sign of next cell
                        }
                        if(!grid.Valid_Cell(icell))continue;
                        int idx = grid.Cell_Index(icell);
                        ++intersection_count[idx];//ray along x coordinate from (0,j,k)'s center intersect with this tri at icell
                    }
                }
            }
        }

        //update i0's phi using i1
        auto check_neighbor = [&] (const Veci<3>& i0, const Veci<3>& i1) {
            if(closest_tri[grid.Cell_Index(i1)] >= 0){
                int ctri = closest_tri[grid.Cell_Index(i1)];
                int tp = mesh.elements[ctri](0);
                int tq = mesh.elements[ctri](1);
                int tr = mesh.elements[ctri](2);
                Vec<3> gx = grid.Center(i0);
                double td = Point_Triangle_Distance(gx, x[tp], x[tq], x[tr]);
                if(td < phi[grid.Cell_Index(i0)]){
                    phi[grid.Cell_Index(i0)] = td;
                    closest_tri[grid.Cell_Index(i0)] = closest_tri[grid.Cell_Index(i1)];
                }
            }
        };

        //sweep in (di,dj,dk) direction
        auto sweep = [&] (int di, int dj, int dk) {
            int ti0, ti1;
            if(di>0) {ti0=1;ti1=grid.cell_counts(0);}
            else {ti0=grid.cell_counts(0)-2;ti1=-1;}
            int tj0, tj1;
            if(dj>0) {tj0=1;tj1=grid.cell_counts(1);}
            else {tj0=grid.cell_counts(1)-2;tj1=-1;}
            int tk0, tk1;
            if(dk>0) {tk0=1;tk1=grid.cell_counts(2);}
            else {tk0=grid.cell_counts(2)-2;tk1=-1;}
            for(int tk=tk0; tk!=tk1; tk+=dk){
                for(int tj=tj0; tj!=tj1; tj+=dj){
                    for(int ti=ti0; ti!=ti1; ti+=di){
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti-di,tj,tk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti,tj-dj,tk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti-di,tj-dj,tk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti,tj,tk-dk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti-di,tj,tk-dk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti,tj-dj,tk-dk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti-di,tj-dj,tk-dk));
                    }
                }
            }
        };

        //sweep in all direction
        for(int pass=0; pass<2; ++pass){
            sweep(+1,+1,+1);
            sweep(-1,-1,-1);
            sweep(+1,+1,-1);
            sweep(-1,-1,+1);
            sweep(+1,-1,+1);
            sweep(-1,+1,-1);
            sweep(+1,-1,-1);
            sweep(-1,+1,+1);
        }

        //determine sign from intersection
        for(int k=0;k<grid.cell_counts(ak);++k){
            for(int j=0;j<grid.cell_counts(aj);++j){
                int total_count=0;
                for(int i=0; i<grid.cell_counts(ai); ++i){
                    Veci<3> cell = Veci<3>::Zero(); cell(ai) = i; cell(aj) = j; cell(ak) = k;
                    int idx = grid.Cell_Index(cell);
                    total_count += intersection_count[idx];
                    if(total_count%2 == 0){//(0,0,0)'s sign is negative
                        phi[idx] = -phi[idx];
                    }
                }
            }
        }
    }

    template<>
    void Mesh_To_SDF_Directional<2>(const Codimension::CodimMesh<2>& mesh, const Codimension::CodimMesher<2>& mesher, const Grid<2>& grid, std::vector<double>& phi, int axis, int narrow_band){
        if(axis==2){std::cerr << "no z axis in 2D" << std::endl;}
        //initialize phi to infinity
        double init_value = grid.cell_counts.sum()*grid.dx;
        phi.resize(grid.cell_counts.prod());
        for(int i=0;i<phi.size();++i)
            phi[i] = init_value;
        //trace the closest segment of each cell
        std::vector<int> closest_seg;closest_seg.resize(grid.cell_counts.prod(), -1);
        //trace the ray intersection at cell(i,j), the ray start from (0,j) in ith direction
        std::vector<int> intersection_count;intersection_count.resize(grid.cell_counts.prod(), 0);


        const std::vector<Vec<2> >& x = mesh.Vertices();
        for(int s=0;s<mesh.elements.size();++s){
            if(!mesher.Valid_Element(s))continue;
            int p,q;
            p = mesh.elements[s](0);
            q = mesh.elements[s](1);
            Veci<2> fp = grid.Cell_Coord(x[p]);
            Veci<2> fq = grid.Cell_Coord(x[q]);
            int i0 = std::clamp(std::min(fp(0), fq(0))-narrow_band, 0, grid.cell_counts(0)-1);
            int i1 = std::clamp(std::max(fp(0), fq(0))+narrow_band, 0, grid.cell_counts(0)-1);
            int j0 = std::clamp(std::min(fp(1), fq(1))-narrow_band, 0, grid.cell_counts(1)-1);
            int j1 = std::clamp(std::max(fp(1), fq(1))+narrow_band, 0, grid.cell_counts(1)-1);

            for(int j=j0;j<=j1;++j){
                for(int i=i0;i<=i1;++i){
                    Vec<2> gx = grid.Center(Veci<2>(i,j));
                    double dist = Point_Segment_Distance<2>(gx, x[p], x[q]);
                    int idx = grid.Cell_Index(Veci<2>(i,j));
                    if(dist < phi[idx]){
                        phi[idx] = dist;
                        closest_seg[idx] = s;
                    }
                }
            }

            // intersection count, ray send in axis direction
            int tj0 = std::clamp(std::min(fp(1-axis),fq(1-axis)), 0, grid.cell_counts(1-axis)-1);
            int tj1 = std::clamp(std::max(fp(1-axis),fq(1-axis)), 0, grid.cell_counts(1-axis)-1);
            for(int tj=tj0; tj<=tj1; ++tj){
                double a;
                Veci<2> ftj = Veci<2>::Zero();
                ftj(1-axis) = tj;
                double ptj = grid.Center(ftj)(1-axis);//get tj's position
                if((ptj>=x[p](1-axis) && ptj<x[q](1-axis)) || (ptj>=x[q](1-axis) && ptj<x[p](1-axis))){
                    a = std::fabs((ptj-x[q](1-axis))/(x[p](1-axis)-x[q](1-axis)));
                    a = std::clamp(a, 0.0, 1.0);
                    double pti = a*x[p](axis) + (1-a)*x[q](axis);//get intersection point's ith coordinate
                    if(pti > grid.domain_min(axis) && pti < grid.domain_max(axis)){
                        Vec<2> gt = Vec<2>::Zero(); gt(axis) = pti; gt(1-axis)=ptj;
                        Veci<2> icell = grid.Cell_Coord(gt);
                        if(pti > grid.Center(icell)(axis)){
                            //intersection happened after the center of icell
                            icell(axis) += 1;//change the sign of next cell
                        }
                        if(!grid.Valid_Cell(icell))continue;
                        int idx = grid.Cell_Index(icell);
                        ++intersection_count[idx];
                    }
                }
            }
        }

        //update i0's phi using i1
        auto check_neighbor = [&] (const Veci<2>& i0, const Veci<2>& i1) {
            if(closest_seg[grid.Cell_Index(i1)] >= 0){
                int cseg = closest_seg[grid.Cell_Index(i1)];
                int tp = mesh.elements[cseg](0);
                int tq = mesh.elements[cseg](1);
                Vec<2> gx = grid.Center(i0);
                double td = Point_Segment_Distance<2>(gx, x[tp], x[tq]);
                if(td < phi[grid.Cell_Index(i0)]){
                    phi[grid.Cell_Index(i0)] = td;
                    closest_seg[grid.Cell_Index(i0)] = closest_seg[grid.Cell_Index(i1)];
                }
            }
        };

        //sweep in (di,dj,dk) direction
        auto sweep = [&] (int di, int dj) {
            int ti0, ti1;
            if(di>0) {ti0=1;ti1=grid.cell_counts(0);}
            else {ti0=grid.cell_counts(0)-2;ti1=-1;}
            int tj0, tj1;
            if(dj>0) {tj0=1;tj1=grid.cell_counts(1);}
            else {tj0=grid.cell_counts(1)-2;tj1=-1;}
            for(int tj=tj0; tj!=tj1; tj+=dj){
                for(int ti=ti0; ti!=ti1; ti+=di){
                    check_neighbor(Veci<2>(ti,tj), Veci<2>(ti-di,tj));
                    check_neighbor(Veci<2>(ti,tj), Veci<2>(ti,tj-dj));
                    check_neighbor(Veci<2>(ti,tj), Veci<2>(ti-di,tj-dj));
                }
            }
        };

        //sweep in all direction
        for(int pass=0; pass<2; ++pass){
            sweep(+1,+1);
            sweep(-1,-1);
            sweep(+1,-1);
            sweep(-1,+1);
        }

        //determine sign from intersection
        for(int j=0; j<grid.cell_counts(1-axis); ++j){
            int total_count=0;
            for(int i=0; i<grid.cell_counts(axis); ++i){
                Veci<2> cell = Veci<2>::Zero(); cell(axis)=i; cell(1-axis)=j;
                int idx = grid.Cell_Index(cell);
                total_count += intersection_count[idx];
                //from bottom to top, bottom's sign is negative
                if(total_count%2 == 0){
                    phi[idx] = -phi[idx];
                }
            }
        }
    }

    template<>
    void Mesh_To_SDF_Directional<3>(const Codimension::CodimMesh<3>& mesh, const Codimension::CodimMesher<3>& mesher, const Grid<3>& grid, std::vector<double>& phi, int axis, int narrow_band){
        //initialize phi to infinity
        double init_value = grid.cell_counts.sum()*grid.dx;
        phi.resize(grid.cell_counts.prod());
        for(int i=0;i<phi.size();++i)
            phi[i] = init_value;
        //trace the closest triangle of each cell
        std::vector<int> closest_tri;closest_tri.resize(grid.cell_counts.prod(), -1);
        //trace the ray intersection at cell(i,j,k), the ray start from (0,j,k) in ith direction
        std::vector<int> intersection_count;intersection_count.resize(grid.cell_counts.prod(), 0);

        int ai = axis;
        int aj = (axis+1)%3;
        int ak = (axis+2)%3;
        const std::vector<Vec<3> >& x = mesh.Vertices();
        for(int t=0;t<mesh.elements.size();++t){
            if(!mesher.Valid_Element(t))continue;
            int p,q,r;
            p = mesh.elements[t](0);
            q = mesh.elements[t](1);
            r = mesh.elements[t](2);
            Veci<3> fp = grid.Cell_Coord(x[p]);
            Veci<3> fq = grid.Cell_Coord(x[q]);
            Veci<3> fr = grid.Cell_Coord(x[r]);
            int i0 = std::clamp(Min(fp(0), fq(0), fr(0))-narrow_band, 0, grid.cell_counts(0)-1);
            int i1 = std::clamp(Max(fp(0), fq(0), fr(0))+narrow_band, 0, grid.cell_counts(0)-1);
            int j0 = std::clamp(Min(fp(1), fq(1), fr(1))-narrow_band, 0, grid.cell_counts(1)-1);
            int j1 = std::clamp(Max(fp(1), fq(1), fr(1))+narrow_band, 0, grid.cell_counts(1)-1);
            int k0 = std::clamp(Min(fp(2), fq(2), fr(2))-narrow_band, 0, grid.cell_counts(2)-1);
            int k1 = std::clamp(Max(fp(2), fq(2), fr(2))+narrow_band, 0, grid.cell_counts(2)-1);

            for(int k=k0;k<=k1;++k){
                for(int j=j0;j<=j1;++j){
                    for(int i=i0;i<=i1;++i){
                        Vec<3> gx = grid.Center(Veci<3>(i,j,k));
                        double dist = Point_Triangle_Distance(gx, x[p], x[q], x[r]);
                        int idx = grid.Cell_Index(Veci<3>(i,j,k));
                        if(dist < phi[idx]){
                            phi[idx] = dist;
                            closest_tri[idx] = t;
                        }
                    }
                }
            }

            // intersection count
            j0 = std::clamp(Min(fp(aj),fq(aj),fr(aj)), 0, grid.cell_counts(aj)-1);
            j1 = std::clamp(Max(fp(aj),fq(aj),fr(aj)), 0, grid.cell_counts(aj)-1);
            k0 = std::clamp(Min(fp(ak),fq(ak),fr(ak)), 0, grid.cell_counts(ak)-1);
            k1 = std::clamp(Max(fp(ak),fq(ak),fr(ak)), 0, grid.cell_counts(ak)-1);
            for(int k=k0; k<=k1; ++k){
                for(int j=j0; j<=j1; ++j){
                    double a, b, c;
                    Veci<3> fjk = Veci<3>::Zero(); fjk(aj) = j; fjk(ak) = k;
                    Vec<3> ljk = grid.Center(fjk);//get j,k's position
                    if(Point_In_Triangle_2D(Vec<2>(ljk(aj), ljk(ak)),
                                            Vec<2>(x[p](aj),x[p](ak)), Vec<2>(x[q](aj),x[q](ak)), Vec<2>(x[r](aj),x[r](ak)), a, b, c)){
                        double li = a*x[p](ai) + b*x[q](ai) + c*x[r](ai);//get intersection point's ith coordinate
                        Vec<3> gt = Vec<3>::Zero(); gt(ai) = li; gt(aj) = ljk(aj); gt(ak) = ljk(ak);
                        Veci<3> icell = grid.Cell_Coord(gt);
                        if(li > grid.Center(icell)(ai)){
                            //intersection happened after the center of icell
                            icell(ai) += 1;//change the sign of next cell
                        }
                        if(!grid.Valid_Cell(icell))continue;
                        int idx = grid.Cell_Index(icell);
                        ++intersection_count[idx];//ray along x coordinate from (0,j,k)'s center intersect with this tri at icell
                    }
                }
            }
        }

        //update i0's phi using i1
        auto check_neighbor = [&] (const Veci<3>& i0, const Veci<3>& i1) {
            if(closest_tri[grid.Cell_Index(i1)] >= 0){
                int ctri = closest_tri[grid.Cell_Index(i1)];
                int tp = mesh.elements[ctri](0);
                int tq = mesh.elements[ctri](1);
                int tr = mesh.elements[ctri](2);
                Vec<3> gx = grid.Center(i0);
                double td = Point_Triangle_Distance(gx, x[tp], x[tq], x[tr]);
                if(td < phi[grid.Cell_Index(i0)]){
                    phi[grid.Cell_Index(i0)] = td;
                    closest_tri[grid.Cell_Index(i0)] = closest_tri[grid.Cell_Index(i1)];
                }
            }
        };

        //sweep in (di,dj,dk) direction
        auto sweep = [&] (int di, int dj, int dk) {
            int ti0, ti1;
            if(di>0) {ti0=1;ti1=grid.cell_counts(0);}
            else {ti0=grid.cell_counts(0)-2;ti1=-1;}
            int tj0, tj1;
            if(dj>0) {tj0=1;tj1=grid.cell_counts(1);}
            else {tj0=grid.cell_counts(1)-2;tj1=-1;}
            int tk0, tk1;
            if(dk>0) {tk0=1;tk1=grid.cell_counts(2);}
            else {tk0=grid.cell_counts(2)-2;tk1=-1;}
            for(int tk=tk0; tk!=tk1; tk+=dk){
                for(int tj=tj0; tj!=tj1; tj+=dj){
                    for(int ti=ti0; ti!=ti1; ti+=di){
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti-di,tj,tk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti,tj-dj,tk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti-di,tj-dj,tk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti,tj,tk-dk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti-di,tj,tk-dk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti,tj-dj,tk-dk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti-di,tj-dj,tk-dk));
                    }
                }
            }
        };

        //sweep in all direction
        for(int pass=0; pass<2; ++pass){
            sweep(+1,+1,+1);
            sweep(-1,-1,-1);
            sweep(+1,+1,-1);
            sweep(-1,-1,+1);
            sweep(+1,-1,+1);
            sweep(-1,+1,-1);
            sweep(+1,-1,-1);
            sweep(-1,+1,+1);
        }

        //determine sign from intersection
        for(int k=0;k<grid.cell_counts(ak);++k){
            for(int j=0;j<grid.cell_counts(aj);++j){
                int total_count=0;
                for(int i=0; i<grid.cell_counts(ai); ++i){
                    Veci<3> cell = Veci<3>::Zero(); cell(ai) = i; cell(aj) = j; cell(ak) = k;
                    int idx = grid.Cell_Index(cell);
                    total_count += intersection_count[idx];
                    if(total_count%2 == 0){//(0,0,0)'s sign is negative
                        phi[idx] = -phi[idx];
                    }
                }
            }
        }
    }

    template<>
    void Mesh_To_SDF_Directional_Narrow_Band<2>(const Codimension::CodimMesh<2>& mesh, const Codimension::CodimMesher<2>& mesher, const Grid<2>& grid, std::vector<double>& phi, int axis, int narrow_band){
        if(axis==2){std::cerr << "no z axis in 2D" << std::endl;}
        //initialize phi to infinity
        double init_value = grid.cell_counts.sum()*grid.dx;
        phi.resize(grid.cell_counts.prod());
        for(int i=0;i<phi.size();++i)
            phi[i] = init_value;
        //trace the closest segment of each cell
        std::vector<int> closest_seg;closest_seg.resize(grid.cell_counts.prod(), -1);
        //trace the ray intersection at cell(i,j), the ray start from (0,j) in ith direction
        std::vector<int> intersection_count;intersection_count.resize(grid.cell_counts.prod(), 0);


        const std::vector<Vec<2> >& x = mesh.Vertices();
        for(int s=0;s<mesh.elements.size();++s){
            if(!mesher.Valid_Element(s))continue;
            int p,q;
            p = mesh.elements[s](0);
            q = mesh.elements[s](1);
            Veci<2> fp = grid.Cell_Coord(x[p]);
            Veci<2> fq = grid.Cell_Coord(x[q]);
            int i0 = std::clamp(std::min(fp(0), fq(0))-narrow_band, 0, grid.cell_counts(0)-1);
            int i1 = std::clamp(std::max(fp(0), fq(0))+narrow_band, 0, grid.cell_counts(0)-1);
            int j0 = std::clamp(std::min(fp(1), fq(1))-narrow_band, 0, grid.cell_counts(1)-1);
            int j1 = std::clamp(std::max(fp(1), fq(1))+narrow_band, 0, grid.cell_counts(1)-1);

            for(int j=j0;j<=j1;++j){
                for(int i=i0;i<=i1;++i){
                    Vec<2> gx = grid.Center(Veci<2>(i,j));
                    double dist = Point_Segment_Distance<2>(gx, x[p], x[q]);
                    int idx = grid.Cell_Index(Veci<2>(i,j));
                    if(dist < phi[idx]){
                        phi[idx] = dist;
                        closest_seg[idx] = s;
                    }
                }
            }

            // intersection count, ray send in axis direction
            int tj0 = std::clamp(std::min(fp(1-axis),fq(1-axis)), 0, grid.cell_counts(1-axis)-1);
            int tj1 = std::clamp(std::max(fp(1-axis),fq(1-axis)), 0, grid.cell_counts(1-axis)-1);
            for(int tj=tj0; tj<=tj1; ++tj){
                double a;
                Veci<2> ftj = Veci<2>::Zero();
                ftj(1-axis) = tj;
                double ptj = grid.Center(ftj)(1-axis);//get tj's position
                if((ptj>=x[p](1-axis) && ptj<x[q](1-axis)) || (ptj>=x[q](1-axis) && ptj<x[p](1-axis))){
                    a = std::fabs((ptj-x[q](1-axis))/(x[p](1-axis)-x[q](1-axis)));
                    a = std::clamp(a, 0.0, 1.0);
                    double pti = a*x[p](axis) + (1-a)*x[q](axis);//get intersection point's ith coordinate
                    if(pti > grid.domain_min(axis) && pti < grid.domain_max(axis)){
                        Vec<2> gt = Vec<2>::Zero(); gt(axis) = pti; gt(1-axis)=ptj;
                        Veci<2> icell = grid.Cell_Coord(gt);
                        if(pti > grid.Center(icell)(axis)){
                            //intersection happened after the center of icell
                            icell(axis) += 1;//change the sign of next cell
                        }
                        if(!grid.Valid_Cell(icell))continue;
                        int idx = grid.Cell_Index(icell);
                        ++intersection_count[idx];
                    }
                }
            }
        }

        //determine sign from intersection
        for(int j=0; j<grid.cell_counts(1-axis); ++j){
            int total_count=0;
            for(int i=0; i<grid.cell_counts(axis); ++i){
                Veci<2> cell = Veci<2>::Zero(); cell(axis)=i; cell(1-axis)=j;
                int idx = grid.Cell_Index(cell);
                total_count += intersection_count[idx];
                //from bottom to top, bottom's sign is negative
                if(total_count%2 == 0){
                    phi[idx] = -phi[idx];
                }
            }
        }
    }

    template<>
    void Mesh_To_SDF_Directional_Narrow_Band<3>(const Codimension::CodimMesh<3>& mesh, const Codimension::CodimMesher<3>& mesher, const Grid<3>& grid, std::vector<double>& phi, int axis, int narrow_band){
        //initialize phi to infinity
        double init_value = grid.cell_counts.sum()*grid.dx;
        phi.resize(grid.cell_counts.prod());
        for(int i=0;i<phi.size();++i)
            phi[i] = init_value;
        //trace the closest triangle of each cell
        std::vector<int> closest_tri(grid.cell_counts.prod(), -1);
        //trace the ray intersection at cell(i,j,k), the ray start from (0,j,k) in ith direction
        std::vector<int> intersection_count(grid.cell_counts.prod(), 0);

        int ai = axis;
        int aj = (axis+1)%3;
        int ak = (axis+2)%3;
        const std::vector<Vec<3> >& x = mesh.Vertices();
        for(int t=0;t<mesh.elements.size();++t){
            if(!mesher.Valid_Element(t))continue;
            int p,q,r;
            p = mesh.elements[t](0);
            q = mesh.elements[t](1);
            r = mesh.elements[t](2);
            Veci<3> fp = grid.Cell_Coord(x[p]);
            Veci<3> fq = grid.Cell_Coord(x[q]);
            Veci<3> fr = grid.Cell_Coord(x[r]);
            int i0 = std::clamp(Min(fp(0), fq(0), fr(0))-narrow_band, 0, grid.cell_counts(0)-1);
            int i1 = std::clamp(Max(fp(0), fq(0), fr(0))+narrow_band, 0, grid.cell_counts(0)-1);
            int j0 = std::clamp(Min(fp(1), fq(1), fr(1))-narrow_band, 0, grid.cell_counts(1)-1);
            int j1 = std::clamp(Max(fp(1), fq(1), fr(1))+narrow_band, 0, grid.cell_counts(1)-1);
            int k0 = std::clamp(Min(fp(2), fq(2), fr(2))-narrow_band, 0, grid.cell_counts(2)-1);
            int k1 = std::clamp(Max(fp(2), fq(2), fr(2))+narrow_band, 0, grid.cell_counts(2)-1);

            for(int k=k0;k<=k1;++k){
                for(int j=j0;j<=j1;++j){
                    for(int i=i0;i<=i1;++i){
                        Vec<3> gx = grid.Center(Veci<3>(i,j,k));
                        double dist = Point_Triangle_Distance(gx, x[p], x[q], x[r]);
                        int idx = grid.Cell_Index(Veci<3>(i,j,k));
                        if(dist < phi[idx]){
                            phi[idx] = dist;
                            closest_tri[idx] = t;
                        }
                    }
                }
            }
        }

        // intersection count
        for(int t=0;t<mesh.elements.size();++t){
            if(!mesher.Valid_Element(t))continue;
            int p,q,r;
            p = mesh.elements[t](0);
            q = mesh.elements[t](1);
            r = mesh.elements[t](2);
            Veci<3> fp = grid.Cell_Coord(x[p]);
            Veci<3> fq = grid.Cell_Coord(x[q]);
            Veci<3> fr = grid.Cell_Coord(x[r]);
            int j0 = std::clamp(Min(fp(aj),fq(aj),fr(aj)), 0, grid.cell_counts(aj)-1);
            int j1 = std::clamp(Max(fp(aj),fq(aj),fr(aj)), 0, grid.cell_counts(aj)-1);
            int k0 = std::clamp(Min(fp(ak),fq(ak),fr(ak)), 0, grid.cell_counts(ak)-1);
            int k1 = std::clamp(Max(fp(ak),fq(ak),fr(ak)), 0, grid.cell_counts(ak)-1);
            for(int k=k0; k<=k1; ++k){
                for(int j=j0; j<=j1; ++j){
                    double a, b, c;
                    Veci<3> fjk = Veci<3>::Zero(); fjk(aj) = j; fjk(ak) = k;
                    Vec<3> ljk = grid.Center(fjk);//get j,k's position
                    if(Point_In_Triangle_2D(Vec<2>(ljk(aj), ljk(ak)),
                                            Vec<2>(x[p](aj),x[p](ak)), Vec<2>(x[q](aj),x[q](ak)), Vec<2>(x[r](aj),x[r](ak)), a, b, c)){
                        double li = a*x[p](ai) + b*x[q](ai) + c*x[r](ai);//get intersection point's ith coordinate
                        Vec<3> gt = Vec<3>::Zero(); gt(ai) = li; gt(aj) = ljk(aj); gt(ak) = ljk(ak);
                        Veci<3> icell = grid.Cell_Coord(gt);
                        if(li > grid.Center(icell)(ai)){
                            //intersection happened after the center of icell
                            icell(ai) += 1;//change the sign of next cell
                        }
                        if(!grid.Valid_Cell(icell))continue;
                        int idx = grid.Cell_Index(icell);
                        ++intersection_count[idx];//ray along x coordinate from (0,j,k)'s center intersect with this tri at icell
                    }
                }
            }
        }

        //determine sign from intersection
#ifdef _USE_OMP_
#pragma omp parallel for collapse(2)
#endif
        for(int k=0;k<grid.cell_counts(ak);++k){
            for(int j=0;j<grid.cell_counts(aj);++j){
                int total_count=0;
                for(int i=0; i<grid.cell_counts(ai); ++i){
                    Veci<3> cell = Veci<3>::Zero(); cell(ai) = i; cell(aj) = j; cell(ak) = k;
                    int idx = grid.Cell_Index(cell);
                    total_count += intersection_count[idx];
                    if(total_count%2 == 0){//(0,0,0)'s sign is negative
                        phi[idx] = -phi[idx];
                    }
                }
            }
        }
    }

    template<>
    void Mesh_To_SDF_Thin_Shell<2>(const SurfaceMesh<2>& mesh, const Grid<2>& grid, std::vector<double>& phi, int axis, int sign, int narrow_band){};

    template<>
    void Mesh_To_SDF_Thin_Shell<3>(const SurfaceMesh<3>& mesh, const Grid<3>& grid, std::vector<double>& phi, int axis, int sign, int narrow_band){
        //initialize phi to infinity
        double init_value = grid.cell_counts.sum()*grid.dx;
        phi.resize(grid.cell_counts.prod());
        for(int i=0;i<phi.size();++i)
            phi[i] = init_value;
        //trace the closest triangle of each cell
        std::vector<int> closest_tri;closest_tri.resize(grid.cell_counts.prod(), -1);


        const std::vector<Vec<3> >& x = mesh.Vertices();
        for(int t=0;t<mesh.elements.size();++t){
            int p,q,r;
            p = mesh.elements[t](0);
            q = mesh.elements[t](1);
            r = mesh.elements[t](2);
            Veci<3> fp = grid.Cell_Coord(x[p]);
            Veci<3> fq = grid.Cell_Coord(x[q]);
            Veci<3> fr = grid.Cell_Coord(x[r]);
            int i0 = std::clamp(Min(fp(0), fq(0), fr(0))-narrow_band, 0, grid.cell_counts(0)-1);
            int i1 = std::clamp(Max(fp(0), fq(0), fr(0))+narrow_band, 0, grid.cell_counts(0)-1);
            int j0 = std::clamp(Min(fp(1), fq(1), fr(1))-narrow_band, 0, grid.cell_counts(1)-1);
            int j1 = std::clamp(Max(fp(1), fq(1), fr(1))+narrow_band, 0, grid.cell_counts(1)-1);
            int k0 = std::clamp(Min(fp(2), fq(2), fr(2))-narrow_band, 0, grid.cell_counts(2)-1);
            int k1 = std::clamp(Max(fp(2), fq(2), fr(2))+narrow_band, 0, grid.cell_counts(2)-1);

            for(int k=k0;k<=k1;++k){
                for(int j=j0;j<=j1;++j){
                    for(int i=i0;i<=i1;++i){
                        Vec<3> gx = grid.Center(Veci<3>(i,j,k));
                        double dist = Point_Triangle_Distance(gx, x[p], x[q], x[r]);
                        int idx = grid.Cell_Index(Veci<3>(i,j,k));
                        if(dist < phi[idx]){
                            phi[idx] = dist;
                            closest_tri[idx] = t;
                        }
                    }
                }
            }
        }

        //update i0's phi using i1
        auto check_neighbor = [&] (const Veci<3>& i0, const Veci<3>& i1) {
            if(closest_tri[grid.Cell_Index(i1)] >= 0){
                int ctri = closest_tri[grid.Cell_Index(i1)];
                int tp = mesh.elements[ctri](0);
                int tq = mesh.elements[ctri](1);
                int tr = mesh.elements[ctri](2);
                Vec<3> gx = grid.Center(i0);
                double td = Point_Triangle_Distance(gx, x[tp], x[tq], x[tr]);
                if(td < phi[grid.Cell_Index(i0)]){
                    phi[grid.Cell_Index(i0)] = td;
                    closest_tri[grid.Cell_Index(i0)] = closest_tri[grid.Cell_Index(i1)];
                }
            }
        };

        //sweep in (di,dj,dk) direction
        auto sweep = [&] (int di, int dj, int dk) {
            int ti0, ti1;
            if(di>0) {ti0=1;ti1=grid.cell_counts(0);}
            else {ti0=grid.cell_counts(0)-2;ti1=-1;}
            int tj0, tj1;
            if(dj>0) {tj0=1;tj1=grid.cell_counts(1);}
            else {tj0=grid.cell_counts(1)-2;tj1=-1;}
            int tk0, tk1;
            if(dk>0) {tk0=1;tk1=grid.cell_counts(2);}
            else {tk0=grid.cell_counts(2)-2;tk1=-1;}
            for(int tk=tk0; tk!=tk1; tk+=dk){
                for(int tj=tj0; tj!=tj1; tj+=dj){
                    for(int ti=ti0; ti!=ti1; ti+=di){
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti-di,tj,tk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti,tj-dj,tk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti-di,tj-dj,tk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti,tj,tk-dk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti-di,tj,tk-dk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti,tj-dj,tk-dk));
                        check_neighbor(Veci<3>(ti,tj,tk), Veci<3>(ti-di,tj-dj,tk-dk));
                    }
                }
            }
        };

        //sweep in all direction
        for(int pass=0; pass<2; ++pass){
            sweep(+1,+1,+1);
            sweep(-1,-1,-1);
            sweep(+1,+1,-1);
            sweep(-1,-1,+1);
            sweep(+1,-1,+1);
            sweep(-1,+1,-1);
            sweep(+1,-1,-1);
            sweep(-1,+1,+1);
        }

        //determine sign from closest_tri
        for(int k=0;k<grid.cell_counts(2);++k){
            for(int j=0;j<grid.cell_counts(1);++j){
                for(int i=0; i<grid.cell_counts(0); ++i){
                    int idx = grid.Cell_Index(Veci<3>(i,j,k));
                    Vec<3> center = grid.Center(Veci<3>(i,j,k));
                    const Veci<3>& tri = mesh.Elements()[closest_tri[idx]];
                    Vec<3> d1 = mesh.Vertices()[tri[1]] - mesh.Vertices()[tri[0]];
                    Vec<3> d2 = mesh.Vertices()[tri[2]] - mesh.Vertices()[tri[0]];
                    Vec<3> norm = d1.cross(d2);
                    if((center-mesh.Vertices()[tri[0]]).dot(norm)>0){
                        phi[idx] = -phi[idx];
                    }
                }
            }
        }
    };
}