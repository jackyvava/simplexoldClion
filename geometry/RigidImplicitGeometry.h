#pragma once

#include "RigidGeometry.h"

/**
 * @brief implicit geometry class, note that surface mesh should be precomputed while initialization
 * @tparam d 
 */
template <int d> class RigidImplicitGeometry : public RigidGeometry<d>{};

template <> class RigidImplicitGeometry<2> : public RigidGeometry<2>
{
public:
    virtual void Precompute_Mesh(double dx);
    virtual double Get_Phi(const Vec<2>& pos) = 0;
    virtual double Get_Mass(const double density) = 0;
    virtual double Get_Inertia_Body(const double density) = 0;
    virtual void Project_Out(Vec<2>& pos) = 0;
    virtual void Project_On(Vec<2>& pos) = 0;
    virtual void Project_Out_Directional(Vec<2>& pos, const Vec<2>& direction) = 0;
    virtual Vec<2> Get_Normal(const Vec<2>& pos) = 0;
};

template <> class RigidImplicitGeometry<3> : public RigidGeometry<3>
{
public:
    virtual void Precompute_Mesh(double dx);
    virtual double Get_Phi(const Vec<3>& pos) = 0;
    virtual double Get_Mass(const double density) = 0;
    virtual Mat<3> Get_Inertia_Body(const double density) = 0;
    virtual void Project_Out(Vec<3>& pos) = 0;
    virtual void Project_On(Vec<3>& pos) = 0;
    virtual void Project_Out_Directional(Vec<3>& pos, const Vec<3>& direction) = 0;
    virtual Vec<3> Get_Normal(const Vec<3>& pos) = 0;
};

template<int d> class SphereGeometry : public RigidImplicitGeometry<d>{};

template<> class SphereGeometry<2> : public RigidImplicitGeometry<2>
{
public:
    double r;
    SphereGeometry(const double _r = 1.0);
    void Initialize(const double _r);
    SphereGeometry& operator = (const SphereGeometry& copy);
    SphereGeometry(const SphereGeometry& copy);
    virtual double Get_Phi(const Vec<2>& pos) override;
    virtual double Get_Mass(const double density)override;
    virtual double Get_Inertia_Body(const double density) override;
    virtual void Project_Out(Vec<2>& pos) override;
    virtual void Project_On(Vec<2>& pos) override;
    virtual void Project_Out_Directional(Vec<2>& pos, const Vec<2>& direction) override;
    virtual Vec<2> Get_Normal(const Vec<2>& pos) override;
};

template<> class SphereGeometry<3> : public RigidImplicitGeometry<3>
{
public:
    double r;
    SphereGeometry(const double _r = 1.0);
    void Initialize(const double _r);
    SphereGeometry& operator = (const SphereGeometry& copy);
    virtual double Get_Phi(const Vec<3>& pos) override;
    virtual double Get_Mass(const double density)override;
    virtual Mat<3> Get_Inertia_Body(const double density)override;
    virtual void Project_Out(Vec<3>& pos) override;
    virtual void Project_On(Vec<3>& pos) override;
    virtual void Project_Out_Directional(Vec<3>& pos, const Vec<3>& direction) override;
    virtual Vec<3> Get_Normal(const Vec<3>& pos) override;
};


template <int d> class BoxGeometry : RigidImplicitGeometry<d>{};

template<> class BoxGeometry<2> : public RigidImplicitGeometry<2>
{
public:
    //a=x b=y
    double a, b;
    BoxGeometry(const double _a = 1.0, const double _b = 1.0);
    void Initialize(const double _a = 1.0, const double _b = 1.0);
    BoxGeometry& operator =(const BoxGeometry& copy);
    virtual double Get_Phi(const Vec<2>& pos) override;
    virtual double Get_Mass(const double density) override;
    virtual double Get_Inertia_Body(const double density) override;
    virtual void Project_Out(Vec<2>& pos) override;
    virtual void Project_On(Vec<2>& pos) override;
    virtual void Project_Out_Directional(Vec<2>& pos, const Vec<2>& direction) override;
    virtual Vec<2> Get_Normal(const Vec<2>& pos) override;
};

template<> class BoxGeometry<3> : public RigidImplicitGeometry<3>
{
public:
    //a=x b=y c=z
    double a, b, c;
    BoxGeometry(const double _a = 1.0, const double _b = 1.0, const double _c = 1.0);
    void Initialize(const double _a = 1.0, const double _b = 1.0, const double _c = 1.0);
    BoxGeometry& operator =(const BoxGeometry& copy);
    virtual double Get_Phi(const Vec<3>& pos) override;
    virtual double Get_Mass(const double density) override;
    virtual Mat<3> Get_Inertia_Body(const double density) override;
    virtual void Project_Out(Vec<3>& pos) override;
    virtual void Project_On(Vec<3>& pos) override;
    virtual void Project_Out_Directional(Vec<3>& pos, const Vec<3>& direction) override;
    virtual Vec<3> Get_Normal(const Vec<3>& pos) override;
};

template <int d> class CylinderGeometry : public RigidImplicitGeometry<d>{};

template<> class CylinderGeometry<2> : public RigidImplicitGeometry<2>
{
public:
    double r, h;
    CylinderGeometry(double _r = 1.0, double _h = 1.0);
    void Initialize(double _r = 1.0, double _h = 1.0);
    CylinderGeometry& operator = (const CylinderGeometry& copy);
    virtual double Get_Phi(const Vec<2>& pos) override;
    virtual double Get_Mass(const double density) override;
    virtual double Get_Inertia_Body(const double density) override;
    virtual void Project_Out(Vec<2>& pos) override;
    virtual void Project_On(Vec<2>& pos) override;
    virtual void Project_Out_Directional(Vec<2>& pos, const Vec<2>& direction) override;
    virtual Vec<2> Get_Normal(const Vec<2>& pos) override;
};

template<> class CylinderGeometry<3> : public RigidImplicitGeometry<3>
{
public:
    double r, h;
    CylinderGeometry(double _r = 1.0, double _h = 1.0);
    void Initialize(double _r = 1.0, double _h = 1.0);
    CylinderGeometry& operator = (const CylinderGeometry& copy);
    virtual double Get_Phi(const Vec<3>& pos) override;
    virtual double Get_Mass(const double density) override;
    virtual Mat<3> Get_Inertia_Body(const double density) override;
    virtual void Project_Out(Vec<3>& pos) override;
    virtual void Project_On(Vec<3>& pos) override;
    virtual void Project_Out_Directional(Vec<3>& pos, const Vec<3>& direction) override;
    virtual Vec<3> Get_Normal(const Vec<3>& pos) override;
};