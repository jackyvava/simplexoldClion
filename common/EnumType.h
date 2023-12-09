#pragma once

enum class IntegratorType : unsigned char
{
    Euler,
    RK1 = Euler,
    Midpoint,
    RK2 = Midpoint,
    Heun,
    Modified_Euler = Heun,
    RK4,
    Backward_Euler,
    Trapezoidal_Rule
};

enum class SolverType : unsigned char
{
    CG,				//conjugate gradient
    DPCG,			//diagonal preconditioned conjugate gradient
    ICPCG,			//incomplete cholesky preconditioned conjugate gradient
    GMG,			//geometry multigrid
    MGPCG,			//multigrid preconditioned conjugate gradient
    BiCGSTAB,		//biconjugate gradient stabilized
};

enum class GridType : unsigned short
{

    Air,
    Fluid,
    Solid,			// if use multiple solids, can assign value (unsigned char)GridType::Solid + i to solid i
};

enum class GeometryType : unsigned short
{
    Duplicated,
    Explicit,
    Shell,
    Box,
    Sphere,
    Cylinder, // in 2d it's just a box, in 3d it's cylinder
};

enum class JointType : unsigned short
{
    Fixed,
    Revolution,
    Sphere,
    Translation,
    Free,
};

enum class Ptype : unsigned short
{
    Free,
    Contact,
    Inside,
    Triple
};