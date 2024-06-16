#pragma once

//#define USE_OpenMesh

#ifdef USE_OpenMesh
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

typedef OpenMesh::PolyMesh_ArrayKernelT<> PolyMesh;
typedef OpenMesh::TriMesh_ArrayKernelT<> TriMesh;

#elif
struct HE_vert;
struct HE_edge;
struct HE_face;
struct HE_edge
{
    HE_vert* vert = nullptr; // vertex at the start of the half-edge
    HE_edge* pair = nullptr; // oppositely oriented adjacent half-edge
    HE_face* face = nullptr; // face the half-edge borders
    HE_edge* next = nullptr; // next half-edge around the face
};
struct HE_vert
{
    double x;
    double y;
    double z;
    double u;
    double v;
    HE_edge* edge = nullptr; // one of the half-edges emantating from the vertex
    int index;               // Index of the vertex
    HE_vert(double _x, double _y, double _z, int _index) : x(_x), y(_y), z(_z), index(_index)
    {
    }
    HE_vert(double _x, double _y, double _z) : x(_x), y(_y), z(_z)
    {
    }
    HE_vert& operator+=(const HE_vert& other)
    {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }
    HE_vert& operator-=(const HE_vert& other)
    {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }
    friend HE_vert operator*(int scalar, const HE_vert& vertex)
    {
        return HE_vert(scalar * vertex.x, scalar * vertex.y, scalar * vertex.z);
    }
    friend HE_vert operator*(double scalar, const HE_vert& vertex)
    {
        return HE_vert(scalar * vertex.x, scalar * vertex.y, scalar * vertex.z);
    }
    friend std::ostream& operator<<(std::ostream& os, const HE_vert& vertex)
    {
        os << "Vertex: (" << vertex.x << ", " << vertex.y << ", " << vertex.z << ")";
        return os;
    }
};
struct HE_face
{
    HE_edge* edge = nullptr; // one of the half-edges bordering the face
};

using namespace std;

#endif // USE_OpenMesh