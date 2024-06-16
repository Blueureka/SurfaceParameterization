#include <iostream>
#include <string>

#include <HalfEdge.h>
#include <Eigen/Eigen>

#include <igl/opengl/glfw/Viewer.h>

igl::opengl::glfw::Viewer viewer;

typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd Vuv;

using namespace std;

enum ParametrizationMethod
{
    TutteUniform = 0,
    TutteCotan,
    ABF,
    LSCM,
    ARAP,

};

void Tutte(TriMesh& mesh)
{
    int nf = mesh.n_faces();
    int nv = mesh.n_vertices();

    // calc surface area
    double area_sum = 0;
    for (const auto& face : mesh.faces())
    {
        auto itfv = mesh.fv_iter(face);
        auto v0 = mesh.point(itfv);
        itfv++;
        auto v1 = mesh.point(itfv);
        itfv++;
        auto v2 = mesh.point(itfv);

        auto e0 = v1 - v0;
        auto e1 = v2 - v0;

        auto avec = cross(e0, e1);
        area_sum += avec.norm() / 2.0;
    }

    // set the boundary vertices to circle
    int boundary_num = 0;
    auto heit = mesh.halfedges_begin();
    while (!mesh.is_boundary(*heit))
        heit++;
    auto he_start = *heit;
    auto he_cur = he_start;
    do
    {
        he_cur = he_cur.next();
        boundary_num++;
    } while (he_cur != he_start);

    double delta_angle = 2 * M_PI / boundary_num;
    double area_1_factor = sqrt(area_sum / M_PI);

    Eigen::SparseMatrix<double> coffMatrix(nv,nv);
    std::vector<Eigen::Triplet<double>> tripletlist;
    Eigen::VectorXd bu = Eigen::VectorXd::Zero(nv);
    Eigen::VectorXd bv = Eigen::VectorXd::Zero(nv);

    for (int i = 0; i < boundary_num; i++)
    {
        auto v = he_start.to();
        bu(v.idx()) = area_1_factor * cos(i * delta_angle);
        bv(v.idx()) = area_1_factor * sin(i * delta_angle);
        he_start = he_start.next();
    }

    for (const auto& v : mesh.vertices())
    {
        int index = v.idx();
        if (mesh.is_boundary(v))
        {
            tripletlist.emplace_back(index, index, 1);
        }
        else
        {
            int cnt = 0;
            for (const auto& neighbor : mesh.vv_range(v))
            {
                tripletlist.emplace_back(index, neighbor.idx(), -1);
                cnt++;
            }
            tripletlist.emplace_back(index, index, cnt);
        }
    }

    coffMatrix.setFromTriplets(tripletlist.begin(), tripletlist.end());
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(coffMatrix);

    Eigen::VectorXd xu = solver.solve(bu);
    Eigen::VectorXd xv = solver.solve(bv);

    for (auto& v : mesh.vertices())
    {
        int index=v.idx();
        mesh.point(v) = {xu(index), xv(index), 0};
    }
}

TriMesh SurfaceParametrization(const TriMesh& mesh, const ParametrizationMethod method)
{
    TriMesh paraMesh;
    paraMesh.assign(mesh);

    switch (method)
    {
    case ParametrizationMethod::TutteUniform:
        Tutte(paraMesh);
        break;
    case ParametrizationMethod::TutteCotan:
        break;
    case ParametrizationMethod::LSCM:
        break;
    case ParametrizationMethod::ARAP:
        break;
    default:
        break;
    }

    return paraMesh;
}

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
    std::cout << "Key: " << key << " " << (unsigned int)key << std::endl;
    switch (key)
    {
    case '1':
        viewer.data().clear();        // 清空屏幕上的网格
        viewer.data().set_mesh(V, F); // 显示修改后的网格
        viewer.data().set_uv(Vuv);
        viewer.data().show_texture = true;
        // viewer.data().set_colors(C);
        break;
    case '2':
        viewer.data().clear();          // 清空屏幕上的网格
        viewer.data().set_mesh(Vuv, F); // 显示修改后的网格
        viewer.data().set_uv(Vuv);
        viewer.data().show_texture = true;
        // viewer.data().set_colors(C);
        break;
    default:
        break;
    }
    return false;
}

int main(int argc, char** argv)
{
    viewer.callback_key_down = &key_down;
    TriMesh mesh;
    
    // read mesh from stdin
    // std::string readPath = std::string(OFF_DATA_WRITE_PATH) + "hexagon.obj";
    std::string readPath = std::string(OFF_DATA_WRITE_PATH) + "BsplineSurface.obj";
    if (!OpenMesh::IO::read_mesh(mesh, readPath))
    {
        std::cerr << "Error: Cannot read mesh from " << readPath << std::endl;
        return 1;
    }

    //Parametrization
    TriMesh paraMesh = SurfaceParametrization(mesh, ParametrizationMethod::TutteUniform);

    // Count number of vertices and faces
    int num_vertices = mesh.n_vertices();
    int num_faces = mesh.n_faces();

    V.resize(num_vertices, 3);
    Vuv.resize(num_vertices, 3);
    F.resize(num_faces, 3);
    // Fill Eigen MatrixXd with vertex coordinates
    int vertex_index = 0;
    for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
    {
        TriMesh::Point p = mesh.point(*v_it);
        V.row(vertex_index) << p[0], p[1], p[2];
        ++vertex_index;
    }
    vertex_index = 0;
    for (TriMesh::VertexIter v_it = paraMesh.vertices_begin(); v_it != paraMesh.vertices_end(); ++v_it)
    {
        TriMesh::Point p = paraMesh.point(*v_it);
        Vuv.row(vertex_index) << p[0], p[1], p[2];
        ++vertex_index;
    }
    

    // Fill Eigen MatrixXi with face indices
    int face_index = 0;
    for (TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
    {
        int i = 0;
        for (TriMesh::FaceVertexIter fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
        {
            F(face_index, i) = fv_it->idx();
            ++i;
        }
        ++face_index;
    }

    

    viewer.data().set_mesh(V, F);
    viewer.data().double_sided = true;
    viewer.data().show_texture = true;
    viewer.launch();

    // write mesh to stdout
    /*std::string writePath = std::string(OFF_DATA_WRITE_PATH) + "hexagonOutput.off";
    if (!OpenMesh::IO::write_mesh(mesh, writePath))
    {
        std::cerr << "Error: cannot write mesh to " << writePath << std::endl;
        return 1;
    }*/

    return 0;
}
