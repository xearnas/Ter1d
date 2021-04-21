#include "Mesh.h"
#include <fstream>

using namespace std;
using namespace Eigen;


//////////----------------------------NOS CLASSES MESH--------------------------------


Mesh1D::Mesh1D(DataFile* data_file)
: _df(data_file)
{
}

Mesh2D::Mesh2D(DataFile* data_file)
: _df(data_file)
{
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//--------------------------MESH 1D-------------------------------
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void Mesh1D::Initialize1D()
{
  cout << "Mesh got 1 dimension" << endl;
  _noeuds.resize(_df->Get_Nx()+1,1);
  _centres.resize(_df->Get_Nx(),1);
  for (int i=0; i<_df->Get_Nx() ; i++)
  {
    _noeuds(i,0)=_df->Get_xmin() + i* _df->Get_dx();
    _centres(i,0)= _noeuds(i,0) + _df->Get_dx()/2.;
  }
  _noeuds(_df->Get_Nx(),0) = _df->Get_xmax();
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//--------------------------MESH 2D-------------------------------
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void Mesh2D::Initialize2D()
{
  cout << "Mesh got 2 dimensions" << endl;
    ifstream mesh_file(_df->Get_mesh_root().data());
    if (!mesh_file.is_open())
    {
        cout << "Unable to open file " << _df->Get_mesh_root() << endl;
        exit(0);
    }
    else
    {
        cout << "-------------------------------------------------" << endl;
        cout << "Reading mesh: " << _df->Get_mesh_root() << endl;
    }

    string file_line;
    vector<Edge> edges_boundary;
    int dim = 3;

    while (!mesh_file.eof())
    {
        getline(mesh_file, file_line);
        if (file_line.find("Dimension") != std::string::npos)
        {
            mesh_file >> dim;
        }
        else if (file_line.find("Vertices") != std::string::npos)
        {
          //ON LIT LE NOMBRE DE SOMMETS
            int nb_vertices(0);
            mesh_file >> nb_vertices;
            cout << "Number of vertices  (" << nb_vertices << ")" << endl;
            _vertices.resize(nb_vertices);
            for (int i = 0 ; i < nb_vertices ; ++i)
            {
                double x,y,z;
                int ref;
                mesh_file >> x >> y >> z >> ref;
                _vertices[i] = Vertex(x, y, ref);
            }
        }
        else if (file_line.find("Edges") != std::string::npos)
        {
          // ON LIT LE NOMBRE DE BORDS
            int nb_edges(0);
            mesh_file >> nb_edges;
            cout << "Number of edges (" << nb_edges << ")" << endl;
            edges_boundary.resize(nb_edges);
            int n1, n2, ref;
            for (int i = 0 ; i < nb_edges ; ++i)
            {
                mesh_file >> n1 >> n2 >> ref;
                n1--;
                n2--;
                /*std::string BC_type("none");
                cout << "BC_Ref.size() = " << _BC_ref.size() << endl;
                for (int i=0 ; i < _BC_ref.size() ; i++)
                {
                    if (ref == _BC_ref[i])
                    {
                        BC_type = _BC_type[i];
                    }
                }
                if (BC_type == "none")
                {
                    cout << "Problem with BC in your mesh (reference or type are wrong)" << endl;
                    exit(0);
                }*/
                edges_boundary[i] = Edge(n1, n2, ref, "BC_free");
            }
        }
        else if (file_line.find("Triangles") != std::string::npos)
        {
          // ON LIT LE NOMBRE DE TRIANGLES
            int nb_triangles(0);
            mesh_file >> nb_triangles;
            cout << "Number of triangles (" << nb_triangles << ")" << endl;
            _triangles.resize(nb_triangles);
            for (int i = 0 ; i < nb_triangles ; ++i)
            {
                int vertex1, vertex2, vertex3, ref;
                mesh_file >> vertex1 >> vertex2 >> vertex3 >> ref;
                vertex1--;
                vertex2--;
                vertex3--;
                _triangles[i] = Triangle(vertex1, vertex2, vertex3, ref);
            }
        }
    }

    cout << "---------Edges and Associated Triangles----------" << endl;
    // Toutes les aretes exterieures du maillage sont presentes
    ///NOMBRE TOTAL D'ARETES
    int nb_edges = (3*_triangles.size() + edges_boundary.size())/2;
    _edges.resize(nb_edges);

    int nb_vertices = _vertices.size();
    vector<int> head_minv(nb_vertices, -1);
    vector<int> next_edge(nb_edges, -1);

    // on rajoute d'abord les aretes du bord
    nb_edges = 0;
    for (int i = 0; i < edges_boundary.size(); i++)
    {
        this->Add_single_edge(edges_boundary[i], -1, head_minv, next_edge, nb_edges);
    }

    // ensuite les aretes interieures
    for (int i = 0; i < _triangles.size(); i++)
    {
        const Eigen::Vector3i& nv = _triangles[i].Get_vertices();
        for (int j = 0; j < 3; j++)
        {
            Edge edge(nv(j), nv((j+1)%3), 0, "none");
            Add_single_edge(edge, i, head_minv, next_edge, nb_edges);
        }
    }

    cout << "-----------Triangles center and area-------------" << endl;
    Build_triangles_center_and_area();

    cout << "------------Edges Normal --------------" << endl;
    Build_edges_normal_length_and_center();

    cout << "-------------------------------------------------" << endl;
}

void Mesh2D::Build_triangles_center_and_area()
{
    _tri_center.resize(_triangles.size(),2);
    _tri_area.resize(_triangles.size());
    _tri_h.resize(_triangles.size());

    for (int i = 0; i < _triangles.size(); i++)
    {
        int n1 = _triangles[i].Get_vertices()(0);
        int n2 = _triangles[i].Get_vertices()(1);
        int n3 = _triangles[i].Get_vertices()(2);

        double x1 = _vertices[n1].Get_coor()(0), y1 = _vertices[n1].Get_coor()(1);
        double x2 = _vertices[n2].Get_coor()(0), y2 = _vertices[n2].Get_coor()(1);
        double x3 = _vertices[n3].Get_coor()(0), y3 = _vertices[n3].Get_coor()(1);

        _tri_area(i) = 0.5*fabs((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
        _tri_h(i) = (sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))
                     +sqrt((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))
                     +sqrt((x3-x1)*(x3-x1)+(y3-y1)*(y3-y1)))/3.;
        // centre du triangle
        _tri_center(i,0) = (x1 + x2 + x3)/3.0;
        _tri_center(i,1) = (y1 + y2 + y3)/3.0;
    }
}

void Mesh2D::Build_edges_normal_length_and_center()
{
    _edg_coord.resize(2);
    _edg_coord[0].resize(_edges.size(),2);
    _edg_coord[1].resize(_edges.size(),2);

    _edg_center.resize(_edges.size(),2);
    _edg_normal.resize(_edges.size(),2);
    _edg_length.resize(_edges.size());

    Eigen::Vector2d diff;

    for (int i = 0; i < _edges.size(); i++)
    {
        int t1 = _edges[i].Get_T1();
        int n1 = _edges[i].Get_vertices()(0);
        int n2 = _edges[i].Get_vertices()(1);

        if (t1 >= 0)
        {
            double x1 = _vertices[n1].Get_coor()(0), y1 = _vertices[n1].Get_coor()(1);
            double x2 = _vertices[n2].Get_coor()(0), y2 = _vertices[n2].Get_coor()(1);

            _edg_coord[0](i,0) = x1;
            _edg_coord[0](i,1) = y1;
            _edg_coord[1](i,0) = x2;
            _edg_coord[1](i,1) = y2;
            _edg_length(i) = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));

            // centre de l'arete
            _edg_center(i,0) = 0.5*(x1+x2);
            _edg_center(i,1) = 0.5*(y1+y2);

            // vecteur entre le centre de l'arete et le centre de l'element e1
            diff = _edg_center.row(i) - _tri_center.row(t1);

            // normale suivant un sens arbitraire
            _edg_normal(i,0) = y1 - y2;
            _edg_normal(i,1) = x2 - x1;

            double scal = diff(0)*_edg_normal(i,0) + diff(1)*_edg_normal(i,1);
            if (scal < 0)
            {
                // on change le signe de la normale
                // pour avoir une normale sortante de e1 vers e2
                _edg_normal.row(i) = -_edg_normal.row(i);
            }
            _edg_normal.row(i) /= _edg_length(i);
        }
    }
}

// methode interne qui rajoute une arete
void Mesh2D::Add_single_edge(const Edge& edge, int ne, vector<int>& head_minv, vector<int>& next_edge, int& nb_edges)
{
    int n1 = edge.Get_vertices()(0);
    int n2 = edge.Get_vertices()(1);
    int ref = edge.Get_reference();
    std::string BC_type = edge.Get_BC();

    bool exist = false;
    // we look at the list of edges leaving from n1
    // if we find the same edge than n1->n2 we add the edge
    for (int e = head_minv[n1]; e != -1; e = next_edge[e])
    {
        if (_edges[e].Get_vertices()(1) == n2)
        {
            if (ne >= 0)
            {
                _edges[e].Add_triangle(ne);
            }
            exist = true;
        }
    }

    // if the edge has not been found, we create it
    if (!exist)
    {
        // we initialize the edge
        _edges[nb_edges] = Edge(n1, n2, ref, BC_type);
        if (ne >= 0)
        {
            _edges[nb_edges].Add_triangle(ne);
        }
        // we update the arrays next_edge and head_minv
        next_edge[nb_edges] = head_minv[n1];
        head_minv[n1] = nb_edges;
        nb_edges++;
    }
}

//::::::::::::::::::::::::::CLASSES VERTEX ET TRIANGLES POUR LE 2D:::::::::::::::::::::::::::::
//-------------------------------------------SOMMETS------------------------------------
Vertex::Vertex()
{
    _v_coor[0] = -10000;
    _v_coor[1] = -10000;
    _ref = -1;
}

Vertex::Vertex(double x, double y, int ref) : _ref(ref)
{
    _v_coor[0] = x;
    _v_coor[1] = y;
}

void Vertex::Print() const
{
    cout << "[x, y] = [" << _v_coor[0] << " " << _v_coor[1] << "];" << endl;
    cout << "ref = " << _ref << endl;
}

//-------------------------------------------BORDS------------------------------------

Edge::Edge()
{
    _v_edge[0] = -1;
    _v_edge[1] = -1;
    _ref = -1;
}

Edge::Edge(int vertex1, int vertex2, int ref, std::string BC) : _ref(ref), _BC(BC)
{
    // sort
    if (vertex1 > vertex2)
    {
        _v_edge[0] = vertex2;
        _v_edge[1] = vertex1;
    }
    else
    {
        _v_edge[0] = vertex1;
        _v_edge[1] = vertex2;
    }
    _t1 = -1;
    _t2 = -1;
}


void Edge::Print() const
{
    cout << "[pt1, pt2] = [" << _v_edge[0] << " " << _v_edge[1] << "];" << endl;
    cout << "[t1, t2] = [" << _t1 << " " << _t2 << "];" << endl;
    cout << "ref = " << _ref << endl;
}

//-------------------------------------------TRIANGLES------------------------------------

Triangle::Triangle()
{
    _v_triangle[0] = -1;
    _v_triangle[1] = -1;
    _v_triangle[2] = -1;
    _ref = -1;
}

Triangle::Triangle(int vertex1, int vertex2, int vertex3, int ref) : _ref(ref)
{
    _v_triangle[0] = vertex1;
    _v_triangle[1] = vertex2;
    _v_triangle[2] = vertex3;
}

void Triangle::Print() const
{
    cout << "[pt1, pt2, pt3] = [" << _v_triangle[0] << " " << _v_triangle[1] << " " << _v_triangle[2] << "];" << endl;
    cout << "ref = " << _ref << endl;
}
