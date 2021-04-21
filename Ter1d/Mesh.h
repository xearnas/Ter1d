#include "Function.h"
#include "DataFile.h"

using namespace std;
using namespace Eigen;
#pragma once

class Vertex
{
private:
    Vector2d _v_coor;
    int _ref;
public:
    Vertex();
    Vertex(double x, double y, int ref);
    void Print() const;
    const Vector2d Get_coor() const
    {
        return _v_coor;
    };
};

class Edge
{
private:
    Vector2i _v_edge;
    int _ref;
    int _t1, _t2;
    string _BC;
public:
    Edge();
    Edge(int vertex1, int vertex2, int ref, string BC);
    void Print() const;
    void Add_triangle(int t)
    {
        if (_t1 == -1)
            _t1 = t;
        else
            _t2 = t;
    }
    const Vector2i& Get_vertices() const
    {
        return _v_edge;
    }
    int Get_T1() const
    {
        return _t1;
    };
    int Get_T2() const
    {
        return _t2;
    };
    int Get_reference() const
    {
        return _ref;
    };
    string Get_BC() const
    {
        return _BC;
    };
};

class Triangle
{
private:
    Vector3i _v_triangle;
    int _ref;
public:
    Triangle();
    Triangle(int vertex1, int vertex2, int vertex3, int ref);
    void Print() const;
    const Vector3i& Get_vertices() const
    {
        return _v_triangle;
    }
};


/*class Mesh
{
protected:
  DataFile* _df;
  string _name_mesh;
  string _file_name;

public:
  Mesh(DataFile* data_file);
  virtual ~Mesh();



  virtual void Initialize()=0;
  virtual MatrixXd Get_noeuds()=0;
  virtual MatrixXd Get_centres()=0;
  virtual const vector<Vertex> & Get_vertices() const =0;
  virtual const vector<Triangle> & Get_triangles() const=0;
  virtual const VectorXd & Get_triangles_area() const=0;
  virtual const VectorXd & Get_triangles_length() const=0;
  virtual const vector<Edge> & Get_edges() const=0;
  virtual const VectorXd & Get_edges_length() const=0;
  virtual const MatrixXd & Get_edges_normal() const=0;
  virtual const vector<Matrix<double, Dynamic, 2> > & Get_edges_coord() const=0;

};*/

class Mesh1D
{
private:
  DataFile* _df;
  VectorXd _noeuds;
  VectorXd _centres;
public:
  Mesh1D(DataFile* data_file);
  void Initialize1D();
  MatrixXd Get_noeuds()
  {
    return _noeuds;
  };
  MatrixXd Get_centres()
  {
    return _centres;
  };
};

class Mesh2D
{
private:
  DataFile* _df;
  // liste de tous les sommets
  vector<Vertex> _vertices;
  // liste de tous les triangles
  vector<Triangle> _triangles;
  // centre de tous les triangles
  MatrixXd _tri_center;
  // aire de tous les triangles
  VectorXd _tri_area;
  // dx de tous les triangles
  VectorXd _tri_h;
  // liste de toutes les arêtes
  vector<Edge> _edges;
  // liste de toutes les normales unitaires !!!
  MatrixXd _edg_normal;
  // liste de toutes les longueurs d'arêtes
  VectorXd _edg_length;
  // centre des aretes
  MatrixXd _edg_center;
  // coordonnées des arrêtes
  vector<Matrix<double, Dynamic, 2> > _edg_coord;
  const vector<int> _BC_ref;
  // vecteur de type
  const vector<string> _BC_type;
public:
  Mesh2D(DataFile* data_file);
  void Initialize2D();
  void Build_triangles_center_and_area();
  void Build_edges_normal_length_and_center();
  void Add_single_edge(const Edge& edge, int ne, vector<int>& head_minv,vector<int>& next_edge, int& nb_edges);
  ///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  //--------------------------------GET TOUTES LES INFOS 2D--------------------------------------
  ///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  const vector<Vertex> & Get_vertices() const
  {
      return _vertices;
  };

  const vector<Triangle> & Get_triangles() const
  {
      return _triangles;
  };
  const MatrixXd & Get_triangles_center() const
  {
      return _tri_center;
  };
  const VectorXd & Get_triangles_area() const
  {
      return _tri_area;
  };
  const VectorXd & Get_triangles_length() const
  {
      return _tri_h;
  };
  const vector<Edge> & Get_edges() const
  {
      return _edges;
  };
  const VectorXd & Get_edges_length() const
  {
      return _edg_length;
  };
  const MatrixXd & Get_edges_normal() const
  {
      return _edg_normal;
  };
  const MatrixXd & Get_edges_center() const
  {
      return _edg_center;
  };
  const vector<Matrix<double, Dynamic, 2> > & Get_edges_coord() const
  {
      return _edg_coord;
  };
};
