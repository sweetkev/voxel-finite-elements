#include "voxel_image.hpp"
#include "ppm.hpp"
#include "mfem.hpp"
#include <iostream>
#include <unordered_map>
#include <tuple>

using namespace mfem;
using namespace std;

VoxelImage::VoxelImage(string file_) : file(file_) {}

Mesh VoxelImage::MakeMesh() {
    string file_type = GetFileType();

    if (file_type == "pgm") {
        return ReadPgm();   
    }
    else if (file_type == "nhdr") {
        return ReadNrrd();
    }
    else if (file_type == "nrrd") {
        MFEM_ABORT("Must pass in .nhdr header file for .nrrd files");
    }
    else {
        MFEM_ABORT("Unsupported file type");
    }
}

bool VoxelImage::operator()(int i, int j) const
{
   MFEM_ASSERT(cell_dims.size() == 2, "2D access on a non-2D image");
   MFEM_ASSERT(0 <= i && i < cell_dims[0] && 0 <= j && j < cell_dims[1],
               "cell index out of range");
   return data[i + cell_dims[0]*j];
}

bool VoxelImage::operator()(int i, int j, int k) const
{
   MFEM_ASSERT(cell_dims.size() == 3, "3D access on a non-3D image");
   MFEM_ASSERT(0 <= i && i < cell_dims[0] && 0 <= j && j < cell_dims[1] &&
               0 <= k && k < cell_dims[2], "cell index out of range");
   return data[i + cell_dims[0]*(j + cell_dims[1]*k)];
}

bool VoxelImage::operator[](int i) const
{
   MFEM_ASSERT(0 <= i && i < (int)data.size(), "cell index out of range");
   return data[i];
}


string VoxelImage::GetFileType() {
   const auto slash = file.find_last_of('/');
   const auto dot = file.find_last_of('.');
   // No extension, or the last '.' is in a directory name (e.g. "./data/img")
   if (dot == string::npos || (slash != string::npos && dot < slash))
   {
      return "";
   }
   return file.substr(dot + 1);    
}

void ReadComments(istream &in)
{
   string buf;
   in >> std::ws; // absorb any white space
   while (in.peek() == '#')
   {
      std::getline(in,buf);
   }
   in >> std::ws; // absorb any white space
}

Mesh VoxelImage::ReadPgm() {
    cell_dims.resize(2);

    ifstream in(file);
    if (!in.good())
    {
        MFEM_ABORT("Image file not found: " << file);
    }

    // Read magic number.
    char c;
    int p;
    in >> c >> p; // Read magic number which should be P2 or P5
    MFEM_VERIFY(c == 'P' && (p == 2 || p == 5),
               "Invalid PGM file! Unrecognized magic number\""
               << c << p << "\".");
    ReadComments(in);

    // Read dimensions.
    in >> cell_dims[0];
    ReadComments(in);
    in >> cell_dims[1];
    ReadComments(in);

    // Read depth.
    int depth;
    in >> depth;
    ReadComments(in);

    // Read pgm
    const int n = cell_dims[0]*cell_dims[1];
    data.assign(n, false);

    int value;
    for (int i = 0; i < n; ++i)
    {
        in >> value;
        if (value > 0) { data[i] = true; }
    }

    // Create mesh
    //dimension of domain and ambient space
    int dim = 2, sdim = 2;

    //number of vertices, elements, and boundary elements. Will be allocated by FinalizeMesh() later
    int nv = (cell_dims[0] + 1)*(cell_dims[1] + 1);
    int ne = 0;
    for (int i = 0; i < cell_dims[0]*cell_dims[1]; ++i)
    {
        if (data[i]) { ++ne; }
    }
    int nb = 0;

    //initialize mesh
    Mesh mesh(dim, nv, ne, nb, sdim);

    // Add all vertices in the grid, even if they are unused. We will remove
    // unused vertices later.
    real_t h = 1.0 / cell_dims[0];
    for (int j = 0; j < cell_dims[1] + 1; ++j)
    {
        for (int i = 0; i < cell_dims[0] + 1; ++i)
        {
            mesh.AddVertex(i*h, j*h);
        }
    }

    // Add the elements.
    int e = 0;
    for (int j = 0; j < cell_dims[1]; ++j)
    {
        for (int i = 0; i < cell_dims[0]; ++i)
        {
            if (operator()(i,j))
            {
                const int v1 = j*(cell_dims[0] + 1) + i;
                const int v2 = v1 + 1;
                const int v3 = v1 + cell_dims[0] + 1 + 1;
                const int v4 = v3 - 1;
                mesh.AddQuad(v1, v2, v3, v4);
                e++;
            }
        }
   }

   mesh.RemoveUnusedVertices();
   mesh.FinalizeMesh(); // The boundary elements will be generated here.

   return mesh;
}

struct LexIndex
{
   using Idx = array<int, 3>;

   int ndim = 0;
   Idx dims = {};
   int size;

   LexIndex(const vector<int> &dims_) : ndim(dims_.size())
   {
      copy(dims_.begin(), dims_.end(), dims.begin());
      size = 1;
      for (int d = 0; d < ndim; ++d) { size *= dims[d]; }
   }

   Idx CartesianIndex(int i) const
   {
      Idx cart;
      for (int d = 0; d < ndim; ++d)
      {
         cart[d] = i % dims[d];
         i /= dims[d];
      }
      return cart;
   }

   int LinearIndex(const Idx &i) const
   {
      int stride = 1;
      int idx = 0;
      for (int d = 0; d < ndim; ++d)
      {
         idx += i[d] * stride;
         stride *= dims[d];
      }
      return idx;
   }

   int Size() const { return size; }
};

Mesh VoxelImage::ReadNrrd()
{
   string header_file = file;
   map<string,string> header;

   ifstream f(header_file);
   string line;
   getline(f, line); // Skip first line
   while (getline(f, line))
   {
      if (line.front() == '#') { continue; }
      auto pos = line.find(':');
      if (pos != string::npos)
      {
         const string field_name = line.substr(0, pos);
         string value = line.substr(pos + 2, string::npos);
         if (value.back() == '\r')
         {
            value.resize(value.size() - 1);
         }
         header[field_name] = value;
      }
   }

   // Data file path is relative to the directory containing the header
   string data_file = header["data file"];
   if (data_file.front() != '/')
   {
      const auto slash = header_file.find_last_of('/');
      if (slash != string::npos)
      {
         data_file = header_file.substr(0, slash + 1) + data_file;
      }
   }
   const int dim = stoi(header["dimension"]);
   const string sizes_str = header["sizes"];
   string directions = header["space directions"];

   vector<int> vert_dims(dim);
   cell_dims.resize(dim);
   {
      istringstream s(sizes_str);
      for (int d = 0; d < dim; ++d)
      {
         s >> cell_dims[d];
         vert_dims[d] = cell_dims[d] + 1;
      }
   }

   Vector phys_dims(dim);
   {
      for (char &c : directions)
      {
         if (c == '(' || c == ',' || c == ')')
         {
            c = ' ';
         }
      }
      istringstream s(directions);

      for (int d = 0; d < dim; ++d)
      {
         double val;
         for (int d2 = 0; d2 < d; ++d2) { s >> val; }
         s >> phys_dims[d];
         for (int d2 = d + 1; d2 < dim; ++d2) { s >> val; }
      }
   }

   int total_cells = 1;
   for (int d = 0; d < dim; ++d) { total_cells *= cell_dims[d]; }
   data.assign(total_cells, false);
   int total_verts = 1;
   for (int d = 0; d < dim; ++d) { total_verts *= vert_dims[d]; }

   vector<uint8_t> cells(total_cells);
   {
      ifstream f(data_file, ios::binary);
      f.read(reinterpret_cast<char*>(cells.data()), total_cells);
   }

   LexIndex cell_idx(cell_dims);
   LexIndex vert_idx(vert_dims);
   LexIndex offset_idx(vector<int>(dim, 2));

   vector<int> lex2cell(total_cells, -1);
   vector<int> lex2vert(total_verts, -1);
   vector<int> cell2lex;
   vector<int> vert2lex;

   cell2lex.reserve(total_cells);
   vert2lex.reserve(total_verts);

   int cell_i = 0;
   int vert_i = 0;
   for (int i = 0; i < total_cells; ++i)
   {
      if (cells[i] != 0)
      {
         data[i] = true;
         lex2cell[i] = cell_i;
         cell2lex.push_back(i);
         ++cell_i;

         auto cart = cell_idx.CartesianIndex(i);
         for (int o = 0; o < offset_idx.Size(); ++o)
         {
            auto cart_offset = cart;
            const auto offset = offset_idx.CartesianIndex(o);
            for (int d = 0; d < dim; ++d) { cart_offset[d] += offset[d]; }
            const int iv_lex = vert_idx.LinearIndex(cart_offset);
            if (lex2vert[iv_lex] < 0)
            {
               lex2vert[iv_lex] = vert_i;
               vert2lex.push_back(iv_lex);
               ++vert_i;
            }
         }
      }
   }

   Mesh mesh(dim, vert_i, cell_i);

   for (int i = 0; i < vert_i; ++i)
   {
      const int lex = vert2lex[i];
      const auto cart = vert_idx.CartesianIndex(lex);
      array<double,3> coords;
      for (int d = 0; d < dim; ++d)
      {
         coords[d] = cart[d]*phys_dims[d] / double(cell_dims[d]);
      }
      mesh.AddVertex(coords.data());
   }

   for (int i = 0; i < cell_i; ++i)
   {
      if (dim == 2)
      {
         MFEM_ABORT("");
      }
      else if (dim == 3)
      {
         array<int,8> verts;
         const int cell_lex = cell2lex[i];
         const auto cell_cart = cell_idx.CartesianIndex(cell_lex);
         const int vert_lex = vert_idx.LinearIndex(cell_cart);
         verts[0] = lex2vert[vert_lex];
         verts[1] = lex2vert[vert_lex + 1];
         verts[2] = lex2vert[vert_lex + 1 + vert_dims[0]];
         verts[3] = lex2vert[vert_lex + vert_dims[0]];

         const int off = vert_dims[0] * vert_dims[1];

         verts[4] = lex2vert[vert_lex + off];
         verts[5] = lex2vert[vert_lex + off + 1];
         verts[6] = lex2vert[vert_lex + off + 1 + vert_dims[0]];
         verts[7] = lex2vert[vert_lex + off + vert_dims[0]];

         mesh.AddHex(verts.data());
      }
      else
      {
         MFEM_ABORT("");
      }
   }

   mesh.FinalizeMesh();

   return mesh;
}