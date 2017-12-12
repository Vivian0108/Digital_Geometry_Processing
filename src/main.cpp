#include "Mesh.hpp"
#include "Viewer.hpp"
#include <cstdlib>
#include <math.h>
#include <vector>
#include <list>
#include <algorithm>
#include <string>
#include <stringstream>
#include <fstream>
#include "DGP/Polygon3.hpp"
#include "DGP/Math.hpp"
#include "DGP/Matrix3.hpp"

#define PI 3.1415926536

void readFromCsv (std::string filename, std::vector< std::vector<float> > & rows) {
    std::fstream inFile (filename.c_str(), std::ios::in);
    std::string line;
    while (std::getline(inFile, line)) {
        std::stringstream lStream(line);
        std::vector<float> temp;
        std::string cell;
        while (std::getline(lStream, cell, ','))
            temp.push_back(::atof(cell.c_str()));
        rows.push_back(temp);
    }
}

void samplePoints (Mesh const & mesh, size_t num_points, std::vector<Mesh::Vertex *> & points) {
    points.clear();
    std::vector<float> tri_areas;
    for (Mesh::FaceConstIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi) {
        std::vector<Vector3> face_points;
        for (Mesh::Face::VertexConstIterator fvi = (*fi).verticesBegin(); fvi != (*fi).verticesEnd(); ++fvi)
            face_points.push_back((*fvi)->getPosition());

        Polygon3 tri; // triangle made up of vertices of the given face
        for (std::vector<Vector3>::iterator i = face_points.begin(); i != face_points.end(); ++i)
          tri.addVertex(*i);

        float area = tri.area();
        tri_areas.push_back(area);
    }

    for (size_t i = 0; i < num_points; ++i) {
        float max_area = 0;
        size_t max_area_index1 = 0;
        for (size_t j = 0; j < tri_areas.size(); ++j) {
            if(tri_areas[j] > max_area) {
                max_area = tri_areas[j];
                max_area_index1 = j;
            }
        }

        Mesh::FaceConstIterator fi1 = mesh.facesBegin();
        for (size_t k = 0; k < max_area_index1 ; ++k, ++fi1) { }

        Mesh::Vertex * pi = *((*fi1).verticesBegin());
        points.push_back(pi);
        tri_areas.erase(tri_areas.begin() + max_area_index1);
    }
}

void curvature (std::vector<Mesh::Vertex *> sampledPoints, std::vector< std::vector<float> > & eigenvalues, \
    std::vector< std::vector<Vector3> > & eigenvectors) {

    for (size_t i = 0; i < sampledPoints.size(); ++i) {
        Matrix3 tau = Matrix3::zero();
        float areaSum = 0.0;
        for (Mesh::Vertex::EdgeConstIterator ei = sampledPoints[i]->edgesBegin(); ei != sampledPoints[i]->edgesEnd(); ++ei) {
            std::vector<Vector3> normals;
            for (Mesh::Edge::FaceConstIterator efi = (*ei)->facesBegin(); efi != (*ei)->facesEnd(); ++efi)
                normals.push_back((*efi)->getNormal());

            float beta = normals[0].dot(normals[1]);
            Mesh::Vertex * endpoint_1 = (*ei)->getEndpoint(0);
            Mesh::Vertex * endpoint_2 = (*ei)->getEndpoint(1);
            Vector3 pos_1 = endpoint_1->getPosition();
            Vector3 pos_2 = endpoint_2->getPosition();
            Vector3 edgeVector = pos_2 - pos_1;
            float edgeLength = edgeVector.length();
            Vector3 unitEdgeVector = edgeVector/edgeLength;
            Matrix3 edgeMat = Matrix3::zero();
            edgeMat(0, 0) = unitEdgeVector[0] * unitEdgeVector[0];
            edgeMat(0, 1) = unitEdgeVector[0] * unitEdgeVector[1];
            edgeMat(0, 2) = unitEdgeVector[0] * unitEdgeVector[2];
            edgeMat(1, 0) = unitEdgeVector[1] * unitEdgeVector[0];
            edgeMat(1, 1) = unitEdgeVector[1] * unitEdgeVector[1];
            edgeMat(1, 2) = unitEdgeVector[1] * unitEdgeVector[2];
            edgeMat(2, 0) = unitEdgeVector[2] * unitEdgeVector[0];
            edgeMat(2, 1) = unitEdgeVector[2] * unitEdgeVector[1];
            edgeMat(2, 2) = unitEdgeVector[2] * unitEdgeVector[2];
            edgeMat = beta*edgeLength*edgeMat;
            
            tau = tau + edgeMat;
        }

        for (Mesh::Vertex::FaceConstIterator vfi = sampledPoints[i]->facesBegin(); vfi != sampledPoints[i]->facesEnd(); ++vfi) {
            Polygon3 face;
            for (Mesh::Face::VertexConstIterator i = (*vfi)->verticesBegin(); i != (*vfi)->verticesEnd(); ++i)
              face.addVertex((*i)->getPosition());

            float area = face.area();
            areaSum += area;
        }
        tau = tau/areaSum;

        float egnvalues[3];
        Vector3 egnvectors[3];
        for (size_t i = 0; i < 3; ++i) {
            egnvalues[i] = 0.0f;
            egnvectors[i].set(0.0f, 0.0f, 0.0f);
        }
        tau.eigenSolveSymmetric(egnvalues, egnvectors);

        float minCurvature, maxCurvature;
        Vector3 minCurvDir, maxCurvDir;
        float min = 100000000000.0;
        float max = 0.0;
        for (size_t i = 0; i < 3; ++i) {
            if (abs(egnvalues[i]) < min) 
                min = egnvalues[i];
            if (egnvalues[i] > max) {
                max = egnvalues[i];
                maxCurvature = max;
                minCurvDir = egnvectors[i];
            }
        }
        size_t minCurvIndex = 1;
        for (size_t j = 0; j < 3; ++j) {
            if ((egnvalues[j] == min) || (egnvalues[j] == max))
                continue;
            else
                minCurvIndex = j;
        }
        minCurvature = egnvalues[minCurvIndex];
        maxCurvDir = egnvectors[minCurvIndex];

        Vector3 normal = minCurvDir.cross(maxCurvDir);
        if ((normal.dot(sampledPoints[i]->getNormal())) < 0)
            for (size_t i = 0; i < 3; ++i)
                minCurvDir[i] = -minCurvDir[i];

        std::vector<float> a;
        a.push_back(minCurvature); a.push_back(maxCurvature);
        eigenvalues.push_back(a);

        std::vector<Vector3> b;
        b.push_back(minCurvDir); b.push_back(maxCurvDir);
        eigenvectors.push_back(b);
    }
}

void signatures (std::vector<Mesh::Vertex *> sampledPoints, std::vector< std::vector<float> > eigenvalues, \
    std::vector< std::vector<Vector3> > eigenvectors, std::vector< std::vector <float> > & signatures) {
    signatures.clear();
    for (size_t i = 0; i < sampledPoints.size(); ++i) {
        if (eigenvalues[i][0] == eigenvalues[i][1])
            continue;
        for (size_t j = 0; j < sampledPoints.size(); ++j) {
            if (j == i)
                continue;
            Vector3 normalI = eigenvectors[i][0].cross(eigenvectors[i][1]);
            Vector3 normalJ = eigenvectors[j][0].cross(eigenvectors[j][1]);

            Matrix3 RotMat1 = Matrix3::rotationArc(normalI, normalJ, true);

            Vector3 transformedE1 = RotMat1 * eigenvectors[i][0];
            // Vector3 transformedE2 = RotMat1 * eigenvectors[i][1];

            Matrix3 RotMat2 = Matrix3::rotationArc(transformedE1, eigenvectors[j][0], true);

            Matrix3 R = RotMat2 * RotMat1;

            // for finding out Euler angles, referred to http://nghiaho.com/?page_id=846
            float Rx = atan2 (R(3,2), R(3,3)) * 180 / PI;
            float Ry = atan2 (-R(3,1), sqrt(R(3,2)*R(3,2) + R(3,3)*R(3,3)));
            float Rz = atan2 (R(2,1), R(1,1));

            float s_ij = 0.5 * ((eigenvalues[i][0]/eigenvalues[j][0]) + (eigenvalues[i][1]/eigenvalues[j][1]));

            Vector3 t_ij = sampledPoints[j]->getPosition() - s_ij*(R*(sampledPoints[i]->getPosition()));

            std::vector<float> T_ij;
            T_ij.push_back(s_ij);
            T_ij.push_back(Rx); T_ij.push_back(Ry); T_ij.push_back(Rz);
            T_ij.push_back(t_ij[0]); T_ij.push_back(t_ij[1]); T_ij.push_back(t_ij[2]);

            signatures.push_back(T_ij); 
            

        }
    }
    //write code for reading values from .csv

    readFromCsv(/* filename */, /* 2D vector in which CSV data is to be stored*/std::vector<int> labels);
    

                //check one ring neighbours of Vi1

    
        
    
}

void patching(std::vector<std::vector<float> > cluster_centers,std::vector<Mesh::Vertex *> sampledPoints,\
    std::vector<std::vector<Mesh::Vertex *> > patches,)


int
usage(int argc, char * argv[])
{
  DGP_CONSOLE << "";
  DGP_CONSOLE << "Usage: " << argv[0] << " <mesh-in> [<target-num-faces> [<mesh-out>]]";
  DGP_CONSOLE << "";

  return -1;
}

int
main(int argc, char * argv[])
{
  if (argc < 2)
    return usage(argc, argv);

  std::string in_path = argv[1];

  long target_num_faces = -1;
  std::string out_path;
  if (argc >= 3)
  {
    target_num_faces = std::atoi(argv[2]);

    if (argc >= 4)
      out_path = argv[3];
  }

  Mesh mesh;
  if (!mesh.load(in_path))
    return -1;

  DGP_CONSOLE << "Read mesh '" << mesh.getName() << "' with " << mesh.numVertices() << " vertices, " << mesh.numEdges()
              << " edges and " << mesh.numFaces() << " faces from " << in_path;

  if (target_num_faces >= 0 && mesh.numFaces() > target_num_faces)
  {
    mesh.decimateQuadricEdgeCollapse(target_num_faces);
    mesh.updateBounds();
  }

  if (!out_path.empty())
  {
    if (!mesh.save(out_path))
      return -1;

    DGP_CONSOLE << "Saved mesh to " << out_path;
  }

  Viewer viewer;
  viewer.setObject(&mesh);
  viewer.launch(argc, argv);

  return 0;
}
