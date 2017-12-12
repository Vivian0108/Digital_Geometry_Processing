#include "MeshEdge.hpp"
#include "MeshFace.hpp"
#include "MeshVertex.hpp"
#include "DGP/Vector4.hpp"

MeshEdge *
MeshEdge::nextAroundEndpoint(int i)
{
  debugAssertM(i == 0 || i == 1, "MeshEdge: Invalid endpoint index");

  if (numFaces() > 2)  // non-manifold
    return NULL;

  // Find which incident face has this endpoint as the origin of the edge when stepping round the face. The required edge
  // is then the predecessor of this edge around the face.
  for (FaceIterator fi = facesBegin(); fi != facesEnd(); ++fi)
  {
    Face * face = *fi;
    MeshEdge * prev = face->getPredecessor(this);
    if (prev->hasEndpoint(endpoints[i]))  // found it!
      return prev;
  }

  return NULL;
}

void
MeshEdge::updateQuadricCollapseError()
{
  // TODO

  // Update both quadric_collapse_error and quadric_collapse_position, using the existing endpoint quadrics and the method of
  // Garland/Heckbert.
  //
  // NOTE: Remember to check if the quadric Q' is invertible. If not, you will have to use a fallback option such as the
  // midpoint of the edge (or in the worst case, set the error to a negative value to indicate this edge should not be
  // collapsed).
  
  DMat4 Q1, Q2, Q, _Qinv;
  Q1 = endpoints[0]->getQuadric();
  Q2 = endpoints[1]->getQuadric();
  Q = Q1 + Q2;
  Vector3 trans(Q(0, 3), Q(1, 3), Q(2, 3));
  MatrixMN< 3, 3, double > tempMat = Q.upper3x3();
  DMat4 _Q(tempMat, trans);
  double _Qdet = _Q.determinant();
  if (abs(_Qdet) > 0.001) {
  	_Qinv = _Q.inverse();

	Vector4 V(0.0, 0.0, 0.0, 1.0);
	Vector4 V_final = _Qinv * V;

	quadric_collapse_position[0] = V_final[0]/V_final[3];
	quadric_collapse_position[1] = V_final[1]/V_final[3];
	quadric_collapse_position[2] = V_final[2]/V_final[3];

	quadric_collapse_error = V_final.dot(Q * V_final);
  }

  else if (abs(_Qdet) < 0.001) {
  	Vector3 p1 = endpoints[0]->getPosition();
  	Vector3 p2 = endpoints[1]->getPosition();

  	// Finding midpoint
  	Vector4 V(0.5*(p1[0] + p2[0]), 0.5*(p1[1] + p2[1]), 0.5*(p1[2] + p2[2]), 1.0);

  	quadric_collapse_position[0] = V[0];
  	quadric_collapse_position[1] = V[1];
  	quadric_collapse_position[2] = V[2];

  	quadric_collapse_error = V.dot(Q * V);
  }

  else
  	quadric_collapse_error = -1.0;
}
