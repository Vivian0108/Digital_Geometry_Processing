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
  double v_new[4] = {0.0,0.0,0.0,0.0};
  double result[4] = {0.0,0.0,0.0,0.0};
  DMat4 quadric_new = endpoints[0]->getQuadric() + endpoints[1]->getQuadric();
  DMat4 A_position = quadric_new;
  A_position.setRow(3,{0.0,0.0,0.0,1.0});
  quadric_collapse_error = 0.0;
  if(abs(A_position.determinant()) > 0.001 ){

    (A_position.inverse()).getColumn(3,v_new);
    quadric_collapse_position[0] = v_new[0]/v_new[3];
    quadric_collapse_position[1] = v_new[1]/v_new[3];
    quadric_collapse_position[2] = v_new[2]/v_new[3];
  }
  else{
    quadric_collapse_position = 0.5*(endpoints[0]->getPosition() + endpoints[1]->getPosition());
    v_new[0] = quadric_collapse_position[0];
    v_new[1] = quadric_collapse_position[1];
    v_new[2] = quadric_collapse_position[2];
    v_new[3] = 1;

  }

  //Vector3 new_position = quadric_collapse_position;
  for (long i = 0; i < 4; ++i){
    
    for (long j = 0; j < 4; ++j){
      result[i] += quadric_new.get(i,j) * v_new[j];
    }
  }
  for (long i =0;i<4;i++){
      quadric_collapse_error += result[i]*v_new[i];
  
  }
  
  //quadric_collapse_error = (new_position.dot(quadric_new*new_position));
  // Update both quadric_collapse_error and quadric_collapse_position, using the existing endpoint quadrics and the method of
  // Garland/Heckbert.
  //
  // NOTE: Remember to check if the quadric Q' is invertible. If not, you will have to use a fallback option such as the
  // midpoint of the edge (or in the worst case, set the error to a negative value to indicate this edge should not be
  // collapsed).

}
