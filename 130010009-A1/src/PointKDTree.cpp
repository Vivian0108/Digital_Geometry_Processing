#include "PointKDTree.hpp"

PointKDTree::Node::~Node()
{
  delete lo;
  delete hi;
}

PointKDTree::PointKDTree(std::vector<Point> const & points)
: root(NULL)
{
  build(points);
}

PointKDTree::PointKDTree(std::vector<Point *> const & points)
: root(NULL)
{
  build(points);
}

PointKDTree::PointKDTree(PointCloud const & pcloud)
: root(NULL)
{
  build(pcloud.getPoints());
}

void
PointKDTree::build(std::vector<Point> const & points)
{
  std::vector<Point *> pp(points.size());
  for (size_t i = 0; i < pp.size(); ++i)
    pp[i] = const_cast<Point *>(&points[i]);  // removing the const is not the greatest thing to do, be careful...

  build(pp);
}
void  PointKDTree::getting_kdtree(Node* root1,std::vector<Point *> const & points){

  AxisAlignedBox3 bbox;

    for (size_t i = 0; i < points.size(); ++i){
      bbox.merge(points[i]->getPosition());//created the bounding box for the points
    }
     
  Node* root = root1;
  long index = (bbox.getHigh() - bbox.getLow()).maxAbsAxis();//index for the largest dimension of the box
  Node* lo;
  Node* hi;
  int lo_count = 0;
  int hi_count = 0;
  
  for (size_t i = 0; i < points.size(); ++i){
    
     if(points[i]->getPosition()[index] <= ((bbox.getHigh()[index]+ bbox.getLow()[index])/2)){

          
          lo->points[lo_count] = points[i];//points lesser than the middle plane go to lo node
          ++lo_count;
      
    } 
    
    else{
          
          hi->points[hi_count]=points[i];//points greater than the middle plane go to hi node
          ++hi_count;
      
    }
  }
  root->points.clear();

  if(lo->points.size() > 5){
    getting_kdtree(root,lo->points);//Not a leaf node, hence splitting
    lo->points.clear();
  }
  if(hi->points.size() > 5){
    getting_kdtree(root,hi->points);//not a leaf node, hence splitting
    hi->points.clear();
    }
}

void
PointKDTree::build(std::vector<Point *> const & points)
{
  root = new Node;
  root->points = points;//created the root node
  // TODO

  static size_t const MAX_POINTS_PER_LEAF = 5;
  getting_kdtree(root,points);//creating the kd tree
  
  // A kd-tree is just a binary search tree, and is constructed in a near-identical way.
  //
  // - Initially assign (pointers to) all points to the root node.
  // - Recursively subdivide the points, splitting the parent box in half along the longest axis and creating two child nodes
  //   for each split. Stop when number of points in node <= MAX_POINTS_PER_LEAF.
  // - Don't forget to save space by clearing the arrays of points in internal nodes. Only the leaves need to store references
  //   to points.
    
}