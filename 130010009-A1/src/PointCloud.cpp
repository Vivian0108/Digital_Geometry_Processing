#include "PointCloud.hpp"
#include "DGP/Graphics/Shader.hpp"
#include "PointKDTree.hpp"
#include <fstream>
#include <sstream>
// #include "Hyperplane3.hpp"

PointCloud::PointCloud(std::vector<Point> const & points_)
: points(points_)
{
  recomputeAABB();
}

PointCloud::PointCloud(std::vector<Vector3> const & positions, std::vector<Vector3> const & normals)
{
  alwaysAssertM(positions.size() == normals.size(), "PointCloud: Number of positions != number of normals");

  for (size_t i = 0; i < positions.size(); ++i)
    points.push_back(Point(positions[i], normals[i]));

  recomputeAABB();
}

bool
PointCloud::load(std::string const & path)
{
  // Simple format: Each line is either
  //   x y z
  //    OR
  //   x y z nx ny nz
  //
  // where (nx, ny, nz) is the normal

  std::ifstream in(path.c_str());
  if (!in)
  {
    DGP_ERROR << "Could not open file for reading: " << path;
    return false;
  }

  std::string line;
  while (getline(in, line))
  {
    // Skip empty lines
    line = trimWhitespace(line);
    if (line.empty())
      continue;

    std::istringstream line_in(line);
    Vector3 p;
    if (!(line_in >> p[0] >> p[1] >> p[2]))
    {
      DGP_ERROR << "Could not read point " << points.size() << " from line: " << line;
      return false;
    }

    // Normal is optional
    Vector3 n;
    if (!(line_in >> n[0] >> n[1] >> n[2]))  // doesn't properly handle malformed lines, but we'll ignore this for now
      n = Vector3::zero();

    points.push_back(Point(p, n));
  }

  recomputeAABB();

  return true;
}

bool
PointCloud::save(std::string const & path) const
{
  std::ofstream out(path.c_str(), std::ios::binary);
  if (!out)
  {
    DGP_ERROR << "Could not open file for writing: " << path;
    return false;
  }

  for (size_t i = 0; i < points.size(); ++i)
  {
    Vector3 const & p = points[i].getPosition();
    Vector3 const & n = points[i].getNormal();
    out << p[0] << ' ' << p[1] << ' ' << p[2] << ' ' << n[0] << ' ' << n[1] << ' ' << n[2] << '\n';
  }

  return true;
}

Graphics::Shader *
createPointShader(Graphics::RenderSystem & rs)
{
  static std::string const VERTEX_SHADER =
"void main()\n"
"{\n"
"  gl_Position = ftransform();\n"
"  gl_FrontColor = gl_Color;\n"
"  gl_BackColor = gl_Color;\n"
"}\n";

  static std::string const FRAGMENT_SHADER =
"void main()\n"
"{\n"
"  gl_FragColor = gl_Color;\n"
"}\n";

  Graphics::Shader * shader = rs.createShader("Point Graphics::Shader");
  if (!shader)
    throw Error("Could not create point shader");

  // Will throw errors on failure
  shader->attachModuleFromString(Graphics::Shader::ModuleType::VERTEX, VERTEX_SHADER.c_str());
  shader->attachModuleFromString(Graphics::Shader::ModuleType::FRAGMENT, FRAGMENT_SHADER.c_str());

  return shader;
}

void
PointCloud::draw(Graphics::RenderSystem & rs, Real normal_len) const
{
  // Make this static to ensure just one shader is created. Assumes rendersystem is constant, not the best design pattern.
  static Graphics::Shader * shader = createPointShader(rs);

  rs.pushShader();
  rs.pushColorFlags();
  rs.pushShapeFlags();

    rs.setShader(shader);
    rs.setPointSize(2.0f);

    Vector3 lo = bbox.getLow();
    Vector3 ext = bbox.getExtent();

    rs.beginPrimitive(Graphics::RenderSystem::Primitive::POINTS);
      for (size_t i = 0; i < points.size(); ++i)
      {
        Vector3 rel_pos = (points[i].getPosition() - lo) / ext;
        rs.setColor(ColorRGB(rel_pos[0], rel_pos[1], rel_pos[2]));
        rs.sendVertex(points[i].getPosition());
      }
    rs.endPrimitive();

    if (normal_len > 0)
    {
      rs.setColor(ColorRGB(0.5, 0.5, 1.0));  // blue

      rs.beginPrimitive(Graphics::RenderSystem::Primitive::LINES);
        for (size_t i = 0; i < points.size(); ++i)
        {
          Vector3 const & p = points[i].getPosition();
          Vector3 const & n = points[i].getNormal();

          rs.sendVertex(p);
          rs.sendVertex(p + normal_len * n);
        }
      rs.endPrimitive();
    }

  rs.popShapeFlags();
  rs.popColorFlags();
  rs.popShader();
}

void
PointCloud::recomputeAABB()
{
  bbox.setNull();

  for (size_t i = 0; i < points.size(); ++i)
    bbox.merge(points[i].getPosition());
}

long
PointCloud::ransac(long num_iters, Real slab_thickness, long min_points, Slab & slab, std::vector<Point *> & slab_points) const
{
  // TODO
  int enab_count = 0;
  // PointKDTree kd_tree;
  // std::vector<Point *> enab_points(points.size());
  std::vector<Point *> enab_points;             //pointer for the enabled points
  srand(100);
  for (size_t j = 0;j<points.size();j++){
    if(points[j].isEnabled()){
      // enab_count +=1;
      enab_points.push_back(const_cast<Point *>(&points[j]));
    }
  }
  

  // srand((unsigned)time(0));
  long random_trip[3];  //array containing 3 random integers
  long count_pts;
  int max_pts =0;
  Plane3 best_fit_plane;  
  Plane3 fit_plane;

  for(size_t i = 0;i < num_iters ;i++){
    for(size_t j=0;j<3;j++){
      random_trip[j] = rand() % enab_points.size();
    }
    fit_plane = Plane3::fromThreePoints(enab_points[random_trip[0]]->getPosition(),
      enab_points[random_trip[1]]->getPosition(),enab_points[random_trip[2]]->getPosition()); //creating a plane from random triplets
    Slab slab1;
    std::vector<Point *> slab1_points;
    slab1.setPlane(fit_plane);            //creating a slab for that plane
    slab1.setThickness(slab_thickness);
      
    count_pts = 0;
    for(size_t k =0;k<enab_points.size();k++){
      if (slab1.contains(enab_points[k]->getPosition())){
            count_pts+=1;
            slab1_points.push_back(enab_points[k]);  //counting the number of points inside the slab
      }
    }

    /*PointKDTree kd_enab_tree(enab_points);    
    kd_enab_tree.rangeQuery(slab1,slab1_points);
    count_pts = slab1_points.size();*/

    if(count_pts >=min_points && count_pts > max_pts){
        max_pts = count_pts;
        best_fit_plane = fit_plane;
        slab_points.clear();
        for(long j =0;j<slab1_points.size();j++){
          slab_points.push_back(slab1_points[j]);
        }
        
        slab.setPlane(best_fit_plane);    //setting the best fit plane for the slab
        slab.setThickness(slab_thickness);
    }    

  }
slab.updateCorners(slab_points);
// slab.draw(Graphics::RenderSystem rs, ColorRGBA const & plane_color, ColorRGBA const & outline_color);

  //   - Construct a kd-tree on the enabled points (remember to build the kd-tree with pointers to existing points -- you
  //     shouldn't be copying the points themselves, either explicitly or implicitly).
  //   - Generate num_iters random triplets of enabled points and fit a plane to them.
  //   - Using the kd-tree, see how many other enabled points are contained in the slab supported by this plane with thickness
  //     slab_thickness (extends to distance 0.5 * slab_thickness on each side of the plane).
  //   - If this number is >= min_points and > the previous maximum, the plane is the current best fit. Set the 'slab' argument
  //     to be the slab for this plane, and update slab_points to be the set of (enabled) matching points for this plane.
  //   - At the end, for visualization purposes, update the corners of the best slab using its set of matching points, and
  //     return the number of (enabled) matching points.

  return count_pts;
}

long
PointCloud::ransacMultiple(long num_planes, long num_iters, Real slab_thickness, long min_points, std::vector<Slab> & slabs)
const
{
  for (size_t i = 0; i < points.size(); ++i)
    points[i].setEnabled(true);

  slabs.clear();
  for (long i = 0; i < num_planes; ++i)
  {
    Slab slab;
    std::vector<Point *> slab_points;

    long num_matching_pts = ransac(num_iters, slab_thickness, min_points, slab, slab_points);
    if (num_matching_pts <= 0)
      break;

    slabs.push_back(slab);

    // Don't consider these points in subsequent slabs
    for (size_t j = 0; j < slab_points.size(); ++j)
      slab_points[j]->setEnabled(false);
  }

  return (long)slabs.size();
}

void
PointCloud::adaptiveDownsample(std::vector<Slab> const & slabs)
{
  // TODO
}
