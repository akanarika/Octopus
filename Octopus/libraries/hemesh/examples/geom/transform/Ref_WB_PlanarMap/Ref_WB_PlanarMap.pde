import wblut.geom.*;
import wblut.processing.*;

WB_RandomPoint source;
WB_Render3D render;
WB_Point[] points;
int numPoints;

WB_PlanarMap XY;
WB_PlanarMap YZ;
WB_PlanarMap XZ;
WB_PlanarMap P;

WB_Triangulation2D triangulationXY;
WB_Triangulation2D triangulationYZ;
WB_Triangulation2D triangulationXZ;
WB_Triangulation2D triangulationP;

void setup() {
  size(1000, 1000, P3D);
  source=new WB_RandomOnSphere().setRadius(250);
  render=new WB_Render3D(this);
  numPoints=500;
  points=new WB_Point[numPoints];
  for (int i=0; i<numPoints; i++) {
    points[i]=source.nextPoint();
  }

  XY=new WB_PlanarMap(WB_PlanarMap.XY, 270);
  YZ=new WB_PlanarMap(WB_PlanarMap.YZ, 270);
  XZ=new WB_PlanarMap(WB_PlanarMap.XZ, 270);
  WB_Plane plane=new WB_Plane(0, 0, 0, 1, 1, 1);
  P=new WB_PlanarMap(plane);
  
  // WB_PlanarMap is a  coordinate system transform, it converts 3D world coordinates to local 3D coordinates in the system of a plane, and vice versa. 
  // This is for example useful to 2D triangulate a polygon with an arbitrary orientation in space.
  // WB_PlanarMap is not identical to projection since the origin of 3D space is mapped to the origin of the plane.
  // (Projection would map the origin to the closest point on the plane, not always the origin of that plane.)
  // Furthermore, WB_PlanarMap is reversible while a projection is not.
  
  triangulationXY=WB_Triangulate.triangulate2D(points, XY);
  triangulationYZ=WB_Triangulate.triangulate2D(points, YZ);
  triangulationXZ=WB_Triangulate.triangulate2D(points, XZ);
  triangulationP=WB_Triangulate.triangulate2D(points, P);
}


void draw() {
  background(55);
  directionalLight(255, 255, 255, 1, 1, -1);
  directionalLight(127, 127, 127, -1, -1, 1);
  translate(width/2, height/2-50);
  rotateY(mouseX*1.0f/width*PI+QUARTER_PI);

  noFill();
  stroke(0);
  render.drawPoint(points, 2);
  stroke(255, 0, 0);

  //WB_PlanarMap is a example of a WB_Map, more specifically a WB_Map2D, a mapping intended to transform between
  //3D and 2D space.
  
  //Render functions labelled with "mapped" apply a WB_Map before rendering, by definition this will render
  //any WB_Map2D to the XY-plane.
  //Render functions labelled "unmapped" do the reverse, unmapping a point before rendering it.
  //To show the results of a WB_Map2D embedded in 3D space, use render functions labelled with "embedded2D". This is a special case of
  // the "unmapped" functions that ignore the local z-coordinate (for example the elevation in a WB_Planar map). 
  
  
  render.drawPointEmbedded2D(points, 2, XY);
  render.drawTriangulationEmbedded2D(triangulationXY, points, XY);
  stroke(0, 255, 0);
  render.drawPointEmbedded2D(points, 2, YZ);
  render.drawTriangulationEmbedded2D(triangulationYZ, points, YZ);
  stroke(0, 0, 255);
  render.drawPointEmbedded2D(points, 2, XZ);
  render.drawTriangulationEmbedded2D(triangulationXZ, points, XZ);
  stroke(0, 255, 255);
  render.drawPointEmbedded2D(points, 2, P);
  render.drawTriangulationEmbedded2D(triangulationP, points, P);
}