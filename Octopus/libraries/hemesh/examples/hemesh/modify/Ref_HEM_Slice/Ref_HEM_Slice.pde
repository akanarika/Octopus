import wblut.math.*;
import wblut.processing.*;
import wblut.core.*;
import wblut.hemesh.*;
import wblut.geom.*;

HE_Mesh mesh;
WB_Render render;
WB_Plane P;
HEM_Slice modifier;
void setup() {
  size(1000, 1000, P3D);
  createMesh();
  modifier=new HEM_Slice();
  P=new WB_Plane(0,0,0,1,1,1); 
  modifier.setPlane(P);// Cut plane 
  //you can also pass directly as origin and normal:  modifier.setPlane(0,0,-200,0,0,1)
  modifier.setOffset(0);// shift cut plane along normal
  modifier.setCap(true);// cap holes
  modifier.setReverse(false);// keep other side of plane
  mesh.modify(modifier);
  
  render=new WB_Render(this);
}

void draw() {
  background(120);
  directionalLight(255, 255, 255, 1, 1, -1);
  directionalLight(127, 127, 127, -1, -1, 1);
  translate(width/2, height/2);
  rotateY(mouseX*1.0f/width*TWO_PI);
  rotateX(0.25*TWO_PI);
  fill(255);
  noStroke();
  render.drawFaces(mesh);
  fill(255,0,0);
  render.drawFaces(modifier.capFaces);
  fill(0,0,255);
  render.drawFaces(modifier.cutFaces);
  noFill();
  stroke(0);
  render.drawEdges(mesh);
  stroke(255,0,0);
  render.drawPlane(P,300);

}


void createMesh(){
  HEC_Cylinder creator=new HEC_Cylinder();
  creator.setFacets(32).setSteps(16).setRadius(50).setHeight(400);
  mesh=new HE_Mesh(creator);
  
}