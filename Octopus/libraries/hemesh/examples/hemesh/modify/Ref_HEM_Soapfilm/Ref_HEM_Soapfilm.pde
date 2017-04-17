//work in progress, artifact-prone. Small and degenerate faces are not being removed yet...

import wblut.math.*;
import wblut.processing.*;
import wblut.core.*;
import wblut.hemesh.*;
import wblut.geom.*;
HE_Mesh mesh;
WB_Render render;
HE_Selection sel;
double A;
double f;
void setup() {
  size(1000,1000,P3D);
  smooth(8);
  HEC_Creator creator=new HEC_Cube(200,1,1,1);

  mesh=new HE_Mesh(creator); 

 HEM_Extrude ext=new HEM_Extrude().setDistance(300);
 mesh.modify(ext);

sel=ext.extruded;
mesh.deleteFaces(sel);
 mesh.subdivide(new HES_Planar(),4);
f=0;
A=mesh.getArea();
  render=new WB_Render(this);
}

void draw() {
  background(55);
  directionalLight(255, 255, 255, 1, 1, -1);
  directionalLight(127, 127, 127, -1, -1, 1);
  translate(width/2, height/2);
  rotateY(mouseX*1.0f/width*TWO_PI);
  rotateX(mouseY*1.0f/height*TWO_PI);
  stroke(0);
  render.drawEdges(mesh);
  noStroke();
  render.drawFaces(mesh);
  if(f<0.9999) update();
}

void update(){
 
 mesh.modify(new HEM_Soapfilm().setIterations(1)); 
 double nA=mesh.getArea();
 f=nA/A;
 println(f);
 A=nA;
}