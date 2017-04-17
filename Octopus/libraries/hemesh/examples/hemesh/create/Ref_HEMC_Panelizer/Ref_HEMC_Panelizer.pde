import wblut.math.*;
import wblut.processing.*;
import wblut.core.*;
import wblut.hemesh.*;
import wblut.geom.*;


HE_Mesh hull;
HE_MeshCollection panels;


WB_Render render;

void setup() {
  size(1000,1000,P3D);
  smooth(8);
  
  //create a hull mesh
  hull=new HE_Mesh(new HEC_Geodesic().setB(2).setC(0).setRadius(200));
  HET_MeshOp.splitFacesCenter(hull);
  //panelize the hull
  HEMC_Panelizer multiCreator=new HEMC_Panelizer();
  multiCreator.setMesh(hull);
  multiCreator.setThickness(10);
  multiCreator.setOffset(0,100);
  panels=multiCreator.create();

  
  render=new WB_Render(this);
}

void draw() {
  background(55);
  directionalLight(255, 255, 255, 1, 1, -1);
  directionalLight(127, 127, 127, -1, -1, 1);
  translate(width/2, height/2);
  rotateY(mouseX*1.0f/width*TWO_PI);
  rotateX(mouseY*1.0f/height*TWO_PI);
  drawFaces();
  drawEdges();
}

void drawEdges(){
  stroke(0);
  render.drawEdges(panels);
}

void drawFaces(){
  noStroke();
  fill(255);
  render.drawFaces(panels);
}