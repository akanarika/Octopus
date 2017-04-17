import wblut.math.*;
import wblut.processing.*;
import wblut.hemesh.*;
import wblut.geom.*;

HE_Mesh mesh;
WB_Render3D render;

void setup() {
  size(1000, 1000, P3D);
  smooth(8);
  HEC_Dodecahedron creator=new HEC_Dodecahedron();creator.setCenter(50,0,10);
  creator.setEdge(200); 
  mesh=new HE_Mesh(creator); 
  
  HET_MeshOp.splitFacesCenter(mesh);
 HET_MeshOp.splitFacesTri(mesh,40);
  mesh.smooth();
  render=new WB_Render3D(this);
}

void draw() {
  background(55);
  directionalLight(255, 255, 255, 1, 1, -1);
  directionalLight(127, 127, 127, -1, -1, 1);
  translate(width/2, height/2);
  rotateY(frameCount*0.005);
  stroke(0);
  render.drawEdges(mesh);
  fill(255);
  noStroke();
  render.drawFaces(mesh);
  fill(255, 0, 0);
  HE_Face f=render.pickClosestFace(mesh,mouseX,mouseY);
  if(f!=null) render.drawFace(f);
}