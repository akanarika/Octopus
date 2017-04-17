import wblut.math.*;
import wblut.processing.*;
import wblut.core.*;
import wblut.hemesh.*;
import wblut.geom.*;


float[][] points;
int numpoints;
HE_Mesh container;
HE_MeshCollection cells;
int numcells;

WB_Render render;

void setup() {
  fullScreen(P3D);
  smooth(8);
  
  //generate points
  numpoints=4000;
  points=new float[numpoints][3];
  for(int i=0;i<numpoints;i++) {
    points[i][0]=random(-300,300);
    points[i][1]=random(-300,300);
    points[i][2]=random(-300,300);
  }
  
  // generate voronoi cells
  HEMC_VoronoiSpheres multiCreator=new HEMC_VoronoiSpheres();
  multiCreator.setPoints(points);
  // alternatively points can be WB_Point[], any Collection<WB_Point> and double[][];
  multiCreator.setN(100);//number of points, can be smaller than number of points in input. 
  multiCreator.setLevel(3);// subdivision level for cell spheres
  multiCreator.setCutoff(200);// maximum radius of cell
   multiCreator.setApprox(true);// approximate cells by point expansion or precise cells by sphere slicing
  multiCreator.setNumTracers(1000);// random points per cell in approcimate mode
  multiCreator.setTraceStep(1);// step size for random points expansion
  multiCreator.setOffset(10);
  cells=multiCreator.create();
  numcells=cells.size();
  
  render=new WB_Render(this);
}

void draw() {
  background(55);
  directionalLight(255, 255, 255, 1, 1, -1);
  directionalLight(127, 127, 127, -1, -1, 1);
  translate(width/2, height/2, 0);
  rotateY(mouseX*1.0f/width*TWO_PI);
  rotateX(mouseY*1.0f/height*TWO_PI);
  drawFaces();
  drawEdges();
}

void drawEdges(){
  stroke(0,50);
  render.drawEdges(cells);
}

void drawFaces(){
  
  noStroke();
  fill(255);
  render.drawFaces(cells);

}