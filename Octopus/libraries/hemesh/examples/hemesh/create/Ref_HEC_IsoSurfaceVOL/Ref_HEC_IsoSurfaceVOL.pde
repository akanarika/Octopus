import wblut.math.*;
import wblut.processing.*;
import wblut.core.*;
import wblut.hemesh.*;
import wblut.geom.*;
import java.util.*;

HE_Mesh mesh;
WB_Render render;
 List<WB_Tetrahedron> tetra;
void setup() {
  fullScreen(P3D);
  smooth(8);


  float[][][] values=new float[51][51][51];
  for (int i = 0; i < 51; i++) {
    for (int j = 0; j < 51; j++) {
      for (int k = 0; k < 51; k++) {
        values[i][j][k]=2.5*noise(0.07*i, 0.07*j, 0.07*k);
      }
    }
  }

  HEC_IsoSurfaceVOL creator=new HEC_IsoSurfaceVOL();
  creator.setResolution(50, 50,50);// number of cells in x,y,z direction
  creator.setSize(12,12,12);// cell size
  creator.setValues(values);// values corresponding to the grid points
  // values can also be double[][][]
  creator.setIsolevel(0.9,1.1);// isolevel to mesh


  // use creator.clearBoundary() to rest boundary values to "no value".
  // A boundary value of "no value" results in an open mesh

  mesh=new HE_Mesh(creator);
tetra=creator.tetra;
  render=new WB_Render(this);
}

void draw() {
  background(25);
  lights();
  translate(width/2, height/2);
  rotateY(mouseX*1.0f/width*TWO_PI);
  rotateX(mouseY*1.0f/height*TWO_PI);
  noStroke();
  fill(255);
  render.drawTetrahedron(tetra);
 
}