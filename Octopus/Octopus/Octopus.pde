/****************************************
*
* Octopus
* CMU 15464/15-664 Technical Animation
* main.pde
*
*/

// Libraries

import com.thomasdiewald.pixelflow.java.DwPixelFlow;
import com.thomasdiewald.pixelflow.java.utils.DwBoundingSphere;
import com.thomasdiewald.pixelflow.java.utils.DwVertexRecorder;

import peasy.*;

PShape s;
int screen_width = 800;
int screen_height = 600;

public void settings() {
  // Init the window
  size(screen_width, screen_height, P3D);
}

public void setup() {
  // Set background color
  background(0);
  
  // Load a 3D model
  s = loadShape("sphere.obj");
  
  // record list of vertices of the given shape
  DwVertexRecorder vertex_recorder = new DwVertexRecorder(this, s);
  
  print(vertex_recorder.verts_count);
}

public void draw() {
  background(0);
  camera(mouseX, height/2, (height/2) / tan(PI/6), width/2, height/2, 0, 0, 1, 0);
  lights();
  directionalLight(255, 255, 255, 0, -1, 0);
  translate(width/2, height/2, -100);
  shape(s, 0, 0, 100, 100);
}