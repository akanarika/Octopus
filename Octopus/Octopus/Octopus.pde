import wblut.nurbs.*;
import wblut.hemesh.*;
import wblut.core.*;
import wblut.geom.*;
import wblut.processing.*;
import wblut.math.*;

import wblut.nurbs.*;
import wblut.hemesh.*;
import wblut.core.*;
import wblut.geom.*;
import wblut.processing.*;
import wblut.math.*;

/****************************************
 *
 * Octopus
 * CMU 15464/15-664 Technical Animation
 * main.pde
 *
 */

// Libraries
import java.util.Locale;

import com.thomasdiewald.pixelflow.java.DwPixelFlow;
import com.thomasdiewald.pixelflow.java.render.skylight.DwSceneDisplay;
import com.thomasdiewald.pixelflow.java.render.skylight.DwSkyLight;
import com.thomasdiewald.pixelflow.java.utils.DwBoundingSphere;
import com.thomasdiewald.pixelflow.java.utils.DwVertexRecorder;

import peasy.*;
import processing.core.PApplet;
import processing.core.PGraphics;
import processing.core.PMatrix3D;
import processing.core.PShape;
import processing.opengl.PGraphics3D;
  
  
  
  
  
    int viewport_w = 1280;
    int viewport_h = 720;
    int viewport_x = 230;
    int viewport_y = 0;
   
    
    // Object
    PShape octopus;
    
    // renderer
    DwSkyLight skylight;
    
    // camera
    PeasyCam cam;
    
    // World 
    World world;
    
    DwVertexRecorder vertex_recorder;
    
    HE_Mesh mesh;
    WB_Render3D render;
    
    public void settings() {
      // Init the window
      size(viewport_w, viewport_h, P3D);
      smooth(30);
    }
    
    public void setup() {
      
      surface.setLocation(viewport_x, viewport_y);
      //peasyCam = new PeasyCam(this, -4.083,  -6.096,   7.000, 61);
           
      cam = new PeasyCam(this, 0.000,  1.000,   0.000, 0.2);
    
      // projection
      perspective(60 * DEG_TO_RAD, width/(float)height, 2, 5000);
      
      render=new WB_Render3D(this);
      
      // Load a 3D model
      //HE_MESH
      HEC_FromOBJFile creator = new HEC_FromOBJFile();
      creator.setPath("C://Users//sijiah//Documents//GitHub//Octopus//Octopus//Octopus//data//sphere.obj").setScale(1.0);
      mesh=new HE_Mesh(creator);
      octopus = WB_PShapeFactory.createWireframePShape(mesh, this);
 
      //octopus = createShapeQuad(loadShape("sphere.obj"));
      //setup_pixel_flow_monitor();
      
      // set up world
      world = new World(new Vector3D(8, 8, 8), 0.5f, mesh);
      
      // setup skylight renderer
    }
    
    void setup_pixel_flow_monitor()
    {
      // record list of vertices of the given shape
      vertex_recorder = new DwVertexRecorder(this, octopus); 
      DwBoundingSphere scene_bs = new DwBoundingSphere();
      scene_bs.compute(vertex_recorder.verts, vertex_recorder.verts_count);
      //print(vertex_recorder.verts_count);
    }
    
    void setup_skylight_renderer()
    {
      
    }
    
    
    
    
    
    
    
    
    
    int phyxelIdx = 0;
    public void draw() {
      
      //// Render useing the skylight renderer
      //render_skylight();
      world.Update();
      render();
    
    }
    
    void render()
    {
      background(255);
      lights();
      directionalLight(255, 255, 255, -1, -1, -1);

      shape(octopus, 0, 0);
      strokeWeight(1);
      beginShape(POINTS);
      //world.draw();
      world.draw(phyxelIdx);
      if(keyPressed)
      {
        if(key == ' ')
           phyxelIdx ++;
      }
      endShape();
      
      translate(world.worldSize.x/2, world.worldSize.y/2, world.worldSize.z/2);
      noFill();
      box(world.worldSize.x, world.worldSize.y, world.worldSize.z);
    }
    
    void render_skylight()
    {
    }
    
   
        