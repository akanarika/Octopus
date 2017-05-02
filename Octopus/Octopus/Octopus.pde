import controlP5.*;

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
    DeformableModel model;
    
    DwVertexRecorder vertex_recorder;
    
    HE_Mesh mesh;
    WB_Render3D render;
    ControlP5 cp5;
    
     float Poisson_Ratio = 0.07f;
     float Young_Variable = 4.45f;
     int Young_Modulus = 5;  // Young's Modulus
     float Mass_Variable = 1.0f;
     float Mass_Density = 2.18f;
     float Damping = 6.3f;
     float kv = 10000.0f;  // spring constant   ( Energy = 1/2 kx^2, x = delta_distance)
     int Neighbour_Count = 20;
     boolean MOVE_CAM = false;
     int Time_Step = 24;
    public void settings() {
      // Init the window
      size(viewport_w, viewport_h, P3D);
      smooth(30);
    }
    
    public void setup() {
      
      surface.setLocation(viewport_x, viewport_y);
      //peasyCam = new PeasyCam(this, -4.083,  -6.096,   7.000, 61);
           
      cam = new PeasyCam(this, 0,  0,  0, 0.2);
    
      // projection
      perspective(60 * DEG_TO_RAD, width/(float)height, 2, 5000);
      
      render=new WB_Render3D(this);
      
      // Load a 3D model
      //HE_MESH
      HEC_FromOBJFile creator = new HEC_FromOBJFile();
      creator.setPath("//Users//Mercurian//Octopus//Octopus//Octopus//data//sphere.obj").setScale(1.0);
      mesh=new HE_Mesh(creator);
      octopus = WB_PShapeFactory.createWireframePShape(mesh, this);
 
      //octopus = createShapeQuad(loadShape("sphere.obj"));
      //setup_pixel_flow_monitor();
      
      // set up world
      //model = new DeformableModel(new Vector3D(8, 8, 8), 1.20f, mesh, Poisson_Ratio, pow(10, Young_Modulus), kv, pow(10, Mass_Density), Damping, Neighbour_Count);
      model = null;
      
      // setup skylight renderer
      createGUI();
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
      model.Update();
      render();
      cam.setActive(MOVE_CAM);
      displayGUI();
    
    }
    
    void displayGUI()
    {
      noLights();
      cam.beginHUD();
      cp5.draw();
      cam.endHUD();
    }
    
    public void keyPressed(){
    if(key == CODED){
      if(keyCode == ALT){
        MOVE_CAM = true;
        }
      }
    }
    
     public void keyReleased(){
      MOVE_CAM = false; 
    }
    
    void createGUI()
    {
      cp5 = new ControlP5(this);
      cp5.setAutoDraw(false);
  
      int sx, sy, px, py, oy;
      sx = 100; sy = 50; oy = (int)(sy*0.2f);
      Group group_physics = cp5.addGroup("global");
      {
        cp5.addSlider("Poisson_Ratio")
         .setPosition(sx, sy)
         .setRange(-1, 0.5)
         ;
         
        cp5.addSlider("Young_Variable")
         .setPosition(sx, sy * 2 + oy)
         .setRange(1, 10)
         ;
         
         cp5.addSlider("Young_Modulus")
         .setPosition(sx, sy * 2 + oy * 2)
         .setRange(1, 8)
         ;
         
         cp5.addSlider("Mass_Variable")
         .setPosition(sx, sy * 3 + oy)
         .setRange(1, 10);
         ;
         
         cp5.addSlider("Mass_Density")
         .setPosition(sx, sy * 3 + oy * 2)
         .setRange(1, 10);
         ;
         
                  
         cp5.addSlider("kv")
         .setPosition(sx, sy * 4 + oy * 3)
         .setRange(0, 100000);
         ;
         
         cp5.addSlider("Damping")
         .setPosition(sx, sy * 5 + oy * 4)
         .setRange(0, 10);
         ;
         cp5.addSlider("Neighbour_Count")
         .setPosition(sx, sy * 6 + oy * 5)
         .setRange(1, 20);
         ;
         
         cp5.addButton("restart")
         .setPosition(sx, sy * 7 + oy * 6)
         .setSize(200,19)
         .setValue(0);
         

     ;
         
      }
    }
    
    
    void restart()
    {
      model = new DeformableModel(new Vector3D(8, 8, 8), 1.20f, mesh, Poisson_Ratio, Young_Variable * pow(10, Young_Modulus), 
          kv, Mass_Variable * pow(10, Mass_Density), Damping, Neighbour_Count);
    }
    
    void render()
    {
      background(0);
      lights();
      directionalLight(255, 255, 255, -1, -1, -1);

      //shape(octopus, 0, 0);
      strokeWeight(1);
      beginShape(POINTS);
      //world.draw();
      if(model != null){
        model.draw(phyxelIdx);
      }
      if(keyPressed)
      {
        if(key == ' ')
           phyxelIdx ++;
      }
      endShape();
      
      translate((float)model.worldSize.x/2, (float)model.worldSize.y/2, (float)model.worldSize.z/2);
      noFill();
      box((float)model.worldSize.x, (float)model.worldSize.y, (float)model.worldSize.z);
    }
    
    void render_skylight()
    {
    }

        