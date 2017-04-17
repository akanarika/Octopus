import wblut.processing.*;
import wblut.hemesh.*;
import wblut.geom.*;
import wblut.math.*;
import java.util.List;

WB_GeometryFactory gf=new WB_GeometryFactory();
WB_Render2D render;
ArrayList<WB_Point> shell;
ArrayList<WB_Point>[] holes;
WB_Polygon polygon;
List<WB_Polygon> boundary;

void setup() {
  size(800, 800,OPENGL);
  smooth(8);
  render=new WB_Render2D(this);
  shell= new ArrayList<WB_Point>();
  for (int i=0; i<20; i++) {
    shell.add(gf.createPointFromPolar(50+75*(i%2+1), TWO_PI/20.0*i));
  }
  holes= new ArrayList[2];
  ArrayList<WB_Point> hole=new ArrayList<WB_Point>();
  for (int i=0; i<10; i++) {
    hole.add(gf.createPointFromPolar(20*(i%2+1), -TWO_PI/10.*i).addSelf(-50, 0, 0));
  } 
  holes[0]=hole;
  hole=new ArrayList<WB_Point>();
  for (int i=0; i<10; i++) {
    hole.add(gf.createPointFromPolar(20*(i%2+1), PI-TWO_PI/10.*i).addSelf(50, 0, 0));
  } 
  holes[1]=hole;
  polygon=gf.createPolygonWithHoles(shell, holes);
  boundary=gf.createBoundaryPolygons(polygon);
 
}

void draw() {
  background(55);
  translate(width/2, height/2);
  scale(1, -1);
  fill(255, 0, 0);
  noStroke();
  render.drawPolygon2D(polygon);
  noFill();
  stroke(0);
  strokeWeight(2);
  for (WB_Polygon poly : boundary) {
    render.drawPolygonEdges2D(poly);
  }
}