 PShape createNewShape(PShape r)
    {
       //PImage tex = loadImage("HeartC.JPG");
      PShape s = createShape();
      s.beginShape(TRIANGLES);
      s.noStroke();
      //s.texture(tex);
      //s.textureMode(NORMAL);
      for (int i=0; i<r.getChildCount (); i++) {
        if (r.getChild(i).getVertexCount() ==3) {
          for (int j=0; j<r.getChild (i).getVertexCount(); j++) {
            PVector p = r.getChild(i).getVertex(j);
            PVector n = r.getChild(i).getNormal(j);
            float u = r.getChild(i).getTextureU(j);
            float v = r.getChild(i).getTextureV(j);
            s.normal(n.x, n.y, n.z);
            s.vertex(p.x, p.y, p.z, u, v);
          }
        }
      }
      s.endShape();
      s.beginShape(QUADS);
      s.noStroke();
      //s.texture(tex);
      //s.textureMode(NORMAL);
      for (int i=0; i<r.getChildCount (); i++) {
        if (r.getChild(i).getVertexCount() ==4) {
          for (int j=0; j<r.getChild (i).getVertexCount(); j++) {
            PVector p = r.getChild(i).getVertex(j);
            PVector n = r.getChild(i).getNormal(j);
            float u = r.getChild(i).getTextureU(j);
            float v = r.getChild(i).getTextureV(j);
            s.normal(n.x, n.y, n.z);
            s.vertex(p.x, p.y, p.z, u, v);
          }
        }
      }
      s.endShape();
      return s;
    }
    
    PShape createShapeTri(PShape r) {
      //PImage tex = loadImage("HeartC.JPG");
      PShape s = createShape();
      s.beginShape(TRIANGLES);
      s.noStroke();
      //s.texture(tex);
      //s.textureMode(NORMAL);
      for (int i=0; i<r.getChildCount (); i++) {
        if (r.getChild(i).getVertexCount() ==3) {
          for (int j=0; j<r.getChild (i).getVertexCount(); j++) {
            PVector p = r.getChild(i).getVertex(j);
            PVector n = r.getChild(i).getNormal(j);
            float u = r.getChild(i).getTextureU(j);
            float v = r.getChild(i).getTextureV(j);
            s.normal(n.x, n.y, n.z);
            s.vertex(p.x, p.y, p.z, u, v);
          }
        }
      }
      s.endShape();
      return s;
    }
    
    PShape createShapeQuad(PShape r) {
      //PImage tex = loadImage("HeartC.JPG");
      PShape s = createShape();
      s.beginShape(QUADS);
      s.noStroke();
      //s.texture(tex);
      //s.textureMode(NORMAL);
      for (int i=0; i<r.getChildCount (); i++) {
        if (r.getChild(i).getVertexCount() ==4) {
          for (int j=0; j<r.getChild (i).getVertexCount(); j++) {
            PVector p = r.getChild(i).getVertex(j);
            PVector n = r.getChild(i).getNormal(j);
            float u = r.getChild(i).getTextureU(j);
            float v = r.getChild(i).getTextureV(j);
            s.normal(n.x, n.y, n.z);
            s.vertex(p.x, p.y, p.z, u, v);
          }
        }
      }
      s.endShape();
      return s;
    }