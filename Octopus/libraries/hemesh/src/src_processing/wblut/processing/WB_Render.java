/*
 * This file is part of HE_Mesh, a library for creating and manipulating meshes.
 * It is dedicated to the public domain. To the extent possible under law,
 * I , Frederik Vanhoutte, have waived all copyright and related or neighboring
 * rights.
 * 
 * This work is published from Belgium. (http://creativecommons.org/publicdomain/zero/1.0/)
 * 
 */
package wblut.processing;

import processing.core.PApplet;
import processing.opengl.PGraphics3D;

/**
 * 
 */
public class WB_Render extends WB_Render3D {
    
    /**
     * 
     *
     * @param home 
     */
    public WB_Render(final PApplet home) {
	super(home);
    }

    /**
     * 
     *
     * @param home 
     */
    public WB_Render(final PGraphics3D home) {
	super(home);
    }
}
