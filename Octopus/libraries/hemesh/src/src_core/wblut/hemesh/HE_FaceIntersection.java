/*
 * This file is part of HE_Mesh, a library for creating and manipulating meshes.
 * It is dedicated to the public domain. To the extent possible under law,
 * I , Frederik Vanhoutte, have waived all copyright and related or neighboring
 * rights.
 * 
 * This work is published from Belgium. (http://creativecommons.org/publicdomain/zero/1.0/)
 * 
 */
package wblut.hemesh;

import wblut.geom.WB_Point;

/**
 * 
 */
public class HE_FaceIntersection {
    
    /**
     * 
     */
    public HE_Face face;
    
    /**
     * 
     */
    public WB_Point point;

    /**
     * 
     *
     * @param f 
     * @param p 
     */
    public HE_FaceIntersection(final HE_Face f, final WB_Point p) {
	face = f;
	point = p;
    }
}
