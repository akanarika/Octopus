/*
 * This file is part of HE_Mesh, a library for creating and manipulating meshes.
 * It is dedicated to the public domain. To the extent possible under law,
 * I , Frederik Vanhoutte, have waived all copyright and related or neighboring
 * rights.
 * 
 * This work is published from Belgium. (http://creativecommons.org/publicdomain/zero/1.0/)
 * 
 */

package wblut.geom;


public class WB_IntersectionResult {


	public double t1 = Float.NEGATIVE_INFINITY;


	public double t2 = Float.NEGATIVE_INFINITY;;


	public boolean intersection = false;


	public double sqDist = Float.POSITIVE_INFINITY;


	public Object object;


	public int dimension = -1;
}
