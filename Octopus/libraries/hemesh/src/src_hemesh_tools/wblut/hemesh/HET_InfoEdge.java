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

public  interface HET_InfoEdge<E extends Object>{
	
	/**
	 * 
	 *
	 * @param edge 
	 * @return 
	 */
	public E retrieve(final HE_Halfedge edge);


}