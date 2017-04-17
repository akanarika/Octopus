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

public  interface HET_InfoVertex<E extends Object>{
	
	/**
	 * 
	 *
	 * @param vertex 
	 * @return 
	 */
	public E retrieve(final HE_Vertex vertex);
}