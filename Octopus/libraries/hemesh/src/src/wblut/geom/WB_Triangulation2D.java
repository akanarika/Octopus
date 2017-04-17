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

/**
 *
 */
public class WB_Triangulation2D {

	/**
	 * 
	 */
	private int[] _triangles;

	/**
	 * 
	 */
	private int[] _edges;

	/**
	 * 
	 */
	public WB_Triangulation2D() {
	}

	/**
	 * 
	 *
	 * @param T
	 * @param E
	 */
	public WB_Triangulation2D(final int[] T, final int[] E) {
		_triangles = T;
		_edges = E;
	}

	/**
	 * 
	 *
	 * @return
	 */
	public int[] getTriangles() {
		return _triangles;
	}

	/**
	 * 
	 *
	 * @return
	 */
	public int[] getEdges() {
		return _edges;
	}
}