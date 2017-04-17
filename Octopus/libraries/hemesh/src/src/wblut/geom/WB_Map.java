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
 * WB_Map is an interface for classes that transform between 3D coordinates
 * through some from of mapping.
 *
 *
 */
public interface WB_Map {

	/**
	 * Map 3D point.
	 *
	 * @param p            3D point
	 * @param result            object implementing the WB_MutableCoordinate interface to
	 *            receive the result;
	 */
	public void mapPoint3D(WB_Coord p, WB_MutableCoord result);

	/**
	 * Map 3D point.
	 *
	 * @param x
	 * @param y
	 * @param z
	 * @param result            object implementing the WB_MutableCoordinate interface to
	 *            receive the result;
	 */
	public void mapPoint3D(double x, double y, double z, WB_MutableCoord result);

	/**
	 * Unmap 3D point.
	 *
	 * @param p
	 * @param result
	 *            object implementing the WB_MutableCoordinate interface to
	 *            receive the result;
	 */
	public void unmapPoint3D(WB_Coord p, WB_MutableCoord result);

	/**
	 * Unmap 3D point.
	 *
	 * @param u
	 * @param v
	 * @param w
	 * @param result
	 *            object implementing the WB_MutableCoordinate interface to
	 *            receive the result;
	 */
	public void unmapPoint3D(double u, double v, double w, WB_MutableCoord result);

	/**
	 * Unmap 2D point.
	 *
	 * @param p
	 * @param result
	 *            object implementing the WB_MutableCoordinate interface to
	 *            receive the result;
	 */

	/**
	 * Map 3D vector.
	 *
	 * @param p
	 * @param result
	 *            object implementing the WB_MutableCoordinate interface to
	 *            receive the result;
	 */
	public void mapVector3D(WB_Coord p, WB_MutableCoord result);

	/**
	 * Map 3D vector.
	 *
	 * @param x
	 * @param y
	 * @param z
	 * @param result
	 *            object implementing the WB_MutableCoordinate interface to
	 *            receive the result;
	 */
	public void mapVector3D(double x, double y, double z, WB_MutableCoord result);

	/**
	 * Unmap 3D vector.
	 *
	 * @param p
	 * @param result
	 *            object implementing the WB_MutableCoordinate interface to
	 *            receive the result;
	 */
	public void unmapVector3D(WB_Coord p, WB_MutableCoord result);

	/**
	 * Unmap 3D vector.
	 *
	 * @param u
	 * @param v
	 * @param w
	 * @param result
	 *            object implementing the WB_MutableCoordinate interface to
	 *            receive the result;
	 */
	public void unmapVector3D(double u, double v, double w, WB_MutableCoord result);

}
