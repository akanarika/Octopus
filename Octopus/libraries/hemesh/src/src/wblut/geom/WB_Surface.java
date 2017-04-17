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
 * Interface for parameterized surfaces.
 *
 * @author Frederik Vanhoutte, W:Blut
 */
public interface WB_Surface {
	/**
	 * Retrieve the point at values (u,v).
	 * 
	 * @param u
	 * 
	 * @param v
	 * 
	 * @return WB_Point
	 */
	public WB_Point surfacePoint(double u, double v);

	/**
	 * Get the lower end of the u parameter range.
	 * 
	 * @return u
	 */
	public double loweru();

	/**
	 * Get the upper end of the u parameter range.
	 * 
	 * @return u
	 */
	public double upperu();

	/**
	 * Get the lower end of the v parameter range.
	 * 
	 * @return v
	 */
	public double lowerv();

	/**
	 * Get the upper end of the v parameter range.
	 * 
	 * @return v
	 */
	public double upperv();
}
