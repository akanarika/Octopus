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

public interface WB_Curve {

	/**
	 *
	 *
	 * @param u
	 * @return
	 */
	public WB_Point curvePoint(double u);

	/**
	 *
	 *
	 * @param u
	 * @return
	 */
	public WB_Vector curveDirection(double u);

	/**
	 *
	 *
	 * @param u
	 * @return
	 */
	public WB_Vector curveDerivative(double u);

	/**
	 *
	 *
	 * @return
	 */
	public double getLowerU();

	/**
	 *
	 *
	 * @return
	 */
	public double getUpperU();
}
