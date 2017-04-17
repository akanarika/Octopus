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
 * Interface for implementing metric operations on 3D coordinates.
 *
 * None of the operators change the calling object.
 *
 * @author Frederik Vanhoutte
 *
 */
public interface WB_CoordinateMetric3D extends WB_CoordinateMetric2D {

	/**
	 *
	 *
	 * @return
	 * @deprecated Use {@link #getLength()} instead
	 */
	@Deprecated
	public double getLength3D();

	/**
	 *
	 *
	 * @return
	 * @deprecated Use {@link #getSqLength()} instead
	 */
	@Deprecated
	public double getSqLength3D();

	/**
	 * 
	 * @param p
	 * @return
	 */
	@Deprecated
	public double getDistance3D(final WB_Coord p);

	/**
	 * 
	 * @param p
	 * @return
	 */
	@Deprecated
	public double getSqDistance3D(final WB_Coord p);

	/**
	 *
	 *
	 * @return
	 * @deprecated Use {@link #getOrthoNormal()} instead
	 */
	@Deprecated
	public WB_Coord getOrthoNormal3D();

	/**
	 * Length of 3D coordinate.
	 *
	 * @return
	 */
	public double getLength();

	/**
	 * Square length of 3D coordinate.
	 *
	 * @return
	 */
	public double getSqLength();

	/**
	 * 3D distance to coordinate
	 *
	 * @param p
	 * @return
	 */
	public double getDistance(final WB_Coord p);

	/**
	 * Square 3D distance to coordinate.
	 *
	 * @param p
	 * @return
	 */
	public double getSqDistance(final WB_Coord p);

	/**
	 * Get vector perpendicular and CCW to this one.
	 *
	 * @return
	 */
	public WB_Coord getOrthoNormal();

	/**
	 * Is this a degenerate vector?
	 *
	 * @return
	 */
	@Override
	public boolean isZero();

	/**
	 * Is this point collinear with two other points?
	 *
	 * @param p
	 * @param q
	 * @return
	 */
	public boolean isCollinear(WB_Coord p, WB_Coord q);

	/**
	 * Is this vector parallel with other vector?
	 *
	 * @param p
	 * @return
	 */
	public boolean isParallel(WB_Coord p);

	/**
	 * Is this vector, within a given tolerance, parallel with other vector?
	 *
	 * @param p
	 * @param tol
	 * @return
	 */
	public boolean isParallel(WB_Coord p, double tol);

	/**
	 * Is this normalized vector parallel with other normalized vector?
	 *
	 * @param p
	 * @return
	 */
	public boolean isParallelNorm(WB_Coord p);

	/**
	 * Is this normalized vector, within a given tolerance, parallel with other
	 * normalized vector?
	 *
	 * @param p
	 * @param tol
	 * @return
	 */
	public boolean isParallelNorm(WB_Coord p, double tol);

	/**
	 * Is this vector perpendicular to other vector?
	 *
	 * @param p
	 * @return
	 */
	public boolean isOrthogonal(WB_Coord p);

	/**
	 * Is this vector, within a given tolerance, perpendicular to other vector?
	 *
	 * @param p
	 * @param tol
	 * @return
	 */
	public boolean isOrthogonal(WB_Coord p, double tol);

	/**
	 * Is this normalized vector perpendicular to other normalized vector?
	 *
	 * @param p
	 * @return
	 */
	public boolean isOrthogonalNorm(WB_Coord p);

	/**
	 * Is this normalized vector, within a given tolerance, perpendicular to
	 * other normalized vector?
	 *
	 * @param p
	 * @param tol
	 * @return
	 */
	public boolean isOrthogonalNorm(WB_Coord p, double tol);

}
