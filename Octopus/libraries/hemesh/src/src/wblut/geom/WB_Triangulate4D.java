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

import java.util.List;

import wblut.external.Delaunay.WB_Delaunay;

/**
 *
 */
public class WB_Triangulate4D extends WB_Triangulate3D {

	/**
	 *
	 */
	public WB_Triangulate4D() {
	}

	/**
	 *
	 *
	 * @param points
	 * @param closest
	 * @return
	 */
	public static WB_Triangulation4D triangulate4D(final List<? extends WB_Coord> points, final double closest) {
		return new WB_Triangulation4D(WB_Delaunay.getTriangulation4D(points, closest).Tri);

	}

	/**
	 *
	 *
	 * @param points
	 * @param closest
	 * @param epsilon
	 * @return
	 */
	public static WB_Triangulation4D triangulate4D(final List<? extends WB_Coord> points, final double closest,
			final double epsilon) {
		return new WB_Triangulation4D(WB_Delaunay.getTriangulation4D(points, closest, epsilon).Tri);
	}

	/**
	 *
	 *
	 * @param points
	 * @param closest
	 * @return
	 */
	public static WB_Triangulation4D triangulate4D(final WB_Coord[] points, final double closest) {
		return new WB_Triangulation4D(WB_Delaunay.getTriangulation4D(points, closest).Tri);
	}

	/**
	 *
	 *
	 * @param points
	 * @param closest
	 * @param epsilon
	 * @return
	 */
	public static WB_Triangulation4D triangulate4D(final WB_Coord[] points, final double closest,
			final double epsilon) {
		return new WB_Triangulation4D(WB_Delaunay.getTriangulation4D(points, closest, epsilon).Tri);
	}

}
