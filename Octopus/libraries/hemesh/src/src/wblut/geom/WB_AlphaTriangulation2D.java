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

import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import javolution.util.FastMap;
import javolution.util.FastTable;

/**
 * WB_AlphaTriangulation2D stores the results of
 * WB_Triangulate.alphaTriangulate2D: a 2D Delaunay triangulation with the
 * corresponding circumcircle radii.
 *
 */
public class WB_AlphaTriangulation2D {

	/**
	 *
	 */
	private int[] triangles;
	private double[] alpha;
	private FastTable<WB_Coord> points;

	/**
	 *
	 *
	 * @param tris
	 * @param points
	 */
	public WB_AlphaTriangulation2D(final int[] tris, final Collection<? extends WB_Coord> points) {
		triangles = Arrays.copyOf(tris, tris.length);
		this.points = new FastTable<WB_Coord>();
		this.points.addAll(points);
		setAlpha();
	}

	/**
	 * 
	 * @param tris
	 * @param points
	 */
	public WB_AlphaTriangulation2D(final int[][] tris, final Collection<? extends WB_Coord> points) {
		triangles = new int[tris.length * 4];
		for (int i = 0; i < tris.length; i++) {
			triangles[4 * i] = tris[i][0];
			triangles[4 * i + 1] = tris[i][1];
			triangles[4 * i + 2] = tris[i][2];
			triangles[4 * i + 3] = tris[i][3];
		}
		this.points = new FastTable<WB_Coord>();
		this.points.addAll(points);
		setAlpha();

	}

	/**
	 *
	 * @param tris
	 * @param points
	 */
	public WB_AlphaTriangulation2D(final int[] tris, final WB_Coord[] points) {
		triangles = Arrays.copyOf(tris, tris.length);
		this.points = new FastTable<WB_Coord>();
		for (WB_Coord p : points) {
			this.points.add(p);
		}
		setAlpha();
	}

	/**
	 *
	 * @param tris
	 * @param points
	 */
	public WB_AlphaTriangulation2D(final int[][] tris, final WB_Coord[] points) {
		triangles = new int[tris.length * 3];
		for (int i = 0; i < tris.length; i++) {
			triangles[3 * i] = tris[i][0];
			triangles[3 * i + 1] = tris[i][1];
			triangles[3 * i + 2] = tris[i][2];
		}
		this.points = new FastTable<WB_Coord>();
		for (WB_Coord p : points) {
			this.points.add(p);
		}
		setAlpha();

	}

	private void setAlpha() {

		alpha = new double[triangles.length / 3];
		int index = 0;
		for (int i = 0; i < triangles.length; i += 3) {
			alpha[index++] = WB_Predicates.circumradius2D(points.get(triangles[i]), points.get(triangles[i + 1]),
					points.get(triangles[i + 2]));

		}

	}

	/**
	 * Get the indices to the triangles vertices as a single array of int. 3
	 * indices per triangle.
	 *
	 * @return
	 */
	public int[] getTriangles() {
		return triangles;
	}

	/**
	 * Get the vertices of the triangulation as an unmodifiable List<WB_Coord>.
	 *
	 * @return
	 */
	public List<WB_Coord> getpoints() {
		return points.unmodifiable();
	}

	/**
	 * Get the circumradii of the triangle vertices as a single array of double.
	 * Original values are copied.
	 *
	 * @return
	 */
	public double[] getAlpha() {
		return Arrays.copyOf(alpha, alpha.length);
	}

	/**
	 * Get the indices to the alpha triangle vertices as a single array of int.
	 * 3 indices per triangle. Only triangles with a circumradius smaller or
	 * equal to a are returned.
	 *
	 * @param a
	 *            alpha value
	 *
	 * @return
	 */
	public int[] getAlphaTriangles(final double a) {
		int[] alphaTriangles = new int[triangles.length];
		int index = 0;
		for (int i = 0; i < triangles.length; i += 3) {
			if (alpha[i / 3] <= a) {
				alphaTriangles[index++] = triangles[i + 0];
				alphaTriangles[index++] = triangles[i + 1];
				alphaTriangles[index++] = triangles[i + 2];

			}
		}
		return Arrays.copyOf(alphaTriangles, index);
	}

	/*
	 * @param a alpha value
	 *
	 * @return
	 */
	public int[] getAlphaEdges(final double a) {

		FastMap<Key, Tuple> edges = new FastMap<Key, Tuple>();
		for (int i = 0; i < triangles.length; i += 4) {
			if (alpha[i / 4] <= a) {
				Key key = new Key(triangles[i], triangles[i + 1]);
				Tuple doubleId = new Tuple(triangles[i], triangles[i + 1]);
				if (edges.get(key) == null) {
					edges.put(key, doubleId);
				} else {
					edges.remove(key);
				}

				key = new Key(triangles[i + 1], triangles[i + 2]);
				doubleId = new Tuple(triangles[i + 1], triangles[i + 2]);
				if (edges.get(key) == null) {
					edges.put(key, doubleId);
				} else {
					edges.remove(key);
				}

				key = new Key(triangles[i + 2], triangles[i]);
				doubleId = new Tuple(triangles[i + 2], triangles[i]);
				if (edges.get(key) == null) {
					edges.put(key, doubleId);
				} else {
					edges.remove(key);
				}

			}

		}
		int[] alphaEdges = new int[2 * edges.size()];
		int index = 0;
		for (Tuple T : edges.values()) {
			alphaEdges[index++] = T.i;
			alphaEdges[index++] = T.j;
		}

		return Arrays.copyOf(alphaEdges, index);
	}

	private class Tuple {
		int i, j;

		Tuple(final int i, final int j) {
			this.i = i;
			this.j = j;
			;
		}

	}

	private class Key {
		int a, b;

		public Key(final int i, final int j) {
			a = Math.min(i, j);
			b = Math.max(i, j);

		}

		@Override
		public int hashCode() {
			return a << 5 ^ b;
		}

		@Override
		public boolean equals(final Object o) {
			if (!(o instanceof Key)) {
				return false;
			}

			Key k = (Key) o;

			return k.a == a && k.b == b;
		}

	}
}