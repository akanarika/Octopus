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

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;

import com.vividsolutions.jts.algorithm.ConvexHull;
import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryCollection;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.geom.Polygon;
import com.vividsolutions.jts.triangulate.DelaunayTriangulationBuilder;
import com.vividsolutions.jts.triangulate.quadedge.QuadEdgeSubdivision;

import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import javolution.util.FastTable;
import wblut.external.Delaunay.WB_Delaunay;
import wblut.external.ProGAL.CTetrahedron;
import wblut.external.ProGAL.CVertex;
import wblut.external.ProGAL.DelaunayComplex;

/**
 *
 */
public class WB_Voronoi {

	/**
	 *
	 */
	private static WB_GeometryFactory geometryfactory = new WB_GeometryFactory();
	final static WB_Map2D XY = geometryfactory.createEmbeddedPlane();

	/**
	 *
	 *
	 * @param points
	 * @param n
	 * @param aabb
	 * @param precision
	 * @return
	 */
	public static List<WB_VoronoiCell3D> getVoronoi3D(final WB_Coord[] points, final int n, final WB_AABB aabb,
			final double precision) {
		final WB_Delaunay triangulation = WB_Delaunay.getTriangulation3D(points, precision);
		final int nv = Math.min(n, points.length);
		if (nv <= 4) {
			return getVoronoi3DBruteForce(points, nv, aabb);
		}
		final List<WB_VoronoiCell3D> result = new FastTable<WB_VoronoiCell3D>();
		for (int i = 0; i < nv; i++) {
			final int[] tetras = triangulation.Vertices[i];
			final List<WB_Point> hullpoints = new ArrayList<WB_Point>();
			for (int t = 0; t < tetras.length; t++) {
				hullpoints.add(triangulation.circumcenters[tetras[t]]);
			}
			final WB_VoronoiCell3D vor = new WB_VoronoiCell3D(hullpoints, geometryfactory.createPoint(points[i]), i);
			if (vor.cell != null) {
				vor.constrain(aabb);
			}
			if (vor.cell != null) {
				result.add(vor);
			}
		}
		return result;
	}

	/**
	 *
	 *
	 * @param points
	 * @param n
	 * @param aabb
	 * @param precision
	 * @return
	 */
	public static List<WB_VoronoiCell3D> getVoronoi3D(final List<? extends WB_Coord> points, final int n,
			final WB_AABB aabb, final double precision) {
		final WB_Delaunay triangulation = WB_Delaunay.getTriangulation3D(points, precision);
		final int nv = Math.min(n, points.size());
		if (nv <= 4) {
			return getVoronoi3DBruteForce(points, nv, aabb);
		}
		final List<WB_VoronoiCell3D> result = new FastTable<WB_VoronoiCell3D>();
		for (int i = 0; i < nv; i++) {
			final int[] tetras = triangulation.Vertices[i];
			final List<WB_Point> hullpoints = new ArrayList<WB_Point>();
			for (int t = 0; t < tetras.length; t++) {
				hullpoints.add(triangulation.circumcenters[tetras[t]]);
			}
			hullpoints.add(new WB_Point(points.get(i)));
			final WB_VoronoiCell3D vor = new WB_VoronoiCell3D(hullpoints, geometryfactory.createPoint(points.get(i)),
					i);
			if (vor.cell != null) {
				vor.constrain(aabb);
			}
			if (vor.cell != null) {
				result.add(vor);
			}
		}
		return result;
	}

	/**
	 *
	 *
	 * @param points
	 * @param aabb
	 * @param precision
	 * @return
	 */
	public static List<WB_VoronoiCell3D> getVoronoi3D(final List<? extends WB_Coord> points, final WB_AABB aabb,
			final double precision) {
		return getVoronoi3D(points, points.size(), aabb, precision);
	}

	/**
	 *
	 *
	 * @param points
	 * @param aabb
	 * @param precision
	 * @return
	 */
	public static List<WB_VoronoiCell3D> getVoronoi3D(final WB_Coord[] points, final WB_AABB aabb,
			final double precision) {
		return getVoronoi3D(points, points.length, aabb, precision);
	}

	/**
	 *
	 *
	 * @param points
	 * @param aabb
	 * @return
	 */
	public static List<WB_VoronoiCell3D> getVoronoi3D(final WB_Coord[] points, final WB_AABB aabb) {
		return getVoronoi3D(points, points.length, aabb);
	}

	/**
	 *
	 *
	 * @param points
	 * @param nv
	 * @param aabb
	 * @return
	 */
	public static List<WB_VoronoiCell3D> getVoronoi3D(final WB_Coord[] points, int nv, final WB_AABB aabb) {
		nv = Math.min(nv, points.length);
		if (nv <= 4) {
			return getVoronoi3DBruteForce(points, nv, aabb);
		}
		final int n = points.length;
		final List<wblut.external.ProGAL.Point> tmppoints = new ArrayList<wblut.external.ProGAL.Point>(n);
		final WB_KDTreeInteger<WB_Coord> tree = new WB_KDTreeInteger<WB_Coord>();
		for (int i = 0; i < n; i++) {
			tmppoints.add(new wblut.external.ProGAL.Point(points[i].xd(), points[i].yd(), points[i].zd()));
			tree.add(points[i], i);
		}
		final DelaunayComplex dc = new DelaunayComplex(tmppoints);
		final List<CVertex> vertices = dc.getVertices();
		final List<WB_VoronoiCell3D> result = new FastTable<WB_VoronoiCell3D>();
		for (int i = 0; i < nv; i++) {
			final CVertex v = vertices.get(i);
			final Set<CTetrahedron> vertexhull = dc.getVertexHull(v);
			final List<WB_Point> hullpoints = new ArrayList<WB_Point>();
			for (final CTetrahedron tetra : vertexhull) {
				// if (!tetra.containsBigPoint()) {
				hullpoints.add(toPoint(tetra.circumcenter()));
				// }
			}
			final List<WB_Point> finalpoints = new FastTable<WB_Point>();
			for (int j = 0; j < hullpoints.size(); j++) {
				finalpoints.add(geometryfactory.createPoint(hullpoints.get(j)));
			}
			final int index = tree.getNearestNeighbor(toPoint(v)).value;
			final WB_VoronoiCell3D vor = new WB_VoronoiCell3D(finalpoints, geometryfactory.createPoint(points[index]),
					index);
			if (aabb != null) {
				vor.constrain(aabb);
			}
			if (vor.cell != null) {
				result.add(vor);
			}
		}
		return result;
	}

	/**
	 *
	 *
	 * @param points
	 * @return
	 */
	public static int[][] getVoronoi3DNeighbors(final WB_Coord[] points) {
		final int nv = points.length;
		if (nv == 2) {
			return new int[][] { { 1 }, { 0 } };
		} else if (nv == 3) {
			return new int[][] { { 1, 2 }, { 0, 2 }, { 0, 1 } };
		} else if (nv == 4) {
			return new int[][] { { 1, 2, 3 }, { 0, 2, 3 }, { 0, 1, 3 }, { 0, 1, 2 } };
		}
		final List<wblut.external.ProGAL.Point> tmppoints = new ArrayList<wblut.external.ProGAL.Point>(nv);
		final WB_KDTreeInteger<WB_Coord> tree = new WB_KDTreeInteger<WB_Coord>();
		for (int i = 0; i < nv; i++) {
			tmppoints.add(new wblut.external.ProGAL.Point(points[i].xd(), points[i].yd(), points[i].zd()));
			tree.add(points[i], i);
		}
		final DelaunayComplex dc = new DelaunayComplex(tmppoints);
		final List<CVertex> vertices = dc.getVertices();
		final int[][] ns = new int[nv][];
		for (int i = 0; i < nv; i++) {
			final CVertex v = vertices.get(i);
			final Set<CTetrahedron> vertexhull = dc.getVertexHull(v);
			final TIntSet neighbors = new TIntHashSet();
			for (final CTetrahedron tetra : vertexhull) {
				for (int j = 0; j < 4; j++) {
					if (!tetra.getPoint(j).isBigpoint()) {
						neighbors.add(tree.getNearestNeighbor(toPoint(tetra.getPoint(j))).value);
					}
				}
			}
			ns[i] = neighbors.toArray();
		}
		return ns;
	}

	/**
	 *
	 *
	 * @param points
	 * @return
	 */
	public static int[][] getVoronoi3DNeighbors(final List<? extends WB_Coord> points) {
		final int nv = points.size();
		if (nv == 2) {
			return new int[][] { { 1 }, { 0 } };
		} else if (nv == 3) {
			return new int[][] { { 1, 2 }, { 0, 2 }, { 0, 1 } };
		} else if (nv == 4) {
			return new int[][] { { 1, 2, 3 }, { 0, 2, 3 }, { 0, 1, 3 }, { 0, 1, 2 } };
		}
		final List<wblut.external.ProGAL.Point> tmppoints = new ArrayList<wblut.external.ProGAL.Point>(nv);
		final WB_KDTreeInteger<WB_Coord> tree = new WB_KDTreeInteger<WB_Coord>();
		WB_Coord p;
		for (int i = 0; i < nv; i++) {
			p = points.get(i);
			tmppoints.add(new wblut.external.ProGAL.Point(p.xd(), p.yd(), p.zd()));
			tree.add(p, i);
		}
		final DelaunayComplex dc = new DelaunayComplex(tmppoints);
		final List<CVertex> vertices = dc.getVertices();
		final int[][] ns = new int[nv][];
		for (int i = 0; i < nv; i++) {
			final CVertex v = vertices.get(i);
			final Set<CTetrahedron> vertexhull = dc.getVertexHull(v);
			final TIntSet neighbors = new TIntHashSet();
			for (final CTetrahedron tetra : vertexhull) {
				for (int j = 0; j < 4; j++) {
					if (!tetra.getPoint(j).isBigpoint()) {
						neighbors.add(tree.getNearestNeighbor(toPoint(tetra.getPoint(j))).value);
					}
				}
			}
			ns[i] = neighbors.toArray();
		}
		return ns;
	}

	/**
	 *
	 *
	 * @param points
	 * @param aabb
	 * @return
	 */
	public static List<WB_VoronoiCell3D> getVoronoi3D(final List<? extends WB_Coord> points, final WB_AABB aabb) {
		return getVoronoi3D(points, points.size(), aabb);
	}

	/**
	 *
	 *
	 * @param points
	 * @param nv
	 * @param aabb
	 * @return
	 */
	public static List<WB_VoronoiCell3D> getVoronoi3D(final List<? extends WB_Coord> points, int nv,
			final WB_AABB aabb) {
		nv = Math.min(nv, points.size());
		if (nv <= 4) {
			return getVoronoi3DBruteForce(points, nv, aabb);
		}
		final int n = points.size();
		final List<wblut.external.ProGAL.Point> tmppoints = new ArrayList<wblut.external.ProGAL.Point>(n);
		final WB_KDTreeInteger<WB_Coord> tree = new WB_KDTreeInteger<WB_Coord>();
		int i = 0;
		for (final WB_Coord p : points) {
			tmppoints.add(new wblut.external.ProGAL.Point(p.xd(), p.yd(), p.zd()));
			tree.add(p, i++);
		}
		final DelaunayComplex dc = new DelaunayComplex(tmppoints);
		final List<CVertex> vertices = dc.getVertices();
		final List<WB_VoronoiCell3D> result = new FastTable<WB_VoronoiCell3D>();
		for (i = 0; i < nv; i++) {
			final CVertex v = vertices.get(i);
			final Set<CTetrahedron> vertexhull = dc.getVertexHull(v);
			v.getAdjacentTriangles();
			final List<WB_Point> hullpoints = new ArrayList<WB_Point>();
			for (final CTetrahedron tetra : vertexhull) {
				// if (!tetra.containsBigPoint()) {
				hullpoints.add(toPoint(tetra.circumcenter()));
				// }
			}
			final List<WB_Point> finalpoints = new FastTable<WB_Point>();
			for (int j = 0; j < hullpoints.size(); j++) {
				finalpoints.add(geometryfactory.createPoint(hullpoints.get(j)));
			}
			final int index = tree.getNearestNeighbor(toPoint(v)).value;
			final WB_VoronoiCell3D vor = new WB_VoronoiCell3D(finalpoints,
					geometryfactory.createPoint(points.get(index)), index);
			if (vor.cell != null) {
				vor.constrain(aabb);
			}
			if (vor.cell != null) {
				result.add(vor);
			}
		}
		return result;
	}

	/**
	 *
	 *
	 * @param points
	 * @param nv
	 * @param aabb
	 * @return
	 */
	public static List<WB_VoronoiCell3D> getVoronoi3DBruteForce(final List<? extends WB_Coord> points, int nv,
			final WB_AABB aabb) {
		nv = Math.min(nv, points.size());
		final int n = points.size();
		final List<WB_VoronoiCell3D> result = new FastTable<WB_VoronoiCell3D>();
		for (int i = 0; i < nv; i++) {
			final ArrayList<WB_Plane> cutPlanes = new ArrayList<WB_Plane>();
			final WB_Point O = new WB_Point();
			WB_Plane P;
			final WB_Mesh cell = geometryfactory.createMesh(aabb);
			for (int j = 0; j < n; j++) {
				if (j != i) {
					final WB_Vector N = new WB_Vector(points.get(i));
					N.subSelf(points.get(j));
					N.normalizeSelf();
					O.set(points.get(i)); // plane origin=point halfway
					// between point i and point j
					O.addSelf(points.get(j));
					O.mulSelf(0.5);
					P = new WB_Plane(O, N);
					cutPlanes.add(P);
				}
			}
			boolean unique;
			final ArrayList<WB_Plane> cleaned = new ArrayList<WB_Plane>();
			for (int j = 0; j < cutPlanes.size(); j++) {
				P = cutPlanes.get(j);
				unique = true;
				for (int k = 0; k < j; k++) {
					final WB_Plane Pj = cutPlanes.get(j);
					if (WB_GeometryOp3D.isEqual(P, Pj)) {
						unique = false;
						break;
					}
				}
				if (unique) {
					cleaned.add(P);
				}
			}
			final WB_VoronoiCell3D vor = new WB_VoronoiCell3D(cell, geometryfactory.createPoint(points.get(i)), i);
			if (vor.cell != null) {
				vor.constrain(cutPlanes);
			}
			if (vor.cell != null) {
				result.add(vor);
			}
			result.add(vor);
		}
		return result;
	}

	/**
	 *
	 *
	 * @param points
	 * @param aabb
	 * @return
	 */
	public static List<WB_VoronoiCell3D> getVoronoi3DBruteForce(final List<? extends WB_Coord> points,
			final WB_AABB aabb) {

		return getVoronoi3DBruteForce(points, points.size(), aabb);

	}

	/**
	 *
	 *
	 * @param points
	 * @param nv
	 * @param aabb
	 * @return
	 */
	public static List<WB_VoronoiCell3D> getVoronoi3DBruteForce(final WB_Coord[] points, int nv, final WB_AABB aabb) {
		nv = Math.min(nv, points.length);
		final int n = points.length;
		final List<WB_VoronoiCell3D> result = new FastTable<WB_VoronoiCell3D>();
		for (int i = 0; i < nv; i++) {
			final ArrayList<WB_Plane> cutPlanes = new ArrayList<WB_Plane>();
			final WB_Point O = new WB_Point();
			WB_Plane P;
			final WB_Mesh cell = geometryfactory.createMesh(aabb);
			for (int j = 0; j < n; j++) {
				if (j != i) {
					final WB_Vector N = new WB_Vector(points[i]);
					N.subSelf(points[j]);
					N.normalizeSelf();
					O.set(points[i]); // plane origin=point halfway
					// between point i and point j
					O.addSelf(points[j]);
					O.mulSelf(0.5);
					P = new WB_Plane(O, N);
					cutPlanes.add(P);
				}
			}
			boolean unique;
			final ArrayList<WB_Plane> cleaned = new ArrayList<WB_Plane>();
			for (int j = 0; j < cutPlanes.size(); j++) {
				P = cutPlanes.get(j);
				unique = true;
				for (int k = 0; k < j; k++) {
					final WB_Plane Pj = cutPlanes.get(j);
					if (WB_GeometryOp3D.isEqual(P, Pj)) {
						unique = false;
						break;
					}
				}
				if (unique) {
					cleaned.add(P);
				}
			}
			final WB_VoronoiCell3D vor = new WB_VoronoiCell3D(cell, geometryfactory.createPoint(points[i]), i);
			if (vor.cell != null) {
				vor.constrain(cutPlanes);
			}
			if (vor.cell != null) {
				result.add(vor);
			}
			result.add(vor);
		}
		return result;
	}

	/**
	 *
	 *
	 * @param points
	 * @param aabb
	 * @return
	 */
	public static List<WB_VoronoiCell3D> getVoronoi3DBruteForce(final WB_Coord[] points, final WB_AABB aabb) {
		return getVoronoi3DBruteForce(points, points.length, aabb);
	}

	/**
	 *
	 *
	 * @param v
	 * @return
	 */
	private static WB_Point toPoint(final wblut.external.ProGAL.Point v) {
		return geometryfactory.createPoint(v.x(), v.y(), v.z());
	}

	/**
	 *
	 *
	 * @param points
	 * @param context
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getVoronoi2D(final WB_Coord[] points, final WB_Map2D context) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			coords.add(toCoordinate(points[i], i, context));
		}
		return voronoi2D(coords, context);
	}

	/**
	 *
	 *
	 * @param points
	 * @param context
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getVoronoi2D(final Collection<? extends WB_Coord> points,
			final WB_Map2D context) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			id++;
		}
		return voronoi2D(coords, context);
	}

	/**
	 *
	 *
	 * @param points
	 * @param d
	 * @param context
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getVoronoi2D(final WB_Coord[] points, final double d, final WB_Map2D context) {
		return getVoronoi2D(points, d, 2, context);
	}

	/**
	 *
	 *
	 * @param points
	 * @param d
	 * @param context
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getVoronoi2D(final Collection<? extends WB_Coord> points, final double d,
			final WB_Map2D context) {
		return getVoronoi2D(points, d, 2, context);
	}

	/**
	 *
	 *
	 * @param points
	 * @param d
	 * @param c
	 * @param context
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getVoronoi2D(final WB_Coord[] points, final double d, final int c,
			final WB_Map2D context) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			coords.add(toCoordinate(points[i], i, context));
		}
		return voronoi2D(coords, d, c, context);
	}

	/**
	 *
	 *
	 * @param points
	 * @param d
	 * @param c
	 * @param context
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getVoronoi2D(final Collection<? extends WB_Coord> points, final double d,
			final int c, final WB_Map2D context) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			id++;
		}
		return voronoi2D(coords, d, c, context);
	}

	/**
	 *
	 *
	 * @param points
	 * @param d
	 * @param c
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getVoronoi2D(final Collection<? extends WB_Coord> points, final double d,
			final int c) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, XY));
			id++;
		}
		return voronoi2D(coords, d, c, XY);
	}

	/**
	 *
	 *
	 * @param points
	 * @param d
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getVoronoi2D(final Collection<? extends WB_Coord> points, final double d) {
		return getVoronoi2D(points, d, 2, XY);
	}

	/**
	 *
	 *
	 * @param points
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getVoronoi2D(final Collection<? extends WB_Coord> points) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, XY));
			id++;
		}
		return voronoi2D(coords, XY);
	}

	/**
	 *
	 *
	 * @param points
	 * @param d
	 * @param c
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getVoronoi2D(final WB_Coord[] points, final double d, final int c) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			coords.add(toCoordinate(points[i], i, XY));
		}
		return voronoi2D(coords, d, c, XY);
	}

	/**
	 *
	 *
	 * @param points
	 * @param d
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getVoronoi2D(final WB_Coord[] points, final double d) {
		return getVoronoi2D(points, d, 2, XY);
	}

	/**
	 *
	 *
	 * @param points
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getVoronoi2D(final WB_Coord[] points) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			coords.add(toCoordinate(points[i], i, XY));
		}
		return voronoi2D(coords, XY);
	}

	/**
	 *
	 *
	 * @param points
	 * @param context
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final WB_Coord[] points, final WB_Map2D context) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			coords.add(toCoordinate(points[i], i, context));
		}
		return clippedVoronoi2D(coords, context);
	}

	/**
	 *
	 *
	 * @param points
	 * @param boundary
	 * @param context
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final WB_Coord[] points, final WB_Coord[] boundary,
			final WB_Map2D context) {
		int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			coords.add(toCoordinate(points[i], i, context));
		}
		n = boundary.length;
		final ArrayList<Coordinate> bdcoords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			bdcoords.add(toCoordinate(boundary[i], i, context));
		}
		if (!bdcoords.get(0).equals(bdcoords.get(n - 1))) {
			bdcoords.add(bdcoords.get(0));
		}
		return clippedVoronoi2D(coords, bdcoords, context);
	}

	/**
	 *
	 *
	 * @param points
	 * @param boundary
	 * @param context
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final WB_Coord[] points, final WB_Polygon boundary,
			final WB_Map2D context) {
		int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			coords.add(toCoordinate(points[i], i, context));
		}

		return clippedVoronoi2D(coords, boundary, context);
	}

	/**
	 *
	 *
	 * @param points
	 * @param context
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final Collection<? extends WB_Coord> points,
			final WB_Map2D context) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			id++;
		}
		return clippedVoronoi2D(coords, context);
	}

	/**
	 *
	 *
	 * @param points
	 * @param boundary
	 * @param context
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final Collection<? extends WB_Coord> points,
			final Collection<? extends WB_Coord> boundary, final WB_Map2D context) {
		int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			id++;
		}
		n = boundary.size();
		final ArrayList<Coordinate> bdcoords = new ArrayList<Coordinate>(n);
		id = 0;
		for (final WB_Coord p : boundary) {
			bdcoords.add(toCoordinate(p, id, context));
			id++;
		}
		if (!bdcoords.get(0).equals(bdcoords.get(n - 1))) {
			bdcoords.add(bdcoords.get(0));
		}
		return clippedVoronoi2D(coords, bdcoords, context);
	}

	/**
	 *
	 *
	 * @param points
	 * @param boundary
	 * @param context
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final Collection<? extends WB_Coord> points,
			final WB_Polygon boundary, final WB_Map2D context) {
		int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			id++;
		}

		return clippedVoronoi2D(coords, boundary, context);
	}

	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final Collection<? extends WB_Coord> points,
			final List<WB_Polygon> boundary, final WB_Map2D context) {
		int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			id++;
		}

		return clippedVoronoi2D(coords, boundary, context);
	}

	/**
	 *
	 *
	 * @param points
	 * @param d
	 * @param context
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final WB_Coord[] points, final double d,
			final WB_Map2D context) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			coords.add(toCoordinate(points[i], i, context));
		}
		return clippedVoronoi2D(coords, d, 2, context);
	}

	/**
	 *
	 *
	 * @param points
	 * @param boundary
	 * @param d
	 * @param context
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final WB_Coord[] points, final WB_Coord[] boundary,
			final double d, final WB_Map2D context) {
		int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			coords.add(toCoordinate(points[i], i, context));
		}
		n = boundary.length;
		final ArrayList<Coordinate> bdcoords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			bdcoords.add(toCoordinate(boundary[i], i, context));
		}
		if (!bdcoords.get(0).equals(bdcoords.get(n - 1))) {
			bdcoords.add(bdcoords.get(0));
		}
		return clippedVoronoi2D(coords, bdcoords, d, 2, context);
	}

	/**
	 *
	 *
	 * @param points
	 * @param boundary
	 * @param d
	 * @param context
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final WB_Coord[] points, final WB_Polygon boundary,
			final double d, final WB_Map2D context) {
		int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			coords.add(toCoordinate(points[i], i, context));
		}

		return clippedVoronoi2D(coords, boundary, d, 2, context);
	}

	/**
	 *
	 *
	 * @param points
	 * @param boundary
	 * @param d
	 * @param context
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final Collection<? extends WB_Coord> points,
			final Collection<? extends WB_Coord> boundary, final double d, final WB_Map2D context) {
		int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			id++;
		}
		n = boundary.size();
		final ArrayList<Coordinate> bdcoords = new ArrayList<Coordinate>(n);
		id = 0;
		for (final WB_Coord p : boundary) {
			bdcoords.add(toCoordinate(p, id, context));
			id++;
		}
		if (!bdcoords.get(0).equals(bdcoords.get(n - 1))) {
			bdcoords.add(bdcoords.get(0));
		}
		return clippedVoronoi2D(coords, bdcoords, d, 2, context);
	}

	/**
	 *
	 *
	 * @param points
	 * @param boundary
	 * @param d
	 * @param context
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final Collection<? extends WB_Coord> points,
			final WB_Polygon boundary, final double d, final WB_Map2D context) {
		int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			id++;
		}

		return clippedVoronoi2D(coords, boundary, d, 2, context);
	}

	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final Collection<? extends WB_Coord> points,
			final List<WB_Polygon> boundary, final double d, final WB_Map2D context) {
		int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			id++;
		}

		return clippedVoronoi2D(coords, boundary, d, 2, context);
	}

	/**
	 *
	 *
	 * @param points
	 * @param d
	 * @param context
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final Collection<? extends WB_Coord> points,
			final double d, final WB_Map2D context) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			id++;
		}
		return clippedVoronoi2D(coords, d, 2, context);
	}

	/**
	 *
	 *
	 * @param points
	 * @param d
	 * @param c
	 * @param context
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final WB_Coord[] points, final double d, final int c,
			final WB_Map2D context) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			coords.add(toCoordinate(points[i], i, context));
		}
		return clippedVoronoi2D(coords, d, c, context);
	}

	/**
	 *
	 *
	 * @param points
	 * @param boundary
	 * @param d
	 * @param c
	 * @param context
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final WB_Coord[] points, final WB_Coord[] boundary,
			final double d, final int c, final WB_Map2D context) {
		int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			coords.add(toCoordinate(points[i], i, context));
		}
		n = boundary.length;
		final ArrayList<Coordinate> bdcoords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			bdcoords.add(toCoordinate(boundary[i], i, context));
		}
		if (!bdcoords.get(0).equals(bdcoords.get(n - 1))) {
			bdcoords.add(bdcoords.get(0));
		}
		return clippedVoronoi2D(coords, bdcoords, d, c, context);
	}

	/**
	 *
	 *
	 * @param points
	 * @param boundary
	 * @param d
	 * @param c
	 * @param context
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final WB_Coord[] points, final WB_Polygon boundary,
			final double d, final int c, final WB_Map2D context) {
		int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			coords.add(toCoordinate(points[i], i, context));
		}

		return clippedVoronoi2D(coords, boundary, d, c, context);
	}

	/**
	 *
	 *
	 * @param points
	 * @param d
	 * @param c
	 * @param context
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final Collection<? extends WB_Coord> points,
			final double d, final int c, final WB_Map2D context) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			id++;
		}
		return clippedVoronoi2D(coords, d, c, context);
	}

	/**
	 *
	 *
	 * @param points
	 * @param boundary
	 * @param d
	 * @param c
	 * @param context
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final Collection<? extends WB_Coord> points,
			final Collection<? extends WB_Coord> boundary, final double d, final int c, final WB_Map2D context) {
		int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			id++;
		}
		n = boundary.size();
		final ArrayList<Coordinate> bdcoords = new ArrayList<Coordinate>(n);
		id = 0;
		for (final WB_Coord p : boundary) {
			bdcoords.add(toCoordinate(p, id, context));
			id++;
		}
		if (!bdcoords.get(0).equals(bdcoords.get(n - 1))) {
			bdcoords.add(bdcoords.get(0));
		}
		return clippedVoronoi2D(coords, bdcoords, d, c, context);
	}

	/**
	 *
	 *
	 * @param points
	 * @param boundary
	 * @param d
	 * @param c
	 * @param context
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final Collection<? extends WB_Coord> points,
			final WB_Polygon boundary, final double d, final int c, final WB_Map2D context) {
		int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			id++;
		}

		return clippedVoronoi2D(coords, boundary, d, c, context);
	}

	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final Collection<? extends WB_Coord> points,
			final List<WB_Polygon> boundary, final double d, final int c, final WB_Map2D context) {
		int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, context));
			id++;
		}

		return clippedVoronoi2D(coords, boundary, d, c, context);
	}

	/**
	 *
	 *
	 * @param points
	 * @param boundary
	 * @param d
	 * @param c
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final Collection<? extends WB_Coord> points,
			final Collection<? extends WB_Coord> boundary, final double d, final int c) {
		int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, XY));
			id++;
		}
		n = boundary.size();
		final ArrayList<Coordinate> bdcoords = new ArrayList<Coordinate>(n);
		id = 0;
		for (final WB_Coord p : boundary) {
			bdcoords.add(toCoordinate(p, id, XY));
			id++;
		}
		if (!bdcoords.get(0).equals(bdcoords.get(n - 1))) {
			bdcoords.add(bdcoords.get(0));
		}
		return clippedVoronoi2D(coords, bdcoords, d, c, XY);
	}

	/**
	 *
	 *
	 * @param points
	 * @param boundary
	 * @param d
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final Collection<? extends WB_Coord> points,
			final Collection<? extends WB_Coord> boundary, final double d) {
		int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, XY));
			id++;
		}
		n = boundary.size();
		final ArrayList<Coordinate> bdcoords = new ArrayList<Coordinate>(n);
		id = 0;
		for (final WB_Coord p : boundary) {
			bdcoords.add(toCoordinate(p, id, XY));
			id++;
		}
		if (!bdcoords.get(0).equals(bdcoords.get(n - 1))) {
			bdcoords.add(bdcoords.get(0));
		}
		return clippedVoronoi2D(coords, bdcoords, d, 2, XY);
	}

	/**
	 *
	 *
	 * @param points
	 * @param boundary
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final Collection<? extends WB_Coord> points,
			final Collection<? extends WB_Coord> boundary) {
		int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, XY));
			id++;
		}
		n = boundary.size();
		final ArrayList<Coordinate> bdcoords = new ArrayList<Coordinate>(n);
		id = 0;
		for (final WB_Coord p : boundary) {
			bdcoords.add(toCoordinate(p, id, XY));
			id++;
		}
		if (!bdcoords.get(0).equals(bdcoords.get(n - 1))) {
			bdcoords.add(bdcoords.get(0));
		}
		return clippedVoronoi2D(coords, bdcoords, XY);
	}

	/**
	 *
	 *
	 * @param points
	 * @param d
	 * @param c
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final Collection<? extends WB_Coord> points,
			final double d, final int c) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, XY));
			id++;
		}
		return clippedVoronoi2D(coords, d, c, XY);
	}

	/**
	 *
	 *
	 * @param points
	 * @param d
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final Collection<? extends WB_Coord> points,
			final double d) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, XY));
			id++;
		}
		return clippedVoronoi2D(coords, d, 2, XY);
	}

	/**
	 *
	 *
	 * @param points
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final Collection<? extends WB_Coord> points) {
		final int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, XY));
			id++;
		}
		return clippedVoronoi2D(coords, XY);
	}

	/**
	 *
	 *
	 * @param points
	 * @param d
	 * @param c
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final WB_Coord[] points, final double d, final int c) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			coords.add(toCoordinate(points[i], i, XY));
		}
		return clippedVoronoi2D(coords, d, c, XY);
	}

	/**
	 *
	 *
	 * @param points
	 * @param d
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final WB_Coord[] points, final double d) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			coords.add(toCoordinate(points[i], i, XY));
		}
		return clippedVoronoi2D(coords, d, 2, XY);
	}

	/**
	 *
	 *
	 * @param points
	 * @param boundary
	 * @param d
	 * @param c
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final WB_Coord[] points, final WB_Coord[] boundary,
			final double d, final int c) {
		int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			coords.add(toCoordinate(points[i], i, XY));
		}
		n = boundary.length;
		final ArrayList<Coordinate> bdcoords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			bdcoords.add(toCoordinate(boundary[i], i, XY));
		}
		if (!bdcoords.get(0).equals(bdcoords.get(n - 1))) {
			bdcoords.add(bdcoords.get(0));
		}
		return clippedVoronoi2D(coords, bdcoords, d, c, XY);
	}

	/**
	 *
	 *
	 * @param points
	 * @param boundary
	 * @param d
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final WB_Coord[] points, final WB_Coord[] boundary,
			final double d) {
		int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			coords.add(toCoordinate(points[i], i, XY));
		}
		n = boundary.length;
		final ArrayList<Coordinate> bdcoords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			bdcoords.add(toCoordinate(boundary[i], i, XY));
		}
		if (!bdcoords.get(0).equals(bdcoords.get(n - 1))) {
			bdcoords.add(bdcoords.get(0));
		}
		return clippedVoronoi2D(coords, bdcoords, d, 2, XY);
	}

	/**
	 *
	 *
	 * @param points
	 * @param boundary
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final WB_Coord[] points, final WB_Coord[] boundary) {
		int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			coords.add(toCoordinate(points[i], i, XY));
		}
		n = boundary.length;
		final ArrayList<Coordinate> bdcoords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			bdcoords.add(toCoordinate(boundary[i], i, XY));
		}
		if (!bdcoords.get(0).equals(bdcoords.get(n - 1))) {
			bdcoords.add(bdcoords.get(0));
		}
		return clippedVoronoi2D(coords, bdcoords, XY);
	}

	/**
	 *
	 *
	 * @param points
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final WB_Coord[] points) {
		final int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			coords.add(toCoordinate(points[i], i, XY));
		}
		return clippedVoronoi2D(coords, XY);
	}

	/**
	 *
	 *
	 * @param points
	 * @param boundary
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final WB_Coord[] points, final WB_Polygon boundary) {
		int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			coords.add(toCoordinate(points[i], i, XY));
		}

		return clippedVoronoi2D(coords, boundary, XY);
	}

	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final WB_Coord[] points, final List<WB_Polygon> boundary) {
		int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			coords.add(toCoordinate(points[i], i, XY));
		}

		return clippedVoronoi2D(coords, boundary, XY);
	}

	/**
	 *
	 *
	 * @param points
	 * @param boundary
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final Collection<? extends WB_Coord> points,
			final WB_Polygon boundary) {
		int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, XY));
			id++;
		}

		return clippedVoronoi2D(coords, boundary, XY);
	}

	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final Collection<? extends WB_Coord> points,
			final List<WB_Polygon> boundary) {
		int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, XY));
			id++;
		}

		return clippedVoronoi2D(coords, boundary, XY);
	}

	/**
	 *
	 *
	 * @param points
	 * @param boundary
	 * @param d
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final WB_Coord[] points, final WB_Polygon boundary,
			final double d) {
		int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			coords.add(toCoordinate(points[i], i, XY));
		}

		return clippedVoronoi2D(coords, boundary, d, 2, XY);
	}

	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final WB_Coord[] points, final List<WB_Polygon> boundary,
			final double d) {
		int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			coords.add(toCoordinate(points[i], i, XY));
		}

		return clippedVoronoi2D(coords, boundary, d, 2, XY);
	}

	/**
	 *
	 *
	 * @param points
	 * @param boundary
	 * @param d
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final Collection<? extends WB_Coord> points,
			final WB_Polygon boundary, final double d) {
		int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, XY));
			id++;
		}

		return clippedVoronoi2D(coords, boundary, d, 2, XY);
	}

	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final Collection<? extends WB_Coord> points,
			final List<WB_Polygon> boundary, final double d) {
		int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, XY));
			id++;
		}

		return clippedVoronoi2D(coords, boundary, d, 2, XY);
	}

	/**
	 *
	 *
	 * @param points
	 * @param boundary
	 * @param d
	 * @param c
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final WB_Coord[] points, final WB_Polygon boundary,
			final double d, final int c) {
		int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			coords.add(toCoordinate(points[i], i, XY));
		}

		return clippedVoronoi2D(coords, boundary, d, c, XY);
	}

	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final WB_Coord[] points, final List<WB_Polygon> boundary,
			final double d, final int c) {
		int n = points.length;
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		for (int i = 0; i < n; i++) {
			coords.add(toCoordinate(points[i], i, XY));
		}

		return clippedVoronoi2D(coords, boundary, d, c, XY);
	}

	/**
	 *
	 *
	 * @param points
	 * @param boundary
	 * @param d
	 * @param c
	 * @return
	 */
	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final Collection<? extends WB_Coord> points,
			final WB_Polygon boundary, final double d, final int c) {
		int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, XY));
			id++;
		}

		return clippedVoronoi2D(coords, boundary, d, c, XY);
	}

	public static List<WB_VoronoiCell2D> getClippedVoronoi2D(final Collection<? extends WB_Coord> points,
			final List<WB_Polygon> boundary, final double d, final int c) {
		int n = points.size();
		final ArrayList<Coordinate> coords = new ArrayList<Coordinate>(n);
		int id = 0;
		for (final WB_Coord p : points) {
			coords.add(toCoordinate(p, id, XY));
			id++;
		}

		return clippedVoronoi2D(coords, boundary, d, c, XY);
	}

	/**
	 *
	 *
	 * @param coords
	 * @param context
	 * @return
	 */
	private static List<WB_VoronoiCell2D> voronoi2D(final ArrayList<Coordinate> coords, final WB_Map2D context) {
		final DelaunayTriangulationBuilder dtb = new DelaunayTriangulationBuilder();
		dtb.setSites(coords);
		final QuadEdgeSubdivision qes = dtb.getSubdivision();
		final GeometryCollection polys = (GeometryCollection) qes.getVoronoiDiagram(new GeometryFactory());
		final int npolys = polys.getNumGeometries();
		final List<WB_VoronoiCell2D> result = new FastTable<WB_VoronoiCell2D>();
		for (int i = 0; i < npolys; i++) {
			final Polygon poly = (Polygon) polys.getGeometryN(i);
			final Coordinate[] polycoord = poly.getCoordinates();
			final List<WB_Coord> polypoints = new FastTable<WB_Coord>();
			for (final Coordinate element : polycoord) {
				polypoints.add(toPoint(element.x, element.y, context));
			}
			final Point centroid = poly.getCentroid();
			final WB_Point pc = centroid == null ? null : toPoint(centroid.getX(), centroid.getY(), context);
			final int index = (int) ((Coordinate) poly.getUserData()).z;
			final double area = poly.getArea();
			result.add(
					new WB_VoronoiCell2D(polypoints, index, geometryfactory.createPoint(coords.get(index)), area, pc));
		}
		return result;
	}

	/**
	 *
	 *
	 * @param coords
	 * @param d
	 * @param c
	 * @param context
	 * @return
	 */
	private static List<WB_VoronoiCell2D> voronoi2D(final ArrayList<Coordinate> coords, final double d, final int c,
			final WB_Map2D context) {
		final DelaunayTriangulationBuilder dtb = new DelaunayTriangulationBuilder();
		dtb.setSites(coords);
		final QuadEdgeSubdivision qes = dtb.getSubdivision();
		final GeometryCollection polys = (GeometryCollection) qes.getVoronoiDiagram(new GeometryFactory());
		final int npolys = polys.getNumGeometries();
		final List<WB_VoronoiCell2D> result = new FastTable<WB_VoronoiCell2D>();
		for (int i = 0; i < npolys; i++) {
			Geometry poly = polys.getGeometryN(i);
			poly = poly.buffer(-d, c);
			poly = polys.getGeometryN(0);
			final Coordinate[] polycoord = poly.getCoordinates();
			final List<WB_Point> polypoints = new FastTable<WB_Point>();
			;
			for (final Coordinate element : polycoord) {
				polypoints.add(toPoint(element.x, element.y, context));
			}
			final Point centroid = poly.getCentroid();
			final WB_Point pc = centroid == null ? null : toPoint(centroid.getX(), centroid.getY(), context);
			final int index = (int) ((Coordinate) poly.getUserData()).z;
			final double area = poly.getArea();
			result.add(
					new WB_VoronoiCell2D(polypoints, index, geometryfactory.createPoint(coords.get(index)), area, pc));
		}
		return result;
	}

	/**
	 *
	 *
	 * @param coords
	 * @param context
	 * @return
	 */
	private static List<WB_VoronoiCell2D> clippedVoronoi2D(final ArrayList<Coordinate> coords, final WB_Map2D context) {
		final DelaunayTriangulationBuilder dtb = new DelaunayTriangulationBuilder();
		dtb.setSites(coords);
		final QuadEdgeSubdivision qes = dtb.getSubdivision();
		final GeometryCollection polys = (GeometryCollection) qes.getVoronoiDiagram(new GeometryFactory());
		final int npolys = polys.getNumGeometries();
		final List<WB_VoronoiCell2D> result = new FastTable<WB_VoronoiCell2D>();
		Coordinate[] coordsArray = new Coordinate[coords.size()];
		coordsArray = coords.toArray(coordsArray);
		final ConvexHull ch = new ConvexHull(coordsArray, new GeometryFactory());
		final Geometry hull = ch.getConvexHull();
		for (int i = 0; i < npolys; i++) {
			Polygon poly = (Polygon) polys.getGeometryN(i);
			final Geometry intersect = poly.intersection(hull.getGeometryN(0));
			final double cellindex = ((Coordinate) poly.getUserData()).z;
			for (int j = 0; j < intersect.getNumGeometries(); j++) {
				if (intersect.getGeometryN(j).getGeometryType().equals("Polygon")
						&& !intersect.getGeometryN(j).isEmpty()) {
					poly = (Polygon) intersect.getGeometryN(j);
					final Coordinate[] polycoord = poly.getCoordinates();
					final List<WB_Point> polypoints = new FastTable<WB_Point>();
					for (final Coordinate element : polycoord) {
						polypoints.add(toPoint(element.x, element.y, context));
					}
					final Point centroid = poly.getCentroid();
					final WB_Point pc = centroid == null ? null : toPoint(centroid.getX(), centroid.getY(), context);
					final int index = (int) cellindex;
					final double area = poly.getArea();
					result.add(new WB_VoronoiCell2D(polypoints, index, geometryfactory.createPoint(coords.get(index)),
							area, pc));
				}
			}
		}
		return result;
	}

	/**
	 *
	 *
	 * @param coords
	 * @param d
	 * @param c
	 * @param context
	 * @return
	 */
	private static List<WB_VoronoiCell2D> clippedVoronoi2D(final ArrayList<Coordinate> coords, final double d,
			final int c, final WB_Map2D context) {
		final DelaunayTriangulationBuilder dtb = new DelaunayTriangulationBuilder();
		dtb.setSites(coords);
		final QuadEdgeSubdivision qes = dtb.getSubdivision();
		final GeometryCollection polys = (GeometryCollection) qes.getVoronoiDiagram(new GeometryFactory());
		final int npolys = polys.getNumGeometries();
		final List<WB_VoronoiCell2D> result = new FastTable<WB_VoronoiCell2D>();
		Coordinate[] coordsArray = new Coordinate[coords.size()];
		coordsArray = coords.toArray(coordsArray);
		final ConvexHull ch = new ConvexHull(coordsArray, new GeometryFactory());
		final Geometry hull = ch.getConvexHull();
		for (int i = 0; i < npolys; i++) {
			Polygon poly = (Polygon) polys.getGeometryN(i);
			Geometry intersect = poly.intersection(hull.getGeometryN(0));
			intersect = intersect.buffer(-d, c);
			final double cellindex = ((Coordinate) poly.getUserData()).z;
			for (int j = 0; j < intersect.getNumGeometries(); j++) {
				if (intersect.getGeometryN(j).getGeometryType().equals("Polygon")
						&& !intersect.getGeometryN(j).isEmpty()) {
					poly = (Polygon) intersect.getGeometryN(j);
					final Coordinate[] polycoord = poly.getCoordinates();
					final List<WB_Point> polypoints = new FastTable<WB_Point>();
					for (final Coordinate element : polycoord) {
						polypoints.add(toPoint(element.x, element.y, context));
					}
					final Point centroid = poly.getCentroid();
					final WB_Point pc = centroid == null ? null : toPoint(centroid.getX(), centroid.getY(), context);
					final int index = (int) cellindex;
					final double area = poly.getArea();
					result.add(new WB_VoronoiCell2D(polypoints, index, geometryfactory.createPoint(coords.get(index)),
							area, pc));
				}
			}
		}
		return result;
	}

	/**
	 *
	 *
	 * @param coords
	 * @param bdcoords
	 * @param context
	 * @return
	 */
	private static List<WB_VoronoiCell2D> clippedVoronoi2D(final ArrayList<Coordinate> coords,
			final ArrayList<Coordinate> bdcoords, final WB_Map2D context) {
		final DelaunayTriangulationBuilder dtb = new DelaunayTriangulationBuilder();
		dtb.setSites(coords);
		final QuadEdgeSubdivision qes = dtb.getSubdivision();
		final GeometryCollection polys = (GeometryCollection) qes.getVoronoiDiagram(new GeometryFactory());
		final int npolys = polys.getNumGeometries();
		final List<WB_VoronoiCell2D> result = new FastTable<WB_VoronoiCell2D>();
		Coordinate[] bdcoordsArray = new Coordinate[bdcoords.size()];
		bdcoordsArray = bdcoords.toArray(bdcoordsArray);
		final Polygon hull = new GeometryFactory().createPolygon(bdcoordsArray);
		for (int i = 0; i < npolys; i++) {
			Polygon poly = (Polygon) polys.getGeometryN(i);
			final Geometry intersect = poly.intersection(hull);
			final double cellindex = ((Coordinate) poly.getUserData()).z;
			for (int j = 0; j < intersect.getNumGeometries(); j++) {
				if (intersect.getGeometryN(j).getGeometryType().equals("Polygon")
						&& !intersect.getGeometryN(j).isEmpty()) {
					poly = (Polygon) intersect.getGeometryN(j);

					final Coordinate[] polycoord = poly.getCoordinates();
					final List<WB_Point> polypoints = new FastTable<WB_Point>();
					for (final Coordinate element : polycoord) {
						polypoints.add(toPoint(element.x, element.y, context));
					}
					final Point centroid = poly.getCentroid();
					final WB_Point pc = centroid == null ? null : toPoint(centroid.getX(), centroid.getY(), context);
					final int index = (int) cellindex;
					final double area = poly.getArea();
					result.add(new WB_VoronoiCell2D(polypoints, index, geometryfactory.createPoint(coords.get(index)),
							area, pc));
				}
			}
		}
		return result;
	}

	/**
	 *
	 *
	 * @param coords
	 * @param constraint
	 * @param context
	 * @return
	 */
	private static List<WB_VoronoiCell2D> clippedVoronoi2D(final ArrayList<Coordinate> coords,
			final WB_Polygon constraint, final WB_Map2D context) {
		final DelaunayTriangulationBuilder dtb = new DelaunayTriangulationBuilder();
		dtb.setSites(coords);
		final QuadEdgeSubdivision qes = dtb.getSubdivision();
		final GeometryCollection polys = (GeometryCollection) qes.getVoronoiDiagram(new GeometryFactory());
		final int npolys = polys.getNumGeometries();
		final List<WB_VoronoiCell2D> result = new FastTable<WB_VoronoiCell2D>();

		final Polygon hull = geometryfactory.toJTSPolygon2D(constraint);
		for (int i = 0; i < npolys; i++) {
			Polygon poly = (Polygon) polys.getGeometryN(i);
			final Geometry intersect = poly.intersection(hull);
			final double cellindex = ((Coordinate) poly.getUserData()).z;
			for (int j = 0; j < intersect.getNumGeometries(); j++) {
				if (intersect.getGeometryN(j).getGeometryType().equals("Polygon")
						&& !intersect.getGeometryN(j).isEmpty()) {
					poly = (Polygon) intersect.getGeometryN(j);

					final Coordinate[] polycoord = poly.getCoordinates();
					final List<WB_Point> polypoints = new FastTable<WB_Point>();
					for (final Coordinate element : polycoord) {
						polypoints.add(toPoint(element.x, element.y, context));
					}
					final Point centroid = poly.getCentroid();
					final WB_Point pc = centroid == null ? null : toPoint(centroid.getX(), centroid.getY(), context);
					final int index = (int) cellindex;
					final double area = poly.getArea();
					result.add(new WB_VoronoiCell2D(polypoints, index, geometryfactory.createPoint(coords.get(index)),
							area, pc));
				}
			}
		}
		return result;
	}

	/**
	 *
	 *
	 * @param coords
	 * @param constraint
	 * @param context
	 * @return
	 */
	private static List<WB_VoronoiCell2D> clippedVoronoi2D(final ArrayList<Coordinate> coords,
			final List<WB_Polygon> constraint, final WB_Map2D context) {
		final DelaunayTriangulationBuilder dtb = new DelaunayTriangulationBuilder();
		dtb.setSites(coords);
		final QuadEdgeSubdivision qes = dtb.getSubdivision();
		final GeometryCollection polys = (GeometryCollection) qes.getVoronoiDiagram(new GeometryFactory());
		final int npolys = polys.getNumGeometries();
		final List<WB_VoronoiCell2D> result = new FastTable<WB_VoronoiCell2D>();

		final Geometry hull = geometryfactory.toJTSMultiPolygon2D(constraint);
		for (int i = 0; i < npolys; i++) {
			Polygon poly = (Polygon) polys.getGeometryN(i);
			final Geometry intersect = poly.intersection(hull);
			final double cellindex = ((Coordinate) poly.getUserData()).z;
			for (int j = 0; j < intersect.getNumGeometries(); j++) {
				if (intersect.getGeometryN(j).getGeometryType().equals("Polygon")
						&& !intersect.getGeometryN(j).isEmpty()) {
					poly = (Polygon) intersect.getGeometryN(j);

					final Coordinate[] polycoord = poly.getCoordinates();
					final List<WB_Point> polypoints = new FastTable<WB_Point>();
					for (final Coordinate element : polycoord) {
						polypoints.add(toPoint(element.x, element.y, context));
					}
					final Point centroid = poly.getCentroid();
					final WB_Point pc = centroid == null ? null : toPoint(centroid.getX(), centroid.getY(), context);
					final int index = (int) cellindex;
					final double area = poly.getArea();
					result.add(new WB_VoronoiCell2D(polypoints, index, geometryfactory.createPoint(coords.get(index)),
							area, pc));
				}
			}
		}
		return result;
	}

	/**
	 *
	 *
	 * @param coords
	 * @param bdcoords
	 * @param d
	 * @param c
	 * @param context
	 * @return
	 */
	private static List<WB_VoronoiCell2D> clippedVoronoi2D(final ArrayList<Coordinate> coords,
			final ArrayList<Coordinate> bdcoords, final double d, final int c, final WB_Map2D context) {
		final DelaunayTriangulationBuilder dtb = new DelaunayTriangulationBuilder();
		dtb.setSites(coords);
		final QuadEdgeSubdivision qes = dtb.getSubdivision();
		final GeometryCollection polys = (GeometryCollection) qes.getVoronoiDiagram(new GeometryFactory());
		final int npolys = polys.getNumGeometries();
		final List<WB_VoronoiCell2D> result = new FastTable<WB_VoronoiCell2D>();
		Coordinate[] bdcoordsArray = new Coordinate[bdcoords.size()];
		bdcoordsArray = bdcoords.toArray(bdcoordsArray);
		final Polygon hull = new GeometryFactory().createPolygon(bdcoordsArray);
		for (int i = 0; i < npolys; i++) {
			Polygon poly = (Polygon) polys.getGeometryN(i);
			Geometry intersect = poly.intersection(hull);
			intersect = intersect.buffer(-d, c);
			final double cellindex = ((Coordinate) poly.getUserData()).z;
			for (int j = 0; j < intersect.getNumGeometries(); j++) {
				if (intersect.getGeometryN(j).getGeometryType().equals("Polygon")
						&& !intersect.getGeometryN(j).isEmpty()) {
					poly = (Polygon) intersect.getGeometryN(j);

					final Coordinate[] polycoord = poly.getCoordinates();
					final List<WB_Point> polypoints = new FastTable<WB_Point>();
					;
					for (final Coordinate element : polycoord) {
						polypoints.add(toPoint(element.x, element.y, context));
					}
					final Point centroid = poly.getCentroid();
					final WB_Point pc = centroid == null ? null : toPoint(centroid.getX(), centroid.getY(), context);
					final int index = (int) cellindex;
					final double area = poly.getArea();
					result.add(new WB_VoronoiCell2D(polypoints, index, geometryfactory.createPoint(coords.get(index)),
							area, pc));
				}
			}
		}
		return result;
	}

	/**
	 *
	 *
	 * @param coords
	 * @param constraint
	 * @param d
	 * @param c
	 * @param context
	 * @return
	 */
	private static List<WB_VoronoiCell2D> clippedVoronoi2D(final ArrayList<Coordinate> coords,
			final WB_Polygon constraint, final double d, final int c, final WB_Map2D context) {
		final DelaunayTriangulationBuilder dtb = new DelaunayTriangulationBuilder();
		dtb.setSites(coords);
		final QuadEdgeSubdivision qes = dtb.getSubdivision();
		final GeometryCollection polys = (GeometryCollection) qes.getVoronoiDiagram(new GeometryFactory());
		final int npolys = polys.getNumGeometries();
		final List<WB_VoronoiCell2D> result = new FastTable<WB_VoronoiCell2D>();
		final Polygon hull = geometryfactory.toJTSPolygon2D(constraint);
		for (int i = 0; i < npolys; i++) {
			Polygon poly = (Polygon) polys.getGeometryN(i);
			Geometry intersect = poly.intersection(hull);
			intersect = intersect.buffer(-d, c);
			final double cellindex = ((Coordinate) poly.getUserData()).z;
			for (int j = 0; j < intersect.getNumGeometries(); j++) {
				if (intersect.getGeometryN(j).getGeometryType().equals("Polygon")
						&& !intersect.getGeometryN(j).isEmpty()) {
					poly = (Polygon) intersect.getGeometryN(j);

					final Coordinate[] polycoord = poly.getCoordinates();
					final List<WB_Point> polypoints = new FastTable<WB_Point>();
					;
					for (final Coordinate element : polycoord) {
						polypoints.add(toPoint(element.x, element.y, context));
					}
					final Point centroid = poly.getCentroid();
					final WB_Point pc = centroid == null ? null : toPoint(centroid.getX(), centroid.getY(), context);
					final int index = (int) cellindex;
					final double area = poly.getArea();
					result.add(new WB_VoronoiCell2D(polypoints, index, geometryfactory.createPoint(coords.get(index)),
							area, pc));
				}
			}
		}
		return result;
	}

	private static List<WB_VoronoiCell2D> clippedVoronoi2D(final ArrayList<Coordinate> coords,
			final List<WB_Polygon> constraint, final double d, final int c, final WB_Map2D context) {
		final DelaunayTriangulationBuilder dtb = new DelaunayTriangulationBuilder();
		dtb.setSites(coords);
		final QuadEdgeSubdivision qes = dtb.getSubdivision();
		final GeometryCollection polys = (GeometryCollection) qes.getVoronoiDiagram(new GeometryFactory());
		final int npolys = polys.getNumGeometries();
		final List<WB_VoronoiCell2D> result = new FastTable<WB_VoronoiCell2D>();
		final Geometry hull = geometryfactory.toJTSMultiPolygon2D(constraint);
		for (int i = 0; i < npolys; i++) {
			Polygon poly = (Polygon) polys.getGeometryN(i);
			Geometry intersect = poly.intersection(hull);
			intersect = intersect.buffer(-d, c);
			final double cellindex = ((Coordinate) poly.getUserData()).z;
			for (int j = 0; j < intersect.getNumGeometries(); j++) {
				if (intersect.getGeometryN(j).getGeometryType().equals("Polygon")
						&& !intersect.getGeometryN(j).isEmpty()) {
					poly = (Polygon) intersect.getGeometryN(j);

					final Coordinate[] polycoord = poly.getCoordinates();
					final List<WB_Point> polypoints = new FastTable<WB_Point>();
					;
					for (final Coordinate element : polycoord) {
						polypoints.add(toPoint(element.x, element.y, context));
					}
					final Point centroid = poly.getCentroid();
					final WB_Point pc = centroid == null ? null : toPoint(centroid.getX(), centroid.getY(), context);
					final int index = (int) cellindex;
					final double area = poly.getArea();
					result.add(new WB_VoronoiCell2D(polypoints, index, geometryfactory.createPoint(coords.get(index)),
							area, pc));
				}
			}
		}
		return result;
	}

	/**
	 *
	 *
	 * @param p
	 * @param i
	 * @param context
	 * @return
	 */
	private static Coordinate toCoordinate(final WB_Coord p, final int i, final WB_Map2D context) {
		final WB_Point tmp = geometryfactory.createPoint();
		context.mapPoint3D(p, tmp);
		final Coordinate c = new Coordinate(tmp.xd(), tmp.yd(), i);
		return c;
	}

	/**
	 *
	 *
	 * @param x
	 * @param y
	 * @param context
	 * @return
	 */
	private static WB_Point toPoint(final double x, final double y, final WB_Map2D context) {
		final WB_Point tmp = geometryfactory.createPoint();
		context.unmapPoint3D(x, y, 0, tmp);
		return tmp;
	}
}
