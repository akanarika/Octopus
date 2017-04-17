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
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import javolution.util.FastTable;
import wblut.external.Delaunay.WB_Delaunay;
import wblut.external.ProGAL.CEdge;
import wblut.external.ProGAL.CTetrahedron;
import wblut.external.ProGAL.CTriangle;
import wblut.external.ProGAL.CVertex;
import wblut.external.ProGAL.DelaunayComplex;
import wblut.external.ProGAL.Point;

/**
 *
 */
public class WB_Triangulate3D extends WB_Triangulate2D {
	/**
	 *
	 */
	private static WB_GeometryFactory geometryfactory = new WB_GeometryFactory();

	/**
	 *
	 */
	public WB_Triangulate3D() {
		super();
	}

	/**
	 *
	 *
	 * @param points
	 * @param closest
	 * @return
	 */
	public static WB_Triangulation3D triangulate3D(final WB_Coord[] points, final double closest) {
		final WB_Triangulation3D result = new WB_Triangulation3D(WB_Delaunay.getTriangulation3D(points, closest).Tri);
		return result;
	}

	/**
	 *
	 *
	 * @param points
	 * @param closest
	 * @return
	 */
	public static WB_Triangulation3D triangulate3D(final Collection<? extends WB_Coord> points, final double closest) {
		final WB_Triangulation3D result = new WB_Triangulation3D(WB_Delaunay.getTriangulation3D(points, closest).Tri);
		return result;
	}

	/**
	 *
	 *
	 * @param points
	 * @return
	 */
	public static WB_Triangulation3D triangulate3D(final WB_Coord[] points) {
		// WB_Predicates predicates = new WB_Predicates();
		final int n = points.length;
		final List<Point> tmppoints = new ArrayList<Point>(n);
		final WB_KDTree<WB_Coord, Integer> tree = new WB_KDTree<WB_Coord, Integer>();
		for (int i = 0; i < n; i++) {
			tmppoints.add(new Point(points[i].xd(), points[i].yd(), points[i].zd()));
			tree.add(points[i], i);
		}
		final DelaunayComplex dc = new DelaunayComplex(tmppoints);
		final List<CTetrahedron> tetras = dc.getTetrahedra();
		final List<CTriangle> tris = dc.getTriangles();
		final List<CEdge> edges = dc.getEdges();
		int nt = tetras.size();
		List<int[]> tmpresult = new ArrayList<int[]>();
		for (int i = 0; i < nt; i++) {
			final int[] tmp = new int[4];
			final CTetrahedron tetra = tetras.get(i);
			int index = tree.getNearestNeighbor(convert(tetra.getPoint(0))).value;
			tmp[0] = index;
			index = tree.getNearestNeighbor(convert(tetra.getPoint(1))).value;
			tmp[1] = index;
			index = tree.getNearestNeighbor(convert(tetra.getPoint(2))).value;
			tmp[2] = index;
			index = tree.getNearestNeighbor(convert(tetra.getPoint(3))).value;
			tmp[3] = index;
			/*
			 * double o = predicates.orientTetra(points[tmp[0]].coords(),
			 * points[tmp[1]].coords(), points[tmp[2]].coords(),
			 * points[tmp[3]].coords()); if (o != 0) {
			 */
			tmpresult.add(tmp);
			/*
			 * }
			 */
		}
		final int[] tetra = new int[4 * tmpresult.size()];
		for (int i = 0; i < tmpresult.size(); i++) {
			for (int j = 0; j < 4; j++) {
				tetra[i * 4 + j] = tmpresult.get(i)[j];
			}
		}
		nt = tris.size();
		tmpresult = new ArrayList<int[]>();
		for (int i = 0; i < nt; i++) {
			final int[] tmp = new int[3];
			final CTriangle tri = tris.get(i);
			int index = tree.getNearestNeighbor(convert(tri.getPoint(0))).value;
			tmp[0] = index;
			index = tree.getNearestNeighbor(convert(tri.getPoint(1))).value;
			tmp[1] = index;
			index = tree.getNearestNeighbor(convert(tri.getPoint(2))).value;
			tmp[2] = index;
			/*
			 * double o = predicates.orientTetra(points[tmp[0]].coords(),
			 * points[tmp[1]].coords(), points[tmp[2]].coords(),
			 * points[tmp[3]].coords()); if (o != 0) {
			 */
			tmpresult.add(tmp);
			/*
			 * }
			 */
		}
		final int[] tri = new int[3 * tmpresult.size()];
		for (int i = 0; i < tmpresult.size(); i++) {
			for (int j = 0; j < 3; j++) {
				tri[3 * i + j] = tmpresult.get(i)[j];
			}
		}
		nt = edges.size();
		tmpresult = new ArrayList<int[]>();
		for (int i = 0; i < nt; i++) {
			final int[] tmp = new int[3];
			final CEdge edge = edges.get(i);
			int index = tree.getNearestNeighbor(convert(edge.getPoint(0))).value;
			tmp[0] = index;
			index = tree.getNearestNeighbor(convert(edge.getPoint(1))).value;
			tmp[1] = index;
			/*
			 * double o = predicates.orientTetra(points[tmp[0]].coords(),
			 * points[tmp[1]].coords(), points[tmp[2]].coords(),
			 * points[tmp[3]].coords()); if (o != 0) {
			 */
			tmpresult.add(tmp);
			/*
			 * }
			 */
		}
		final int[] edge = new int[2 * tmpresult.size()];
		for (int i = 0; i < tmpresult.size(); i++) {
			for (int j = 0; j < 2; j++) {
				edge[2 * i + j] = tmpresult.get(i)[j];
			}
		}
		final List<WB_Coord> pts = new FastTable<WB_Coord>();
		for (final WB_Coord p : points) {
			pts.add(p);
		}
		final WB_Triangulation3D result = new WB_Triangulation3D(tetra, tri, edge);
		return result;
	}

	/**
	 *
	 *
	 * @param points
	 * @return
	 */
	public static WB_Triangulation3D triangulate3D(final Collection<? extends WB_Coord> points) {
		final int n = points.size();
		final List<Point> tmppoints = new ArrayList<Point>(n);
		final WB_KDTree<WB_Coord, Integer> tree = new WB_KDTree<WB_Coord, Integer>();
		int i = 0;
		for (final WB_Coord p : points) {
			tmppoints.add(new Point(p.xd(), p.yd(), p.zd()));
			tree.add(p, i);
			i++;
		}
		final DelaunayComplex dc = new DelaunayComplex(tmppoints);
		final List<CTetrahedron> tetras = dc.getTetrahedra();
		final List<CTriangle> tris = dc.getTriangles();
		final List<CEdge> edges = dc.getEdges();
		int nt = tetras.size();
		List<int[]> tmpresult = new ArrayList<int[]>();
		for (i = 0; i < nt; i++) {
			final int[] tmp = new int[4];
			final CTetrahedron tetra = tetras.get(i);
			int index = tree.getNearestNeighbor(convert(tetra.getPoint(0))).value;
			tmp[0] = index;
			index = tree.getNearestNeighbor(convert(tetra.getPoint(1))).value;
			tmp[1] = index;
			index = tree.getNearestNeighbor(convert(tetra.getPoint(2))).value;
			tmp[2] = index;
			index = tree.getNearestNeighbor(convert(tetra.getPoint(3))).value;
			tmp[3] = index;
			/*
			 * double o = predicates.orientTetra(points[tmp[0]].coords(),
			 * points[tmp[1]].coords(), points[tmp[2]].coords(),
			 * points[tmp[3]].coords()); if (o != 0) {
			 */
			tmpresult.add(tmp);
			/*
			 * }
			 */
		}
		final int[] tetra = new int[4 * tmpresult.size()];
		for (i = 0; i < tmpresult.size(); i += 4) {
			for (int j = 0; j < 4; j++) {
				tetra[4 * i + j] = tmpresult.get(i)[j];
			}
		}
		nt = tris.size();
		tmpresult = new ArrayList<int[]>();
		for (i = 0; i < nt; i++) {
			final int[] tmp = new int[3];
			final CTriangle tri = tris.get(i);
			int index = tree.getNearestNeighbor(convert(tri.getPoint(0))).value;
			tmp[0] = index;
			index = tree.getNearestNeighbor(convert(tri.getPoint(1))).value;
			tmp[1] = index;
			index = tree.getNearestNeighbor(convert(tri.getPoint(2))).value;
			tmp[2] = index;
			/*
			 * double o = predicates.orientTetra(points[tmp[0]].coords(),
			 * points[tmp[1]].coords(), points[tmp[2]].coords(),
			 * points[tmp[3]].coords()); if (o != 0) {
			 */
			tmpresult.add(tmp);
			/*
			 * }
			 */
		}
		final int[] tri = new int[3 * tmpresult.size()];
		for (i = 0; i < tmpresult.size(); i++) {
			for (int j = 0; j < 3; j++) {
				tri[3 * i + j] = tmpresult.get(i)[j];
			}
		}
		nt = edges.size();
		tmpresult = new ArrayList<int[]>();
		for (i = 0; i < nt; i++) {
			final int[] tmp = new int[3];
			final CEdge edge = edges.get(i);
			int index = tree.getNearestNeighbor(convert(edge.getPoint(0))).value;
			tmp[0] = index;
			index = tree.getNearestNeighbor(convert(edge.getPoint(1))).value;
			tmp[1] = index;
			/*
			 * double o = predicates.orientTetra(points[tmp[0]].coords(),
			 * points[tmp[1]].coords(), points[tmp[2]].coords(),
			 * points[tmp[3]].coords()); if (o != 0) {
			 */
			tmpresult.add(tmp);
			/*
			 * }
			 */
		}
		final int[] edge = new int[2 * tmpresult.size()];
		for (i = 0; i < tmpresult.size(); i++) {
			for (int j = 0; j < 2; j++) {
				edge[2 * i + j] = tmpresult.get(i)[j];
			}
		}
		final WB_Triangulation3D result = new WB_Triangulation3D(tetra, tri, edge);
		return result;
	}

	/**
	 *
	 *
	 * @param v
	 * @return
	 */
	private static WB_Point convert(final CVertex v) {
		return geometryfactory.createPoint(v.x(), v.y(), v.z());
	}

	/**
	 *
	 *
	 * @param v
	 * @return
	 */
	private static WB_Point convert(final Point v) {
		return geometryfactory.createPoint(v.x(), v.y(), v.z());
	}

	/**
	 *
	 *
	 * @param points
	 * @return
	 */
	public static WB_AlphaTriangulation3D alphaTriangulate3D(final WB_Coord[] points) {

		final WB_Triangulation3D tri = WB_Triangulate3D.triangulate3D(points);
		return new WB_AlphaTriangulation3D(tri.getTetrahedra(), points);
	}

	/**
	 *
	 *
	 * @param points
	 * @return
	 */
	public static WB_AlphaTriangulation3D alphaTriangulate3D(final Collection<? extends WB_Coord> points) {

		final WB_Triangulation3D tri = WB_Triangulate3D.triangulate3D(points);
		return new WB_AlphaTriangulation3D(tri.getTetrahedra(), points);
	}

	/**
	 *
	 * @param points
	 * @param jitter
	 * @return
	 */
	public static WB_AlphaTriangulation3D alphaTriangulate3D(final WB_Coord[] points, final double jitter) {
		WB_Coord[] jigPoints = Arrays.copyOf(points, points.length);
		WB_RandomOnSphere ros = new WB_RandomOnSphere();
		int i = 0;
		for (WB_Coord p : points) {
			jigPoints[i++] = WB_Point.addMul(p, jitter, ros.nextVector());

		}

		final WB_Triangulation3D tri = WB_Triangulate3D.triangulate3D(jigPoints);
		return new WB_AlphaTriangulation3D(tri.getTetrahedra(), points);
	}

	/**
	 *
	 * @param points
	 * @param jitter
	 * @return
	 */
	public static WB_AlphaTriangulation3D alphaTriangulate3D(final Collection<? extends WB_Coord> points,
			final double jitter) {
		FastTable<WB_Point> jigPoints = new FastTable<WB_Point>();
		WB_RandomOnSphere ros = new WB_RandomOnSphere();
		for (WB_Coord p : points) {
			jigPoints.add(WB_Point.addMul(p, jitter, ros.nextVector()));

		}
		final WB_Triangulation3D tri = WB_Triangulate3D.triangulate3D(jigPoints);
		return new WB_AlphaTriangulation3D(tri.getTetrahedra(), points);
	}
}
