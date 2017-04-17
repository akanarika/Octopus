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

import wblut.math.WB_Epsilon;

/**
 *
 */
public class WB_VoronoiCell3D {

	/**
	 *
	 */
	WB_Point generator;

	/**
	 *
	 */
	int index;

	/**
	 *
	 */
	WB_Mesh cell;

	/**
	 *
	 */
	boolean open;

	/**
	 *
	 */
	boolean sliced;

	/**
	 *
	 */
	boolean[] onBoundary;

	/**
	 *
	 */
	private WB_GeometryFactory geometryfactory = new WB_GeometryFactory();

	/**
	 *
	 *
	 * @param points
	 * @param generator
	 * @param index
	 */
	public WB_VoronoiCell3D(final WB_Coord[] points, final WB_Point generator, final int index) {
		this.generator = generator;
		this.index = index;
		cell = geometryfactory.createConvexHull(points, false);
		if (cell != null) {
			onBoundary = new boolean[cell.getNumberOfVertices()];
		}
	}

	/**
	 *
	 *
	 * @param points
	 * @param generator
	 * @param index
	 */
	public WB_VoronoiCell3D(final List<? extends WB_Coord> points, final WB_Point generator, final int index) {
		this.generator = generator;
		this.index = index;
		cell = geometryfactory.createConvexHull(points, false);
		if (cell != null) {
			onBoundary = new boolean[cell.getNumberOfVertices()];
		}
	}

	/**
	 *
	 *
	 * @param cell
	 * @param generator
	 * @param index
	 */
	public WB_VoronoiCell3D(final WB_Mesh cell, final WB_Point generator, final int index) {
		this.generator = generator;
		this.index = index;
		this.cell = cell;
		if (cell != null) {
			onBoundary = new boolean[cell.getNumberOfVertices()];
		}
	}

	/**
	 *
	 *
	 * @param container
	 */
	public void constrain(final WB_AABB container) {
		final WB_AABB aabb = cell.getAABB();
		if (container.contains(aabb)) {
			return;
		}
		if (aabb.intersects(container)) {
			final double[] min = container._min;
			final double[] max = container._max;
			final WB_Point mmm = geometryfactory.createPoint(min[0], min[1], min[2]);
			final WB_Point Mmm = geometryfactory.createPoint(max[0], min[1], min[2]);
			final WB_Point mMm = geometryfactory.createPoint(min[0], max[1], min[2]);
			final WB_Point mmM = geometryfactory.createPoint(min[0], min[1], max[2]);
			final WB_Point MMM = geometryfactory.createPoint(max[0], max[1], max[2]);
			final WB_Point mMM = geometryfactory.createPoint(min[0], max[1], max[2]);
			final WB_Point MmM = geometryfactory.createPoint(max[0], min[1], max[2]);
			final WB_Point MMm = geometryfactory.createPoint(max[0], max[1], min[2]);
			final List<WB_Plane> planes = new ArrayList<WB_Plane>(6);
			planes.add(geometryfactory.createPlane(mmm, Mmm, mMm));
			planes.add(geometryfactory.createPlane(mmm, mMm, mmM));
			planes.add(geometryfactory.createPlane(mmm, mmM, Mmm));
			planes.add(geometryfactory.createPlane(MMM, MmM, mMM));
			planes.add(geometryfactory.createPlane(MMM, mMM, MMm));
			planes.add(geometryfactory.createPlane(MMM, MMm, MmM));
			constrain(planes);
		} else {
			cell = null;
		}
	}

	/**
	 *
	 *
	 * @param convexMesh
	 * @param d
	 */
	public void constrain(final WB_Mesh convexMesh, final double d) {
		constrain(convexMesh.getPlanes(d));
	}

	/**
	 *
	 *
	 * @param convexMesh
	 */
	public void constrain(final WB_Mesh convexMesh) {
		constrain(convexMesh.getPlanes(0));
	}

	/**
	 *
	 *
	 * @param planes
	 */
	public void constrain(final Collection<? extends WB_Plane> planes) {
		for (final WB_Plane P : planes) {
			if (cell != null) {
				slice(P);
			}
		}
		if (cell != null) {
			onBoundary = new boolean[cell.getNumberOfVertices()];
			double d;
			WB_Coord p;
			pointloop: for (int i = 0; i < cell.getNumberOfVertices(); i++) {
				p = cell.getVertex(i);
				for (final WB_Plane P : planes) {
					d = WB_GeometryOp3D.getDistanceToPlane3D(p, P);
					if (WB_Epsilon.isZero(d)) {
						onBoundary[i] = true;
						continue pointloop;
					}
				}
			}
			final int hfl = cell.getNumberOfFaces();
			for (int i = hfl - 1; i > -1; i--) {
				final int[] face = cell.getFace(i);
				boolean boundary = true;
				for (int j = 0; j < face.length; j++) {
					if (!onBoundary[face[j]]) {
						boundary = false;
						break;
					}
				}
				if (boundary) {
					open = true;
				}
			}
		}
	}

	/**
	 *
	 *
	 * @param P
	 */
	private void slice(final WB_Plane P) {
		final WB_Classification[] classifyPoints = ptsPlane(P);
		final List<WB_Coord> newPoints = new ArrayList<WB_Coord>();
		for (int i = 0; i < classifyPoints.length; i++) {
			if (classifyPoints[i] != WB_Classification.BACK) {
				newPoints.add(cell.getVertex(i));
			}
		}
		final int[][] edges = cell.getEdgesAsInt();
		for (final int[] edge : edges) {
			if (((classifyPoints[edge[0]] == WB_Classification.BACK)
					&& (classifyPoints[edge[1]] == WB_Classification.FRONT))
					|| ((classifyPoints[edge[1]] == WB_Classification.BACK)
							&& (classifyPoints[edge[0]] == WB_Classification.FRONT))) {
				final WB_Coord a = cell.getVertex(edge[0]);
				final WB_Coord b = cell.getVertex(edge[1]);
				newPoints.add((WB_Point) WB_GeometryOp3D.getIntersection3D(a, b, P).object);
				sliced = true;
			}
		}
		cell = geometryfactory.createConvexHull(newPoints, false);

	}

	/**
	 *
	 *
	 * @param P
	 * @return
	 */
	private WB_Classification[] ptsPlane(final WB_Plane P) {
		final WB_Classification[] result = new WB_Classification[cell.getNumberOfVertices()];
		for (int i = 0; i < cell.getNumberOfVertices(); i++) {
			result[i] = WB_GeometryOp3D.classifyPointToPlane3D(cell.getVertex(i), P);
		}
		return result;
	}

	/**
	 *
	 *
	 * @return
	 */
	public int getIndex() {
		return index;
	}

	/**
	 *
	 *
	 * @return
	 */
	public WB_Point getGenerator() {
		return generator;
	}

	/**
	 *
	 *
	 * @return
	 */
	public WB_Mesh getMesh() {
		return cell;
	}

	/**
	 *
	 *
	 * @return
	 */
	public boolean[] boundaryArray() {
		return onBoundary;
	}

	/**
	 *
	 *
	 * @return
	 */
	public boolean isOpen() {
		return open;
	}

	/**
	 *
	 *
	 * @return
	 */
	public boolean isSliced() {
		return sliced;
	}
}
