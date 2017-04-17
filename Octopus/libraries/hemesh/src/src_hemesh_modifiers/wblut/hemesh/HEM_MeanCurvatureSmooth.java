/*
 * This file is part of HE_Mesh, a library for creating and manipulating meshes.
 * It is dedicated to the public domain. To the extent possible under law,
 * I , Frederik Vanhoutte, have waived all copyright and related or neighboring
 * rights.
 *
 * This work is published from Belgium. (http://creativecommons.org/publicdomain/zero/1.0/)
 *
 */
package wblut.hemesh;

import java.util.Iterator;
import java.util.List;

import wblut.core.WB_ProgressCounter;
import wblut.geom.WB_AABB;
import wblut.geom.WB_Coord;
import wblut.geom.WB_Point;
import wblut.geom.WB_Vector;

/**
 *
 */
public class HEM_MeanCurvatureSmooth extends HEM_Modifier {

	/**
	 *
	 */
	private boolean autoRescale;

	/**
	 *
	 */
	private boolean keepBoundary;

	private double lambda;

	/**
	 *
	 */
	private int iter;

	/**
	 *
	 */
	public HEM_MeanCurvatureSmooth() {
		lambda = 0.5;
		iter = 1;
		keepBoundary = false;

	}

	/**
	 *
	 *
	 * @param b
	 * @return
	 */
	public HEM_MeanCurvatureSmooth setAutoRescale(final boolean b) {
		autoRescale = b;
		return this;
	}

	/**
	 *
	 *
	 * @param r
	 * @return
	 */
	public HEM_MeanCurvatureSmooth setIterations(final int r) {
		iter = r;
		return this;
	}

	/**
	 *
	 *
	 * @param b
	 * @return
	 */
	public HEM_MeanCurvatureSmooth setKeepBoundary(final boolean b) {
		keepBoundary = b;
		return this;
	}

	/**
	 *
	 *
	 * @param lambda
	 * @return
	 */
	public HEM_MeanCurvatureSmooth setLambda(final double lambda) {
		this.lambda = lambda;
		return this;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.hemesh.HEM_Modifier#apply(wblut.hemesh.HE_Mesh)
	 */
	@Override
	public HE_Mesh apply(final HE_Mesh mesh) {
		tracker.setStatus(this, "Starting HEM_MeanCurvatureSmooth.", +1);
		WB_AABB box = new WB_AABB();
		if (autoRescale) {
			box = mesh.getAABB();
		}
		final WB_Coord[] newPositions = new WB_Coord[mesh.getNumberOfVertices()];
		if (iter < 1) {
			iter = 1;
		}
		WB_ProgressCounter counter = new WB_ProgressCounter(iter * mesh.getNumberOfVertices(), 10);

		tracker.setStatus(this, "Smoothing vertices.", counter);
		for (int r = 0; r < iter; r++) {
			Iterator<HE_Vertex> vItr = mesh.vItr();
			HE_Vertex v;
			int id = 0;
			WB_Point p;
			while (vItr.hasNext()) {
				v = vItr.next();
				if (v.isBoundary() && keepBoundary) {
					newPositions[id] = v;
				} else {
					p = new WB_Point();
					double factor = 0;
					HE_Halfedge he = v.getHalfedge();
					do {
						double cotana = he.getPrevInFace().getCotan();
						double cotanb = he.getPair().getPrevInFace().getCotan();
						p.addMulSelf(cotana + cotanb, WB_Vector.sub(he.getEndVertex(), v));
						factor += cotana + cotanb;

						he = he.getNextInVertex();
					} while (he != v.getHalfedge());
					newPositions[id] = p.mulSelf(lambda / factor).addSelf(v);

				}
				id++;
				System.out.println();
			}
			vItr = mesh.vItr();
			id = 0;
			while (vItr.hasNext()) {
				vItr.next().set(newPositions[id]);
				id++;
				counter.increment();
			}
		}

		if (autoRescale) {
			mesh.fitInAABB(box);
		}
		tracker.setStatus(this, "Exiting HEM_MeanCurvatureSmooth.", -1);
		return mesh;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see
	 * wblut.hemesh.modifiers.HEB_Modifier#modifySelected(wblut.hemesh.HE_Mesh)
	 */
	@Override
	public HE_Mesh apply(final HE_Selection selection) {
		tracker.setStatus(this, "Starting HEM_MeanCurvatureSmooth.", +1);
		selection.collectVertices();
		WB_AABB box = new WB_AABB();
		if (autoRescale) {
			box = selection.parent.getAABB();
		}
		final WB_Coord[] newPositions = new WB_Coord[selection.getNumberOfVertices()];
		if (iter < 1) {
			iter = 1;
		}
		WB_ProgressCounter counter = new WB_ProgressCounter(iter * selection.getNumberOfVertices(), 10);

		tracker.setStatus(this, "Smoothing vertices.", counter);
		for (int r = 0; r < iter; r++) {
			Iterator<HE_Vertex> vItr = selection.vItr();
			HE_Vertex v;
			HE_Vertex n;
			List<HE_Vertex> neighbors;
			int id = 0;
			while (vItr.hasNext()) {
				v = vItr.next();
				final WB_Point p = new WB_Point();
				if (v.isBoundary() && keepBoundary) {
					newPositions[id] = v;
				} else {
					neighbors = v.getNeighborVertices();
					final Iterator<HE_Vertex> nItr = neighbors.iterator();
					while (nItr.hasNext()) {
						n = nItr.next();
						if (!selection.contains(n)) {
							nItr.remove();
						}
					}

					for (int i = 0; i < neighbors.size(); i++) {
						p.addMulSelf(lambda / neighbors.size(), neighbors.get(i));
					}
					newPositions[id] = p.addMulSelf(1.0 - lambda, v);
				}
				id++;
			}
			vItr = selection.vItr();
			id = 0;
			while (vItr.hasNext()) {
				vItr.next().set(newPositions[id]);
				id++;
				counter.increment();
			}
		}
		if (autoRescale) {
			selection.parent.fitInAABB(box);
		}
		tracker.setStatus(this, "Exiting HEM_MeanCurvatureSmooth.", -1);
		return selection.parent;
	}
}
