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

import java.util.List;
import java.util.Map;

import javolution.util.FastMap;
import wblut.core.WB_ProgressCounter;
import wblut.geom.WB_Coord;
import wblut.geom.WB_GeometryFactory;

/**
 *
 */
public class HEM_Crocodile extends HEM_Modifier {
	/**
	 *
	 */
	private static WB_GeometryFactory gf = new WB_GeometryFactory();
	/**
	 *
	 */
	private double distance;
	/**
	 *
	 */
	public HE_Selection spikes;
	/**
	 *
	 */
	private double chamfer;

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.hemesh.modifiers.HEB_Modifier#modify(wblut.hemesh.HE_Mesh)
	 */
	/**
	 *
	 */
	public HEM_Crocodile() {
		chamfer = 0.5;
	}

	/**
	 *
	 *
	 * @param d
	 * @return
	 */
	public HEM_Crocodile setDistance(final double d) {
		distance = d;
		return this;
	}

	/**
	 *
	 *
	 * @param c
	 * @return
	 */
	public HEM_Crocodile setChamfer(final double c) {
		chamfer = c;
		return this;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.hemesh.HEM_Modifier#apply(wblut.hemesh.HE_Mesh)
	 */
	@Override
	public HE_Mesh apply(final HE_Mesh mesh) {
		final HE_Selection selection = HE_Selection.selectAllVertices(mesh);
		return apply(selection);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see
	 * wblut.hemesh.modifiers.HEB_Modifier#modifySelected(wblut.hemesh.HE_Mesh)
	 */
	@Override
	public HE_Mesh apply(final HE_Selection selection) {
		spikes = new HE_Selection(selection.parent);
		selection.collectVertices();
		tracker.setStatus(this, "Starting HEM_Crocodile.", +1);
		final Map<Long, WB_Coord> umbrellapoints = new FastMap<Long, WB_Coord>();
		HE_VertexIterator vitr = selection.vItr();
		HE_Vertex v;
		if (chamfer == 0) {
			tracker.setStatus(this, "Chamfer is 0, nothing to do. Exiting HEM_Crocodile.", -1);
			return selection.parent;
		}
		if (chamfer < 0) {
			chamfer *= -1;
		}
		if (chamfer > 0.5 && chamfer < 1.0) {
			chamfer = 1.0 - chamfer;
		} else if (chamfer < 0 || chamfer > 1) {
			tracker.setStatus(this, "Chamfer is outside range (0-0.5), nothing to do. Exiting HEM_Crocodile.", -1);
			return selection.parent;
		}
		if (chamfer == 0.5) {
			WB_ProgressCounter counter = new WB_ProgressCounter(selection.getNumberOfVertices(), 10);
			tracker.setStatus(this, "Enumerating vertex umbrellas.", counter);
			List<HE_Halfedge> star;
			while (vitr.hasNext()) {
				v = vitr.next();
				star = v.getEdgeStar();
				for (final HE_Halfedge e : star) {
					umbrellapoints.put(e.key, e.getEdgeCenter());
				}
				counter.increment();
			}
			counter = new WB_ProgressCounter(umbrellapoints.size(), 10);

			tracker.setStatus(this, "Splitting edges.", counter);
			for (final long e : umbrellapoints.keySet()) {
				HET_MeshOp.splitEdge(selection.parent, selection.parent.getEdgeWithKey(e), umbrellapoints.get(e));
				counter.increment();
			}
		} else {
			List<HE_Halfedge> star;
			WB_ProgressCounter counter = new WB_ProgressCounter(selection.getNumberOfVertices(), 10);

			tracker.setStatus(this, "Enumerating vertex umbrellas.", counter);
			while (vitr.hasNext()) {
				v = vitr.next();
				star = v.getHalfedgeStar();
				for (final HE_Halfedge he : star) {
					umbrellapoints.put(he.key, gf.createInterpolatedPoint(he.getVertex(), he.getEndVertex(), chamfer));
				}
				counter.increment();
			}
			counter = new WB_ProgressCounter(umbrellapoints.size(), 10);

			tracker.setStatus(this, "Splitting edges.", counter);
			for (final long he : umbrellapoints.keySet()) {
				HET_MeshOp.splitEdge(selection.parent, selection.parent.getHalfedgeWithKey(he), umbrellapoints.get(he));
				counter.increment();
			}
		}
		WB_ProgressCounter counter = new WB_ProgressCounter(selection.getNumberOfVertices(), 10);

		tracker.setStatus(this, "Splitting faces.", counter);
		vitr = selection.vItr();
		while (vitr.hasNext()) {
			v = vitr.next();
			final HE_VertexHalfedgeOutCirculator vhoc = new HE_VertexHalfedgeOutCirculator(v);
			HE_Halfedge he;
			while (vhoc.hasNext()) {
				he = vhoc.next();
				if (he.getFace() != null) {
					spikes.union(HET_MeshOp.splitFace(selection.parent, he.getFace(), he.getEndVertex(),
							he.getPrevInVertex().getEndVertex()));
				}
			}
			counter.increment();
			v.addMulSelf(distance, v.getVertexNormal());
		}
		tracker.setStatus(this, "Exiting HEM_Crocodile.", -1);
		return selection.parent;
	}
}
