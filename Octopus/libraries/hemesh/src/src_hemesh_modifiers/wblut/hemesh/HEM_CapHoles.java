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

import javolution.util.FastTable;
import wblut.core.WB_ProgressCounter;

/**
 *
 */
public class HEM_CapHoles extends HEM_Modifier {

	public List<HE_Face> caps = new FastTable<HE_Face>();

	/**
	 *
	 */
	public HEM_CapHoles() {
		super();
		caps = new FastTable<HE_Face>();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.hemesh.HE_Modifier#apply(wblut.hemesh.HE_Mesh)
	 */
	@Override
	public HE_Mesh apply(final HE_Mesh mesh) {
		tracker.setStatus(this, "Starting HEM_CapHoles.", +1);
		caps = new FastTable<HE_Face>();
		final List<HE_Halfedge> unpairedEdges = mesh.getUnpairedHalfedges();
		HE_RAS<HE_Halfedge> loopedHalfedges;
		HE_Halfedge start;
		HE_Halfedge he;
		HE_Halfedge hen;
		HE_Face nf;
		HE_RAS<HE_Halfedge> newHalfedges;
		HE_Halfedge phe;
		HE_Halfedge nhe;
		WB_ProgressCounter counter = new WB_ProgressCounter(unpairedEdges.size(), 10);
		tracker.setStatus(this, "Finding loops and closing holes.", counter);
		while (unpairedEdges.size() > 0) {
			boolean abort = false;
			loopedHalfedges = new HE_RAS.HE_RASTrove<HE_Halfedge>();
			start = unpairedEdges.get(0);
			loopedHalfedges.add(start);
			he = start;
			hen = start;
			boolean noNextFound = false;
			do {
				for (int i = 0; i < unpairedEdges.size(); i++) {
					hen = unpairedEdges.get(i);
					if (hen.getVertex() == he.getNextInFace().getVertex()) {
						if (i > 0 && loopedHalfedges.contains(hen)) {// loop
																		// found
																		// but
																		// not
																		// in
																		// start,
																		// throw
																		// start
																		// out
																		// and
																		// redo
							abort = true;
							unpairedEdges.remove(start);
						}
						loopedHalfedges.add(hen);
						break;
					}
				}
				if (hen.getVertex() != he.getNextInFace().getVertex()) {
					noNextFound = true;
				}
				he = hen;
			} while (hen.getNextInFace().getVertex() != start.getVertex() && !noNextFound && !abort);

			if (!abort) {
				nf = new HE_Face();
				boolean noLoopFound = start.getVertex() != loopedHalfedges.get(loopedHalfedges.size() - 1)
						.getNextInFace().getVertex();
				int ii = 0;
				StringBuilder sb = new StringBuilder(100);
				if (noLoopFound) {
					sb.append("Polyline found: ");
					for (ii = 0; ii < loopedHalfedges.size() - 1; ii++) {
						sb.append(unpairedEdges.indexOf(loopedHalfedges.get(ii)) + "-> ");
					}
					sb.append(unpairedEdges.indexOf(loopedHalfedges.get(ii)));
				} else {
					sb.append("Cycle found: ");
					for (ii = 0; ii < loopedHalfedges.size(); ii++) {
						sb.append(unpairedEdges.indexOf(loopedHalfedges.get(ii)) + "-> ");
					}
					sb.append(unpairedEdges.indexOf(loopedHalfedges.get(0)));
				}
				tracker.setStatus(this, sb.toString(), 0);
				unpairedEdges.removeAll(loopedHalfedges);
				if (!noLoopFound) {
					mesh.add(nf);
					caps.add(nf);
				}
				newHalfedges = new HE_RAS.HE_RASTrove<HE_Halfedge>();
				for (int i = 0; i < loopedHalfedges.size(); i++) {
					phe = loopedHalfedges.get(i);
					nhe = new HE_Halfedge();
					newHalfedges.add(nhe);
					mesh.setVertex(nhe, phe.getNextInFace().getVertex());

					mesh.setPair(nhe, phe);
					mesh.add(nhe);
					if (!noLoopFound) {
						mesh.setFace(nhe, nf);
						if (nf.getHalfedge() == null) {
							mesh.setHalfedge(nf, nhe);
						}
					}
				}
				if (!noLoopFound) {
					mesh.cycleHalfedgesReverse(newHalfedges.getObjects());
				} else {
					mesh.orderHalfedgesReverse(newHalfedges.getObjects());

				}
				counter.increment(newHalfedges.size());
			}
		}
		mesh.cleanUnusedElementsByFace();
		mesh.capHalfedges();

		tracker.setStatus(this, "Exiting HEM_CapHoles.", -1);
		return mesh;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.hemesh.HE_Modifier#apply(wblut.hemesh.HE_Mesh)
	 */
	@Override
	public HE_Mesh apply(final HE_Selection selection) {
		return apply(selection.parent);
	}
}
