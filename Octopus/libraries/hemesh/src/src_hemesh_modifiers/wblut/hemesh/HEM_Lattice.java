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

import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import gnu.trove.iterator.TLongLongIterator;
import gnu.trove.map.TLongLongMap;
import wblut.core.WB_ProgressCounter;
import wblut.geom.WB_Coord;
import wblut.geom.WB_GeometryFactory;
import wblut.geom.WB_Point;
import wblut.geom.WB_Polygon;
import wblut.geom.WB_Vector;

/**
 *
 */
public class HEM_Lattice extends HEM_Modifier {
	/**
	 *
	 */
	private static final WB_GeometryFactory gf = new WB_GeometryFactory();
	/**
	 *
	 */
	private double d;
	/**
	 *
	 */
	private double sew;
	/**
	 *
	 */
	private double hew;
	/**
	 *
	 */
	private double thresholdAngle;
	/**
	 *
	 */
	private boolean fuse;
	/**
	 *
	 */
	private double fuseAngle;
	/**
	 *
	 */
	private double ibulge, obulge;

	/**
	 *
	 */
	public HEM_Lattice() {
		super();
		d = 0;
		sew = 0;
		thresholdAngle = -1;
		fuseAngle = Math.PI / 36;
		fuse = false;
		ibulge = obulge = 0;
	}

	/**
	 *
	 *
	 * @param d
	 * @return
	 */
	public HEM_Lattice setDepth(final double d) {
		this.d = d;
		return this;
	}

	/**
	 *
	 *
	 * @param w
	 * @return
	 */
	public HEM_Lattice setWidth(final double w) {
		sew = 0.5 * w;
		hew = w;
		return this;
	}

	/**
	 *
	 *
	 * @param w
	 * @param hew
	 * @return
	 */
	public HEM_Lattice setWidth(final double w, final double hew) {
		sew = 0.5 * w;
		this.hew = hew;
		return this;
	}

	/**
	 *
	 *
	 * @param d
	 * @return
	 */
	public HEM_Lattice setBulge(final double d) {
		ibulge = obulge = d;
		return this;
	}

	/**
	 *
	 *
	 * @param inner
	 * @param outer
	 * @return
	 */
	public HEM_Lattice setBulge(final double inner, final double outer) {
		ibulge = inner;
		obulge = outer;
		return this;
	}

	/**
	 *
	 *
	 * @param b
	 * @return
	 */
	public HEM_Lattice setFuse(final boolean b) {
		fuse = b;
		return this;
	}

	/**
	 *
	 *
	 * @param a
	 * @return
	 */
	public HEM_Lattice setThresholdAngle(final double a) {
		thresholdAngle = a;
		return this;
	}

	/**
	 *
	 *
	 * @param a
	 * @return
	 */
	public HEM_Lattice setFuseAngle(final double a) {
		fuseAngle = a;
		return this;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.hemesh.HE_Modifier#apply(wblut.hemesh.HE_Mesh)
	 */
	@Override
	public HE_Mesh apply(final HE_Mesh mesh) {
		tracker.setStatus(this, "Starting HEM_Lattice.", +1);
		if (d == 0) {
			tracker.setStatus(this,
					"Can't create with zero thickness, use HEM_PunchHoles instead. Exiting HEM_Lattice.", -1);
			return mesh;
		}
		if (sew == 0) {
			tracker.setStatus(this, "Can't create with zero width. Exiting HEM_Lattice.", -1);
			return mesh;
		}
		final HEM_Extrude extm = new HEM_Extrude().setDistance(0).setRelative(false).setChamfer(sew).setFuse(fuse)
				.setHardEdgeChamfer(hew).setFuseAngle(fuseAngle).setThresholdAngle(thresholdAngle);
		mesh.modify(extm);

		tracker.setStatus(this, "Creating inner mesh.", 0);

		HEC_Copy cc = new HEC_Copy().setMesh(mesh);
		final HE_Mesh innerMesh = cc.create();

		TLongLongMap allheCorrelation = cc.halfedgeCorrelation;

		WB_ProgressCounter counter = new WB_ProgressCounter(mesh.getNumberOfFaces(), 10);
		tracker.setStatus(this, "Creating face correlations.", counter);
		final HashMap<Long, Long> faceCorrelation = new HashMap<Long, Long>();
		final Iterator<HE_Face> fItr1 = mesh.fItr();
		final Iterator<HE_Face> fItr2 = innerMesh.fItr();
		HE_Face f1;
		HE_Face f2;
		while (fItr1.hasNext()) {
			f1 = fItr1.next();
			f2 = fItr2.next();
			faceCorrelation.put(f1.key(), f2.key());
			counter.increment();
		}
		counter = new WB_ProgressCounter(mesh.getNumberOfHalfedges(), 10);
		tracker.setStatus(this, "Creating boundary halfedge correlations.", counter);
		final HashMap<Long, Long> heCorrelation = new HashMap<Long, Long>();
		HE_Halfedge he1;
		HE_Halfedge he2;
		for (TLongLongIterator it = allheCorrelation.iterator(); it.hasNext();) {
			it.advance();
			he1 = mesh.getHalfedgeWithKey(it.key());

			if (he1.getFace() == null) {
				he2 = innerMesh.getHalfedgeWithKey(it.value());
				heCorrelation.put(he1.key(), he2.key());
			}
			counter.increment();
		}
		tracker.setStatus(this, "Shrinking inner mesh.", 0);
		final HEM_VertexExpand expm = new HEM_VertexExpand().setDistance(-d);
		innerMesh.modify(expm);
		HET_MeshOp.flipFaces(innerMesh);
		final int nf = mesh.getNumberOfFaces();
		final HE_Face[] origFaces = mesh.getFacesAsArray();
		mesh.addVertices(innerMesh.getVerticesAsArray());
		mesh.addFaces(innerMesh.getFacesAsArray());
		mesh.addHalfedges(innerMesh.getHalfedgesAsArray());
		HE_Face fo;
		HE_Face fi;
		List<HE_Halfedge> hei;
		List<HE_Halfedge> heo;
		WB_Point[] viPos;
		WB_Polygon poly;
		HE_Halfedge heoc, heic, heon, hein, heio, heoi;
		HE_Face fNew;
		WB_Coord ni;
		WB_Coord no;
		counter = new WB_ProgressCounter(nf, 10);

		tracker.setStatus(this, "Connecting outer and inner faces.", counter);
		for (int i = 0; i < nf; i++) {
			fo = origFaces[i];
			final Long innerKey = faceCorrelation.get(fo.key());
			if (extm.extruded.contains(fo)) {
				fi = mesh.getFaceWithKey(innerKey);
				if (obulge != 0) {
					no = fo.getFaceNormal();
					fo.push(WB_Vector.mul(no, obulge));
				}
				if (ibulge != 0) {
					ni = fi.getFaceNormal();
					fi.push(WB_Vector.mul(ni, ibulge));
				}
				final int nvo = fo.getFaceOrder();
				final int nvi = fi.getFaceOrder();
				hei = fi.getFaceHalfedges();
				viPos = new WB_Point[nvi];
				for (int j = 0; j < nvi; j++) {
					viPos[j] = new WB_Point(hei.get(j).getVertex());
				}
				poly = gf.createSimplePolygon(viPos);
				heo = fo.getFaceHalfedges();
				for (int j = 0; j < nvo; j++) {
					heoc = heo.get(j);
					heon = heo.get((j + 1) % nvo);
					final int cic = poly.closestIndex(heoc.getVertex());
					final int cin = poly.closestIndex(heon.getVertex());
					heic = hei.get(cin);
					hein = hei.get(cic);
					heio = new HE_Halfedge();
					heoi = new HE_Halfedge();
					fNew = new HE_Face();
					mesh.setVertex(heoi, heon.getVertex());
					mesh.setVertex(heio, hein.getVertex());
					mesh.setNext(heoc, heoi);

					mesh.setFace(heoc, fNew);
					if (cic == cin) {
						mesh.setNext(heoi, heio);
						mesh.setFace(heoi, fNew);
					} else {
						mesh.setNext(heoi, heic);
						mesh.setFace(heoi, fNew);
						mesh.setNext(heic, heio);
						mesh.setFace(heic, fNew);
					}
					mesh.setNext(heio, heoc);
					mesh.setFace(heio, fNew);
					mesh.setHalfedge(fNew, heoc);
					mesh.add(heio);
					mesh.add(heoi);
					mesh.add(fNew);
					mesh.remove(fo);
					mesh.remove(fi);
				}
			}
			counter.increment();
		}
		counter = new WB_ProgressCounter(heCorrelation.size(), 10);

		tracker.setStatus(this, "Connecting outer and inner boundaries.", counter);
		final Iterator<Map.Entry<Long, Long>> it = heCorrelation.entrySet().iterator();
		while (it.hasNext()) {
			final Map.Entry<Long, Long> pairs = it.next();
			he1 = mesh.getHalfedgeWithKey(pairs.getKey());
			he2 = mesh.getHalfedgeWithKey(pairs.getValue());
			heio = new HE_Halfedge();
			heoi = new HE_Halfedge();
			mesh.add(heio);
			mesh.add(heoi);
			mesh.setVertex(heio, he1.getPair().getVertex());
			mesh.setVertex(heoi, he2.getPair().getVertex());
			mesh.setNext(he1, heio);
			mesh.setNext(heio, he2);
			mesh.setNext(he2, heoi);
			mesh.setNext(heoi, he1);
			fNew = new HE_Face();
			mesh.add(fNew);
			mesh.setHalfedge(fNew, he1);
			mesh.setFace(he1, fNew);
			mesh.setFace(he2, fNew);
			mesh.setFace(heio, fNew);
			mesh.setFace(heoi, fNew);
			counter.increment();
		}
		mesh.pairHalfedges();
		if (d < 0) {
			HET_MeshOp.flipFaces(mesh);
		}
		tracker.setStatus(this, "Exiting HEM_Lattice.", -1);
		return mesh;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.hemesh.HE_Modifier#apply(wblut.hemesh.HE_Mesh)
	 */
	@Override
	public HE_Mesh apply(final HE_Selection selection) {
		tracker.setStatus(this, "Starting HEM_Lattice.", +1);
		if (d == 0) {
			tracker.setStatus(this,
					"Can't create with zero thickness, use HEM_PunchHoles instead. Exiting HEM_Lattice.", -1);
			return selection.parent;
		}
		if (sew == 0) {
			tracker.setStatus(this, "Can't create with zero width. Exiting HEM_Lattice.", -1);
			return selection.parent;
		}
		final HEM_Extrude extm = new HEM_Extrude().setDistance(0).setRelative(false).setChamfer(sew).setFuse(fuse)
				.setHardEdgeChamfer(hew).setFuseAngle(fuseAngle).setThresholdAngle(thresholdAngle);
		selection.modify(extm);
		tracker.setStatus(this, "Creating inner mesh.", 0);
		HEC_Copy cc = new HEC_Copy().setMesh(selection.parent);
		final HE_Mesh innerMesh = cc.create();
		TLongLongMap allheCorrelation = cc.halfedgeCorrelation;

		WB_ProgressCounter counter = new WB_ProgressCounter(selection.parent.getNumberOfFaces(), 10);

		tracker.setStatus(this, "Creating face correlations.", counter);
		final HashMap<Long, Long> faceCorrelation = new HashMap<Long, Long>();
		final Iterator<HE_Face> fItr1 = selection.parent.fItr();
		final Iterator<HE_Face> fItr2 = innerMesh.fItr();
		HE_Face f1;
		HE_Face f2;
		while (fItr1.hasNext()) {
			f1 = fItr1.next();
			f2 = fItr2.next();
			faceCorrelation.put(f1.key(), f2.key());
			counter.increment();
		}
		counter = new WB_ProgressCounter(selection.parent.getNumberOfHalfedges(), 10);

		tracker.setStatus(this, "Creating boundary halfedge correlations.", counter);
		final HashMap<Long, Long> heCorrelation = new HashMap<Long, Long>();
		HE_Halfedge he1;
		HE_Halfedge he2;
		for (TLongLongIterator it = allheCorrelation.iterator(); it.hasNext();) {
			it.advance();
			he1 = selection.parent.getHalfedgeWithKey(it.key());

			if (he1.getFace() == null) {
				he2 = innerMesh.getHalfedgeWithKey(it.value());
				heCorrelation.put(he1.key(), he2.key());
			}
			counter.increment();
		}
		tracker.setStatus(this, "Shrinking inner mesh.", 0);
		final HEM_VertexExpand expm = new HEM_VertexExpand().setDistance(-d);
		innerMesh.modify(expm);
		HET_MeshOp.flipFaces(innerMesh);
		final int nf = selection.parent.getNumberOfFaces();
		final HE_Face[] origFaces = selection.parent.getFacesAsArray();
		selection.parent.addVertices(innerMesh.getVerticesAsArray());
		selection.parent.addFaces(innerMesh.getFacesAsArray());
		selection.parent.addHalfedges(innerMesh.getHalfedgesAsArray());
		HE_Face fo;
		HE_Face fi;
		List<HE_Halfedge> hei;
		List<HE_Halfedge> heo;
		WB_Point[] viPos;
		WB_Polygon poly;
		HE_Halfedge heoc, heic, heon, hein, heio, heoi;
		HE_Face fNew;
		WB_Coord ni, no;
		counter = new WB_ProgressCounter(nf, 10);

		tracker.setStatus(this, "Connecting outer and inner faces.", counter);
		for (int i = 0; i < nf; i++) {
			fo = origFaces[i];
			final Long innerKey = faceCorrelation.get(fo.key());
			if (extm.extruded.contains(fo)) {
				fi = selection.parent.getFaceWithKey(innerKey);
				if (obulge != 0) {
					no = fo.getFaceNormal();
					fo.push(WB_Vector.mul(no, obulge));
				}
				if (ibulge != 0) {
					ni = fi.getFaceNormal();
					fi.push(WB_Vector.mul(ni, ibulge));
				}
				final int nvo = fo.getFaceOrder();
				final int nvi = fi.getFaceOrder();
				hei = fi.getFaceHalfedges();
				viPos = new WB_Point[nvi];
				for (int j = 0; j < nvi; j++) {
					viPos[j] = new WB_Point(hei.get(j).getVertex());
				}
				poly = gf.createSimplePolygon(viPos);
				heo = fo.getFaceHalfedges();
				for (int j = 0; j < nvo; j++) {
					heoc = heo.get(j);
					heon = heo.get((j + 1) % nvo);
					final int cic = poly.closestIndex(heoc.getVertex());
					final int cin = poly.closestIndex(heon.getVertex());
					heic = hei.get(cin);
					hein = hei.get(cic);
					heio = new HE_Halfedge();
					heoi = new HE_Halfedge();
					fNew = new HE_Face();
					selection.parent.setVertex(heoi, heon.getVertex());
					selection.parent.setVertex(heio, hein.getVertex());
					selection.parent.setNext(heoc, heoi);
					selection.parent.setFace(heoc, fNew);
					if (cic == cin) {
						selection.parent.setNext(heoi, heio);
						selection.parent.setFace(heoi, fNew);
					} else {
						selection.parent.setNext(heoi, heic);
						selection.parent.setFace(heoi, fNew);
						selection.parent.setNext(heic, heio);
						selection.parent.setFace(heic, fNew);
					}
					selection.parent.setNext(heio, heoc);
					selection.parent.setFace(heio, fNew);
					selection.parent.setHalfedge(fNew, heoc);
					selection.parent.add(heio);
					selection.parent.add(heoi);
					selection.parent.add(fNew);
					selection.parent.remove(fo);
					selection.parent.remove(fi);
				}
			}
			counter.increment();
		}
		counter = new WB_ProgressCounter(heCorrelation.size(), 10);

		tracker.setStatus(this, "Connecting outer and inner boundaries.", counter);
		final Iterator<Map.Entry<Long, Long>> it = heCorrelation.entrySet().iterator();
		while (it.hasNext()) {
			final Map.Entry<Long, Long> pairs = it.next();
			he1 = selection.parent.getHalfedgeWithKey(pairs.getKey());
			he2 = selection.parent.getHalfedgeWithKey(pairs.getValue());
			heio = new HE_Halfedge();
			heoi = new HE_Halfedge();
			selection.parent.add(heio);
			selection.parent.add(heoi);
			selection.parent.setVertex(heio, he1.getPair().getVertex());
			selection.parent.setVertex(heoi, he2.getPair().getVertex());
			selection.parent.setNext(he1, heio);
			selection.parent.setNext(heio, he2);
			selection.parent.setNext(he2, heoi);
			selection.parent.setNext(heoi, he1);
			fNew = new HE_Face();
			selection.parent.add(fNew);
			selection.parent.setHalfedge(fNew, he1);
			selection.parent.setFace(he1, fNew);
			selection.parent.setFace(he2, fNew);
			selection.parent.setFace(heio, fNew);
			selection.parent.setFace(heoi, fNew);
			counter.increment();
		}
		selection.parent.pairHalfedges();
		if (d < 0) {
			HET_MeshOp.flipFaces(selection.parent);
		}
		tracker.setStatus(this, "Exiting HEM_Lattice.", -1);
		return selection.parent;
	}
}
