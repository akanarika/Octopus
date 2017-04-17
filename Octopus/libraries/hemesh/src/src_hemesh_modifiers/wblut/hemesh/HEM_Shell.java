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

import gnu.trove.iterator.TLongLongIterator;
import gnu.trove.map.TLongLongMap;

/**
 * Turns a solid into a rudimentary shelled structure.
 *
 * @author Frederik Vanhoutte (W:Blut)
 *
 */
public class HEM_Shell extends HEM_Modifier {
	/**
	 *
	 */
	private double d;

	/**
	 *
	 */
	public HEM_Shell() {
		super();
		d = 0;
	}

	/**
	 *
	 *
	 * @param d
	 * @return
	 */
	public HEM_Shell setThickness(final double d) {
		this.d = d;
		return this;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.hemesh.HE_Modifier#apply(wblut.hemesh.HE_Mesh)
	 */
	@Override
	public HE_Mesh apply(final HE_Mesh mesh) {
		if (d == 0) {
			return mesh;
		}

		HEC_Copy cc = new HEC_Copy().setMesh(mesh);
		final HE_Mesh innerMesh = cc.create();

		TLongLongMap heCorrelation = cc.halfedgeCorrelation;

		final HEM_VertexExpand expm = new HEM_VertexExpand().setDistance(-d);
		innerMesh.modify(expm);
		HET_MeshOp.flipFaces(innerMesh);
		mesh.add(innerMesh);
		HE_Halfedge he1, he2, heio, heoi;
		HE_Face fNew;
		for (TLongLongIterator it = heCorrelation.iterator(); it.hasNext();) {
			it.advance();
			he1 = mesh.getHalfedgeWithKey(it.key());
			if (he1.isOuterBoundary()) {
				he2 = mesh.getHalfedgeWithKey(it.value());
				heio = new HE_Halfedge();
				heoi = new HE_Halfedge();
				mesh.setVertex(heio, he1.getPair().getVertex());
				heio.setUVW(he1.getPair().getVertex().getUVW(he1.getPair().getFace()));
				mesh.setVertex(heoi, he2.getPair().getVertex());
				heoi.setUVW(he2.getPair().getVertex().getUVW(he2.getPair().getFace()));
				mesh.setNext(he1, heio);
				mesh.setNext(heio, he2);
				mesh.setNext(he2, heoi);
				mesh.setNext(heoi, he1);
				fNew = new HE_Face();
				fNew.setInternalLabel(1);
				mesh.add(fNew);
				mesh.setHalfedge(fNew, he1);
				mesh.setFace(he1, fNew);
				mesh.setFace(he2, fNew);
				mesh.setFace(heio, fNew);
				mesh.setFace(heoi, fNew);
				mesh.add(heio);
				mesh.add(heoi);
			}

		}
		mesh.pairHalfedges();
		mesh.capHalfedges();
		if (d < 0) {
			HET_MeshOp.flipFaces(mesh);
		}
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
