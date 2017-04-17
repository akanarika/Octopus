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

import wblut.geom.WB_Plane;

/**
 * Planar cut of a mesh. Both parts are returned as separate meshes.
 *
 * @author Frederik Vanhoutte (W:Blut)
 *
 */
public class HEMC_SplitMesh extends HEMC_MultiCreator {
	/** Cutting plane. */
	private WB_Plane P;
	/** Source mesh. */
	private HE_Mesh mesh;
	/** Cap holes?. */
	private boolean cap = true;;
	/** The offset. */
	private double offset;
	private boolean simpleCap = true;

	/**
	 * Set offset.
	 *
	 * @param d
	 *            offset
	 * @return self
	 */
	public HEMC_SplitMesh setOffset(final double d) {
		offset = d;
		return this;
	}

	/**
	 * Instantiates a new HEMC_SplitMesh.
	 *
	 */
	public HEMC_SplitMesh() {
		super();
	}

	/**
	 * Set split plane.
	 *
	 * @param P
	 *            plane
	 * @return self
	 */
	public HEMC_SplitMesh setPlane(final WB_Plane P) {
		this.P = P;
		return this;
	}

	/**
	 * Set source mesh.
	 *
	 * @param mesh
	 *            mesh to split
	 * @return self
	 */
	public HEMC_SplitMesh setMesh(final HE_Mesh mesh) {
		this.mesh = mesh;
		return this;
	}

	/**
	 * Set option to cap holes.
	 *
	 * @param b
	 *            true, false;
	 * @return self
	 */
	public HEMC_SplitMesh setCap(final Boolean b) {
		cap = b;
		return this;
	}

	public HEMC_SplitMesh setSimpleCap(final Boolean b) {
		simpleCap = b;
		return this;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.hemesh.HE_MultiCreator#create()
	 */
	@Override
	public HE_MeshCollection create() {
		final HE_MeshCollection result = new HE_MeshCollection();
		if (mesh == null) {
			_numberOfMeshes = 0;
			return result;
		}
		result.add(mesh.copy());
		if (P == null) {
			_numberOfMeshes = 1;
			return result;
		}
		final HEM_Slice sm = new HEM_Slice();
		sm.setPlane(P).setReverse(false).setCap(cap).setOffset(offset).setSimpleCap(simpleCap);
		sm.apply(result.getMesh(0));
		result.getMesh(0).resetFaceInternalLabels();
		HE_FaceIterator fItr = sm.capFaces.fItr();
		while (fItr.hasNext()) {
			fItr.next().setInternalLabel(1);
		}
		P.flipNormal();
		sm.setPlane(P).setReverse(false).setCap(cap).setOffset(offset).setSimpleCap(simpleCap);
		result.add(mesh.copy());
		sm.apply(result.getMesh(1));
		result.getMesh(0).resetFaceInternalLabels();
		fItr = sm.capFaces.fItr();
		while (fItr.hasNext()) {
			fItr.next().setInternalLabel(1);
		}
		_numberOfMeshes = 2;
		return result;
	}
}
