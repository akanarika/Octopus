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

/**
 *
 */
public class HEM_CenterSplit extends HEM_Modifier {

	/**
	 *
	 */
	private double d;

	/**
	 *
	 */
	private double c;

	/**
	 *
	 */
	private HE_Selection selectionOut;

	/**
	 *
	 */
	public HEM_CenterSplit() {
		super();
		d = 0;
		c = 0.5;
	}

	/**
	 *
	 *
	 * @param d
	 * @return
	 */
	public HEM_CenterSplit setOffset(final double d) {
		this.d = d;
		return this;
	}

	/**
	 *
	 *
	 * @param c
	 * @return
	 */
	public HEM_CenterSplit setChamfer(final double c) {
		this.c = c;
		return this;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.hemesh.HE_Modifier#apply(wblut.hemesh.HE_Mesh)
	 */
	@Override
	public HE_Mesh apply(final HE_Mesh mesh) {
		tracker.setStatus(this, "Starting HEC_CenterSplit.", +1);
		final HEM_Extrude ext = new HEM_Extrude().setChamfer(c).setDistance(d);
		mesh.modify(ext);
		selectionOut = ext.extruded;
		tracker.setStatus(this, "Exiting HEC_CenterSplit.", -1);
		return mesh;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.hemesh.HE_Modifier#apply(wblut.hemesh.HE_Mesh)
	 */
	@Override
	public HE_Mesh apply(final HE_Selection selection) {
		tracker.setStatus(this, "Starting HEC_CenterSplit.", +1);
		final HEM_Extrude ext = new HEM_Extrude().setChamfer(c).setDistance(d);
		selection.modify(ext);
		selectionOut = ext.extruded;
		tracker.setStatus(this, "Exiting HEC_CenterSplit.", -1);
		return selection.parent;
	}

	/**
	 *
	 *
	 * @return
	 */
	public HE_Selection getCenterFaces() {
		return this.selectionOut;
	}
}
