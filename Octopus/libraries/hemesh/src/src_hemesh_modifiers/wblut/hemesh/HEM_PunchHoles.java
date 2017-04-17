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
public class HEM_PunchHoles extends HEM_Modifier {

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
	public HEM_PunchHoles() {
		super();
		sew = 0;
		thresholdAngle = -1;
		fuseAngle = Math.PI / 36;
		fuse = false;
	}

	/**
	 * 
	 *
	 * @param w
	 * @return
	 */
	public HEM_PunchHoles setWidth(final double w) {
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
	public HEM_PunchHoles setWidth(final double w, final double hew) {
		sew = 0.5 * w;
		this.hew = hew;
		return this;
	}

	/**
	 * 
	 *
	 * @param b
	 * @return
	 */
	public HEM_PunchHoles setFuse(final boolean b) {
		fuse = b;
		return this;
	}

	/**
	 * 
	 *
	 * @param a
	 * @return
	 */
	public HEM_PunchHoles setThresholdAngle(final double a) {
		thresholdAngle = a;
		return this;
	}

	/**
	 * 
	 *
	 * @param a
	 * @return
	 */
	public HEM_PunchHoles setFuseAngle(final double a) {
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
		if (sew == 0) {
			return mesh;
		}
		final HEM_Extrude extm = new HEM_Extrude().setDistance(0).setRelative(false).setChamfer(sew).setFuse(fuse)
				.setHardEdgeChamfer(hew).setFuseAngle(fuseAngle).setThresholdAngle(thresholdAngle);
		mesh.modify(extm);
		mesh.deleteFaces(extm.extruded);
		return mesh;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see wblut.hemesh.HE_Modifier#apply(wblut.hemesh.HE_Mesh)
	 */
	@Override
	public HE_Mesh apply(final HE_Selection selection) {
		if (sew == 0) {
			return selection.parent;
		}
		final HEM_Extrude extm = new HEM_Extrude().setDistance(0).setRelative(false).setChamfer(sew).setFuse(fuse)
				.setHardEdgeChamfer(hew).setFuseAngle(fuseAngle).setThresholdAngle(thresholdAngle);
		selection.modify(extm);
		selection.parent.deleteFaces(extm.extruded);
		return selection.parent;
	}
}
