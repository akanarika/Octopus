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

import wblut.math.WB_MTRandom;

/**
 * Planar cut of a mesh. Both parts are returned as separate meshes.
 *
 * @author Frederik Vanhoutte (W:Blut)
 *
 */
public class HEMC_Panelizer extends HEMC_MultiCreator {
	/** Source mesh. */
	private HE_Mesh mesh;
	/** The thickness. */
	private double thickness;
	/** The range. */
	private double trange;
	private WB_MTRandom random;
	private double offset;
	/** The range. */
	private double orange;

	public HEMC_Panelizer() {
		super();
		random = new WB_MTRandom();
	}

	/**
	 * Set thickness.
	 *
	 * @param d
	 *            offset
	 * @return self
	 */
	public HEMC_Panelizer setThickness(final double d) {
		thickness = d;
		trange = 0;

		return this;
	}

	public HEMC_Panelizer setOffset(final double d) {
		offset = d;
		orange = 0;

		return this;
	}

	/**
	 * Sets the thickness.
	 *
	 * @param dmin
	 *            the dmin
	 * @param dmax
	 *            the dmax
	 * @return the hEM c_ panelizer
	 */
	public HEMC_Panelizer setThickness(final double dmin, final double dmax) {
		thickness = dmin;
		trange = dmax - dmin;
		return this;
	}

	public HEMC_Panelizer setOffset(final double dmin, final double dmax) {
		offset = dmin;
		orange = dmax - dmin;
		return this;
	}

	public HEMC_Panelizer setSeed(final long seed) {

		random.setSeed(seed);
		return this;
	}

	/**
	 * Set source mesh.
	 *
	 * @param mesh
	 *            mesh to panelize
	 * @return self
	 */
	public HEMC_Panelizer setMesh(final HE_Mesh mesh) {
		this.mesh = mesh;
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

		int id = 0;
		final HEC_Polygon pc = new HEC_Polygon().setThickness(thickness);
		for (final HE_Face f : mesh.getFaces()) {
			pc.setThickness(thickness + trange > 0 ? thickness + random.nextDouble() * trange : 0);
			pc.setOffset(offset + random.nextDouble() * orange);
			pc.setPolygon(f.toPolygon());
			result.add(new HE_Mesh(pc));
			id++;
		}
		_numberOfMeshes = id;
		return result;
	}
}
