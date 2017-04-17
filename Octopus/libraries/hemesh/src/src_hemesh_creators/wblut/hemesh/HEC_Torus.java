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

import wblut.geom.WB_Point;

/**
 * Torus.
 *
 * @author Frederik Vanhoutte (W:Blut)
 *
 */
public class HEC_Torus extends HEC_Creator {
	/** Tube radius. */
	private double Ri;
	/** Torus Radius. */
	private double Ro;
	/** Facets. */
	private int tubefacets;
	/** Height steps. */
	private int torusfacets;

	private int twist;

	private double tubephase;

	private double torusphase;

	/**
	 *
	 */
	public HEC_Torus() {
		super();
		Ri = 50;
		Ro = 100;
		tubefacets = 6;
		torusfacets = 6;
		tubephase=0.0;
		torusphase=0.0;
	}

	/**
	 * Instantiates a new torus.
	 *
	 * @param Ri
	 *
	 * @param Ro
	 *
	 * @param tubefacets
	 *
	 * @param torusfacets
	 *
	 */
	public HEC_Torus(final double Ri, final double Ro, final int tubefacets,
			final int torusfacets) {
		this();
		this.Ri = Ri;
		this.Ro = Ro;
		this.tubefacets = tubefacets;
		this.torusfacets = torusfacets;
	}

	/**
	 * Sets the radius.
	 *
	 * @param Ri
	 *
	 * @param Ro
	 *
	 * @return
	 */
	public HEC_Torus setRadius(final double Ri, final double Ro) {
		this.Ri = Ri;
		this.Ro = Ro;
		return this;
	}

	/**
	 * Sets the tube facets.
	 *
	 * @param facets
	 *
	 * @return
	 */
	public HEC_Torus setTubeFacets(final int facets) {
		tubefacets = facets;
		return this;
	}

	/**
	 * Sets the torus facets.
	 *
	 * @param facets
	 *
	 * @return
	 */
	public HEC_Torus setTorusFacets(final int facets) {
		torusfacets = facets;
		return this;
	}

	/**
	 * Sets twist.
	 *
	 * @param t
	 *
	 * @return
	 */
	public HEC_Torus setTwist(final int t) {
		twist = t;
		return this;
	}

	/**
	 * Sets torus phase.
	 *
	 * @param p
	 *
	 * @return
	 */
	public HEC_Torus setTorusPhase(final double p) {
		torusphase = p;
		return this;
	}

	/**
	 * Sets tube phase.
	 *
	 * @param p
	 *
	 * @return
	 */
	public HEC_Torus setTubePhase(final double p) {
		tubephase = p;
		return this;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.hemesh.HE_Creator#create()
	 */
	@Override
	protected HE_Mesh createBase() {
		final WB_Point[] vertices = new WB_Point[(tubefacets + 1)
		                                         * (torusfacets + 1)];
		final WB_Point[] uvws = new WB_Point[(tubefacets + 1)
		                                     * (torusfacets + 1)];
		final double dtua = (2 * Math.PI) / tubefacets;
		final double dtoa = (2 * Math.PI) / torusfacets;
		final double dv = 1.0 / tubefacets;
		final double du = 1.0 / torusfacets;
		final double dtwa = (twist * dtoa) / tubefacets;
		int id = 0;
		WB_Point basevertex;
		for (int j = 0; j < (torusfacets + 1); j++) {
			final int lj = (j == torusfacets) ? 0 : j;
			final double ca = Math.cos((lj * dtoa)+torusphase);
			final double sa = Math.sin((lj * dtoa)+torusphase);
			for (int i = 0; i < (tubefacets + 1); i++) {
				final int li = (i == tubefacets) ? 0 : i;
				basevertex = new WB_Point(Ro
						+ (Ri * Math.cos((dtua * li) + (j * dtwa)+tubephase)), 0, Ri
						* Math.sin((dtua * li) + (j * dtwa)+tubephase));
				vertices[id] = new WB_Point(ca * basevertex.xd(), sa
						* basevertex.xd(), basevertex.zd());
				uvws[id] = new WB_Point(j * du, i * dv, 0);
				id++;
			}
		}
		final int nfaces = tubefacets * torusfacets;
		id = 0;
		final int[][] faces = new int[nfaces][];
		int j = 0;
		for (j = 0; j < torusfacets; j++) {
			for (int i = 0; i < tubefacets; i++) {
				faces[id] = new int[4];
				faces[id][0] = i + (j * (tubefacets + 1));
				faces[id][1] = i + ((j + 1) * (tubefacets + 1));
				faces[id][2] = (i + 1) + ((j + 1) * (tubefacets + 1));
				faces[id][3] = (i + 1) + (j * (tubefacets + 1));
				id++;
			}
		}

		final HEC_FromFacelist fl = new HEC_FromFacelist();
		fl.setVertices(vertices).setFaces(faces).setUVW(uvws);
		return fl.createBase();
	}
}
