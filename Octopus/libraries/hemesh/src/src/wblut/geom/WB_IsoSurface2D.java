/*
 * This file is part of HE_Mesh, a library for creating and manipulating meshes.
 * It is dedicated to the public domain. To the extent possible under law,
 * I , Frederik Vanhoutte, have waived all copyright and related or neighboring
 * rights.
 *
 * This work is published from Belgium. (http://creativecommons.org/publicdomain/zero/1.0/)
 *
 */
package wblut.geom;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;

import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import javolution.util.FastTable;
import wblut.math.WB_Epsilon;

/**
 *
 */
public class WB_IsoSurface2D {
	/**
	 *
	 */
	final static int ONVERTEX = 0;
	/**
	 *
	 */
	final static int ONEDGE = 1;
	/**
	 *
	 */
	final static int NEGATIVE = 0;
	/**
	 *
	 */
	final static int EQUAL = 1;
	/**
	 *
	 */
	final static int POSITIVE = 2;
	/**
	 *
	 */
	int[] digits = new int[8];
	/*
	 * VERTICES 000 ijk=0 100 Ijk=1 010 iJk=2 110 IJk=3 001 ijK=4 101 IjK=5 011
	 * iJK=6 111 IJK=7
	 */
	/**
	 *
	 */
	final static WB_Point[] gridvertices = new WB_Point[] { new WB_Point(0, 0), new WB_Point(1, 0), new WB_Point(0, 1),
			new WB_Point(1, 1) };
	// EDGES: 2 vertices per edge
	/**
	 *
	 */
	final static int[][] edges = { { 0, 1 }, // x ijk
			{ 0, 2 }, // y ijk
			{ 1, 3 }, // y Ijk
			{ 2, 3 }, // x iJk
	};
	// ISOVERTICES: 8
	// type=ONVERTEX iso vertex on vertex, index in vertex list
	// type=ONEDGE iso vertex on edge, index in edge list
	/**
	 *
	 */
	final static int[][] isovertices = new int[][] { { 0, 0 }, { 0, 1 }, { 0, 2 }, { 0, 3 }, { 1, 0 }, { 1, 1 },
			{ 1, 2 }, { 1, 3 } };
	/**
	 *
	 */
	int[][] entries;
	/**
	 *
	 */
	private double[][] values;
	/**
	 *
	 */
	private int resx, resy;
	/**
	 *
	 */
	private double cx, cy;
	/**
	 *
	 */
	private double dx, dy;
	/**
	 *
	 */
	private double isolevel;
	/**
	 *
	 */
	private double boundary;
	/**
	 *
	 */
	private TIntObjectMap<WB_Point> xedges;
	/**
	 *
	 */
	private TIntObjectMap<WB_Point> yedges;
	/**
	 *
	 */
	private TIntObjectMap<WB_Point> vertices;
	/**
	 *
	 */
	private TIntObjectMap<VertexRemap> vertexremaps;
	/**
	 *
	 */
	double gamma;

	/**
	 *
	 */
	private boolean invert;
	FastTable<WB_Point> segs;

	/**
	 *
	 */
	public WB_IsoSurface2D() {
		String line = "";
		final String cvsSplitBy = " ";
		BufferedReader br = null;
		InputStream is = null;
		InputStreamReader isr = null;
		entries = new int[6561][];
		try {
			is = this.getClass().getClassLoader().getResourceAsStream("resources/isonepcube2D.txt");
			isr = new InputStreamReader(is);
			br = new BufferedReader(isr);
			int i = 0;
			while ((line = br.readLine()) != null) {
				final String[] cell = line.split(cvsSplitBy);
				final int[] indices = new int[cell.length];
				for (int j = 0; j < cell.length; j++) {
					indices[j] = Integer.parseInt(cell[j]);
				}
				entries[i] = indices;
				i++;
			}
		} catch (final FileNotFoundException e) {
			e.printStackTrace();
		} catch (final IOException e) {
			e.printStackTrace();
		} finally {
			if (br != null) {
				try {
					br.close();
					isr.close();
					is.close();
				} catch (final IOException e) {
					e.printStackTrace();
				}
			}
		}
		gamma = 0;
		boundary = Double.NaN;
	}

	/**
	 *
	 *
	 * @param gamma
	 * @return
	 */
	public WB_IsoSurface2D setGamma(final double gamma) {
		this.gamma = gamma;
		return this;
	}

	/**
	 * Number of cells.
	 *
	 * @param resx
	 *            the resx
	 * @param resy
	 *            the resy
	 *
	 * @return self
	 */
	public WB_IsoSurface2D setResolution(final int resx, final int resy) {
		this.resx = resx;
		this.resy = resy;
		return this;
	}

	/**
	 * Size of cell.
	 *
	 * @param dx
	 * @param dy
	 * @return self
	 */
	public WB_IsoSurface2D setSize(final double dx, final double dy) {
		this.dx = dx;
		this.dy = dy;
		return this;
	}

	/**
	 * Values at grid points.
	 *
	 * @param values
	 *            double[resx+1][resy+1]
	 * @return self
	 */
	public WB_IsoSurface2D setValues(final double[][] values) {
		this.values = new double[resx + 1][resy + 1];
		for (int i = 0; i <= resx; i++) {
			for (int j = 0; j <= resy; j++) {

				this.values[i][j] = values[i][j];

			}
		}
		return this;
	}

	/**
	 * Sets the values.
	 *
	 * @param values
	 *            float[resx+1][resy+1]
	 * @return self
	 */
	public WB_IsoSurface2D setValues(final float[][] values) {
		this.values = new double[resx + 1][resy + 1];
		for (int i = 0; i <= resx; i++) {
			for (int j = 0; j <= resy; j++) {
				this.values[i][j] = values[i][j];
			}
		}
		return this;
	}

	/**
	 * Isolevel to render.
	 *
	 * @param v
	 *            isolevel
	 * @return self
	 */
	public WB_IsoSurface2D setIsolevel(final double v) {
		isolevel = v;
		return this;
	}

	/**
	 * Boundary level.
	 *
	 * @param v
	 *            boundary level
	 * @return self
	 */
	public WB_IsoSurface2D setBoundary(final double v) {
		boundary = v;
		return this;
	}

	/**
	 * Clear boundary level.
	 *
	 * @return self
	 */
	public WB_IsoSurface2D clearBoundary() {
		boundary = Double.NaN;
		return this;
	}

	/**
	 * Invert isosurface.
	 *
	 * @param invert
	 *            true/false
	 * @return self
	 */
	public WB_IsoSurface2D setInvert(final boolean invert) {
		this.invert = invert;
		return this;
	}

	public WB_IsoSurface2D setCenter(final WB_Coord c) {
		cx = c.xd();
		cy = c.yd();

		return this;
	}

	public List<WB_Point> create() {
		vertices = new TIntObjectHashMap<WB_Point>(1024, 0.5f, -1);
		xedges = new TIntObjectHashMap<WB_Point>(1024, 0.5f, -1);
		yedges = new TIntObjectHashMap<WB_Point>(1024, 0.5f, -1);
		if (gamma > 0) {
			mapvertices();
			setvalues();
			polygonise();
			snapvertices();
			resetvalues();
		} else {
			polygonise();

		}
		return segs;
	}

	/**
	 *
	 */
	private void mapvertices() {

		vertexremaps = new TIntObjectHashMap<VertexRemap>(1024, 0.5f, -1);
		final WB_Point offset = new WB_Point(cx - 0.5 * resx * dx, cy - 0.5 * resy * dy);
		if (Double.isNaN(boundary)) {
			for (int i = 0; i < resx; i++) {
				// System.out.println("HEC_IsoSurface: " + (i + 1) + " of " +
				// resx);
				for (int j = 0; j < resy; j++) {

					getPolygons(i, j, classifyCell(i, j), offset, true);

				}
			}
		} else {
			for (int i = -1; i < resx + 1; i++) {
				// System.out.println("HEC_IsoSurface: " + (i + 1) + " of " +
				// resx);
				for (int j = -1; j < resy + 1; j++) {

					getPolygons(i, j, classifyCell(i, j), offset, true);

				}
			}
		}
	}

	/**
	 *
	 */
	void setvalues() {
		VertexRemap vr;
		for (final Object o : vertexremaps.values()) {
			vr = (VertexRemap) o;
			vr.snapvertex.set(vr.p);
			values[vr.i][vr.j] = isolevel;
		}
	}

	/**
	 * Polygonise.
	 */
	private void polygonise() {

		segs = new FastTable<WB_Point>();
		final WB_Point offset = new WB_Point(cx - 0.5 * resx * dx, cy - 0.5 * resy * dy);
		if (Double.isNaN(boundary)) {
			for (int i = 0; i < resx; i++) {
				// System.out.println("HEC_IsoSurface: " + (i + 1) + " of " +
				// resx);
				for (int j = 0; j < resy; j++) {

					getPolygons(i, j, classifyCell(i, j), offset, false);

				}
			}
		} else {
			for (int i = -1; i < resx + 1; i++) {
				// System.out.println("HEC_IsoSurface: " + (i + 1) + " of " +
				// resx);
				for (int j = -1; j < resy + 1; j++) {

					getPolygons(i, j, classifyCell(i, j), offset, false);

				}
			}
		}
	}

	/**
	 *
	 */
	void snapvertices() {
		VertexRemap vr;
		for (final Object o : vertexremaps.values()) {
			vr = (VertexRemap) o;
			vr.snapvertex.set(vr.p);
		}
	}

	/**
	 *
	 */
	void resetvalues() {
		VertexRemap vr;
		for (final Object o : vertexremaps.values()) {
			vr = (VertexRemap) o;
			values[vr.i][vr.j] = vr.originalvalue;
		}
	}

	/**
	 * Gets the polygons.
	 *
	 * @param i
	 *            the i
	 * @param j
	 *            the j
	 *
	 * @param cubeindex
	 *            the cubeindex
	 * @param offset
	 * @param dummyrun
	 * @return the polygons
	 */
	private List<WB_Point> getPolygons(final int i, final int j, final int cubeindex, final WB_Point offset,
			final boolean dummyrun) {
		final int[] indices = entries[cubeindex];
		final int numsegs = indices[0];
		int currentindex = 1;
		for (int t = 0; t < numsegs; t++) {
			final WB_Point v2 = getIsoVertex(indices[currentindex++], i, j, offset, dummyrun);
			final WB_Point v1 = getIsoVertex(indices[currentindex++], i, j, offset, dummyrun);
			if (!dummyrun) {
				segs.add(v2);
				segs.add(v1);
			}
		}
		return segs;
	}

	/**
	 *
	 * @param i
	 * @param j
	 * @return
	 */
	private int classifyCell(final int i, final int j) {
		if (Double.isNaN(boundary)) {
			if (i < 0 || j < 0 || i >= resx || j >= resy) {
				return -1;
			}
		}
		digits = new int[8];
		int cubeindex = 0;
		int offset = 1;
		if (invert) {
			if (value(i, j) < isolevel) {
				cubeindex += 2 * offset;
				digits[0] = POSITIVE;
			} else if (value(i, j) == isolevel) {
				cubeindex += offset;
				digits[0] = EQUAL;
			}
			offset *= 3;
			if (value(i + 1, j) < isolevel) {
				cubeindex += 2 * offset;
				digits[1] = POSITIVE;
			} else if (value(i + 1, j) == isolevel) {
				cubeindex += offset;
				digits[1] = EQUAL;
			}
			offset *= 3;
			if (value(i, j + 1) < isolevel) {
				cubeindex += 2 * offset;
				digits[2] = POSITIVE;
			} else if (value(i, j + 1) == isolevel) {
				cubeindex += offset;
				digits[2] = EQUAL;
			}
			offset *= 3;
			if (value(i + 1, j + 1) < isolevel) {
				cubeindex += 2 * offset;
				digits[3] = POSITIVE;
			} else if (value(i + 1, j + 1) == isolevel) {
				cubeindex += offset;
				digits[3] = EQUAL;
			}

		} else {
			if (value(i, j) > isolevel) {
				cubeindex += 2 * offset;
				digits[0] = POSITIVE;
			} else if (value(i, j) == isolevel) {
				cubeindex += offset;
				digits[0] = EQUAL;
			}
			offset *= 3;
			if (value(i + 1, j) > isolevel) {
				cubeindex += 2 * offset;
				digits[1] = POSITIVE;
			} else if (value(i + 1, j) == isolevel) {
				cubeindex += offset;
				digits[1] = EQUAL;
			}
			offset *= 3;
			if (value(i, j + 1) > isolevel) {
				cubeindex += 2 * offset;
				digits[2] = POSITIVE;
			} else if (value(i, j + 1) == isolevel) {
				cubeindex += offset;
				digits[2] = EQUAL;
			}
			offset *= 3;
			if (value(i + 1, j + 1) > isolevel) {
				cubeindex += 2 * offset;
				digits[3] = POSITIVE;
			} else if (value(i + 1, j + 1) == isolevel) {
				cubeindex += offset;
				digits[3] = EQUAL;
			}

		}
		return cubeindex;
	}

	/**
	 *
	 *
	 * @param isopointindex
	 * @param i
	 * @param j
	 * @param k
	 * @param offset
	 * @param dummyrun
	 * @return
	 */
	WB_Point getIsoVertex(final int isopointindex, final int i, final int j, final WB_Point offset,
			final boolean dummyrun) {
		if (isovertices[isopointindex][0] == ONVERTEX) {
			switch (isovertices[isopointindex][1]) {
			case 0:
				return vertex(i, j, offset);
			case 1:
				return vertex(i + 1, j, offset);
			case 2:
				return vertex(i, j + 1, offset);
			case 3:
				return vertex(i + 1, j + 1, offset);

			default:
				return null;
			}
		} else if (isovertices[isopointindex][0] == ONEDGE) {
			switch (isovertices[isopointindex][1]) {
			case 0:
				return xedge(i, j, offset, dummyrun);
			case 1:
				return yedge(i, j, offset, dummyrun);
			case 2:
				return yedge(i + 1, j, offset, dummyrun);
			case 3:
				return xedge(i, j + 1, offset, dummyrun);

			default:
				return null;
			}
		}
		return null;
	}

	/**
	 *
	 *
	 * @param i
	 * @param j
	 * @param offset
	 * @return
	 */
	private WB_Point vertex(final int i, final int j, final WB_Point offset) {
		WB_Point vertex = vertices.get(index(i, j));
		if (vertex != null) {
			return vertex;
		}
		final WB_Point p0 = new WB_Point(i * dx, j * dy);
		vertex = new WB_Point(p0.addSelf(offset));
		vertices.put(index(i, j), vertex);
		return vertex;
	}

	/**
	 * Xedge.
	 *
	 * @param i
	 *            i: -1 .. resx+1
	 * @param j
	 *            j: -1 .. resy+1
	 *
	 * @param offset
	 * @param dummyrun
	 * @return edge vertex
	 */
	private WB_Point xedge(final int i, final int j, final WB_Point offset, final boolean dummyrun) {
		final WB_Point p0 = new WB_Point(i * dx, j * dy);
		final WB_Point p1 = new WB_Point(i * dx + dx, j * dy);
		final double val0 = value(i, j);
		final double val1 = value(i + 1, j);
		double mu;
		if (dummyrun) {
			mu = (isolevel - val0) / (val1 - val0);
			if (mu < gamma) {
				VertexRemap vr = vertexremaps.get(index(i, j));
				if (vr == null) {
					vr = new VertexRemap();
					vr.closestd = mu * dx;
					vr.i = i;
					vr.j = j;

					vr.originalvalue = values[i][j];
					vr.p = interp(isolevel, p0, p1, val0, val1).addSelf(offset);
					vr.snapvertex = vertex(i, j, offset);
					vertexremaps.put(index(i, j), vr);
				} else {
					if (vr.closestd > mu * dx) {
						vr.closestd = mu * dx;
						vr.i = i;
						vr.j = j;

						vr.originalvalue = values[i][j];
						vr.p = interp(isolevel, p0, p1, val0, val1).addSelf(offset);
						vr.snapvertex = vertex(i, j, offset);
					}
				}
			} else if (mu > 1 - gamma) {
				VertexRemap vr = vertexremaps.get(index(i + 1, j));
				if (vr == null) {
					vr = new VertexRemap();
					vr.closestd = (1 - mu) * dx;
					vr.i = i + 1;
					vr.j = j;

					vr.originalvalue = values[i + 1][j];
					vr.p = interp(isolevel, p0, p1, val0, val1).addSelf(offset);
					vr.snapvertex = vertex(i + 1, j, offset);
					vertexremaps.put(index(i + 1, j), vr);
				} else {
					if (vr.closestd > (1 - mu) * dx) {
						vr.closestd = (1 - mu) * dx;
						vr.i = i + 1;
						vr.j = j;

						vr.originalvalue = values[i + 1][j];
						vr.p = interp(isolevel, p0, p1, val0, val1).addSelf(offset);
						vr.snapvertex = vertex(i + 1, j, offset);
					}
				}
			}
			return null;
		}
		WB_Point xedge = xedges.get(index(i, j));
		if (xedge != null) {
			return xedge;
		}
		xedge = new WB_Point(interp(isolevel, p0, p1, val0, val1));
		xedge.addSelf(offset);
		xedges.put(index(i, j), xedge);
		return xedge;
	}

	/**
	 * Yedge.
	 *
	 * @param i
	 *            i: -1 .. resx+1
	 * @param j
	 *            j: -1 .. resy+1
	 * @param offset
	 * @param dummyrun
	 * @return edge vertex
	 */
	private WB_Point yedge(final int i, final int j, final WB_Point offset, final boolean dummyrun) {
		WB_Point yedge = yedges.get(index(i, j));
		if (yedge != null) {
			return yedge;
		}
		final WB_Point p0 = new WB_Point(i * dx, j * dy);
		final WB_Point p1 = new WB_Point(i * dx, j * dy + dy);
		final double val0 = value(i, j);
		final double val1 = value(i, j + 1);
		double mu;
		if (dummyrun) {
			mu = (isolevel - val0) / (val1 - val0);
			if (mu < gamma) {
				VertexRemap vr = vertexremaps.get(index(i, j));
				if (vr == null) {
					vr = new VertexRemap();
					vr.closestd = mu * dy;
					vr.i = i;
					vr.j = j;

					vr.originalvalue = values[i][j];
					vr.p = interp(isolevel, p0, p1, val0, val1).addSelf(offset);
					vr.snapvertex = vertex(i, j, offset);
					vertexremaps.put(index(i, j), vr);
				} else {
					if (vr.closestd > mu * dy) {
						vr.closestd = mu * dy;
						vr.i = i;
						vr.j = j;

						vr.originalvalue = values[i][j];
						vr.p = interp(isolevel, p0, p1, val0, val1).addSelf(offset);
						vr.snapvertex = vertex(i, j, offset);
					}
				}
			} else if (mu > 1 - gamma) {
				VertexRemap vr = vertexremaps.get(index(i, j + 1));
				if (vr == null) {
					vr = new VertexRemap();
					vr.closestd = (1 - mu) * dy;
					vr.i = i;
					vr.j = j + 1;

					vr.originalvalue = values[i][j + 1];
					vr.p = interp(isolevel, p0, p1, val0, val1).addSelf(offset);
					vr.snapvertex = vertex(i, j + 1, offset);
					vertexremaps.put(index(i, j + 1), vr);
				} else {
					if (vr.closestd > (1 - mu) * dy) {
						vr.closestd = (1 - mu) * dy;
						vr.i = i;
						vr.j = j + 1;

						vr.originalvalue = values[i][j + 1];
						vr.p = interp(isolevel, p0, p1, val0, val1).addSelf(offset);
						vr.snapvertex = vertex(i, j + 1, offset);
					}
				}
			}
			return null;
		}
		yedge = new WB_Point(interp(isolevel, p0, p1, val0, val1));
		yedge.addSelf(offset);
		yedges.put(index(i, j), yedge);
		return yedge;
	}

	/*
	 * Linearly interpolate the position where an isosurface cuts an edge
	 * between two vertices, each with their own scalar value
	 */
	/**
	 * Interp.
	 *
	 * @param isolevel
	 *            the isolevel
	 * @param p1
	 *            the p1
	 * @param p2
	 *            the p2
	 * @param valp1
	 *            the valp1
	 * @param valp2
	 *            the valp2
	 * @return the h e_ vertex
	 */
	private WB_Point interp(final double isolevel, final WB_Point p1, final WB_Point p2, final double valp1,
			final double valp2) {
		double mu;
		if (WB_Epsilon.isEqualAbs(isolevel, valp1)) {
			return p1;
		}
		if (WB_Epsilon.isEqualAbs(isolevel, valp2)) {
			return p2;
		}
		if (WB_Epsilon.isEqualAbs(valp1, valp2)) {
			return p1;
		}
		mu = (isolevel - valp1) / (valp2 - valp1);
		return new WB_Point(p1.xd() + mu * (p2.xd() - p1.xd()), p1.yd() + mu * (p2.yd() - p1.yd()));
	}

	/**
	 *
	 * @param i
	 * @param j
	 * @return
	 */
	private int index(final int i, final int j) {
		return i + 1 + (resx + 2) * (j + 1);
	}

	/**
	 * Value.
	 *
	 * @param i
	 *            i
	 * @param j
	 *            j
	 * @return double
	 */
	private double value(final int i, final int j) {
		if (Double.isNaN(boundary)) { // if no boundary is set i,j,k should
			// always be between o and resx,rey,resz
			return values[i][j];
		}
		if (i < 0 || j < 0 || i > resx || j > resy) {
			return invert ? -boundary : boundary;
		}
		return values[i][j];
	}

	/**
	 *
	 */
	class VertexRemap {
		/**
		 *
		 */
		int i, j;
		/**
		 *
		 */
		double closestd;
		/**
		 *
		 */
		WB_Point p;
		/**
		 *
		 */
		double originalvalue;
		/**
		 *
		 */
		WB_Point snapvertex;
	}

}
