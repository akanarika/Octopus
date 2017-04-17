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

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import javolution.util.FastTable;
import wblut.geom.WB_Coord;
import wblut.geom.WB_Point;
import wblut.geom.WB_Voronoi;

/**
 * Creates the Voronoi cells of a collection of points, constrained by a mesh.
 *
 * @author Frederik Vanhoutte (W:Blut)
 *
 */
public class HEMC_VoronoiCells extends HEMC_MultiCreator {
	/** Points. */
	private List<WB_Coord> points;
	/** Number of points. */
	private int numberOfPoints;
	/** Container. */
	private HE_Mesh container;
	/** Treat container as surface?. */
	private boolean surface;
	/** The simple cap. */
	private boolean simpleCap;
	/**
	 *
	 */
	private boolean bruteForce;
	/** Offset. */
	private double offset;
	/** The inner. */
	public HE_Selection[] inner;
	/** The outer. */
	public HE_Selection[] outer;
	/** Create divided skin of container. */
	private boolean createSkin;

	/**
	 * Instantiates a new HEMC_VoronoiCells.
	 *
	 */
	public HEMC_VoronoiCells() {
		super();
		simpleCap = true;
	}

	/**
	 * Set mesh, defines both points and container.
	 *
	 * @param mesh
	 *            HE_Mesh
	 * @param addCenter
	 *            add mesh center as extra point?
	 * @return self
	 */
	public HEMC_VoronoiCells setMesh(final HE_Mesh mesh, final boolean addCenter) {
		if (addCenter) {
			points = new FastTable<WB_Coord>();
			points.addAll(mesh.getVertices());
			points.add(mesh.getCenter());
		} else {
			points = new FastTable<WB_Coord>();
			points.addAll(mesh.getVertices());
		}
		container = mesh;
		return this;
	}

	/**
	 * Set points that define cell centers.
	 *
	 * @param points
	 *            array of vertex positions
	 * @return self
	 */
	public HEMC_VoronoiCells setPoints(final WB_Coord[] points) {
		this.points = new FastTable<WB_Coord>();
		for (WB_Coord p : points) {
			this.points.add(p);
		}
		return this;
	}

	/**
	 * Set points that define cell centers.
	 *
	 * @param points
	 *            collection of vertex positions
	 * @return self
	 */
	public HEMC_VoronoiCells setPoints(final Collection<? extends WB_Coord> points) {
		this.points = new FastTable<WB_Coord>();
		this.points.addAll(points);
		return this;
	}

	/**
	 * Set points that define cell centers.
	 *
	 * @param points
	 *            2D array of double of vertex positions
	 * @return self
	 */
	public HEMC_VoronoiCells setPoints(final double[][] points) {
		final int n = points.length;
		this.points = new FastTable<WB_Coord>();
		for (int i = 0; i < n; i++) {
			this.points.add(new WB_Point(points[i][0], points[i][1], points[i][2]));
		}
		return this;
	}

	/**
	 * Set points that define cell centers.
	 *
	 * @param points
	 *            2D array of float of vertex positions
	 * @return self
	 */
	public HEMC_VoronoiCells setPoints(final float[][] points) {
		final int n = points.length;
		this.points = new FastTable<WB_Coord>();
		for (int i = 0; i < n; i++) {
			this.points.add(new WB_Point(points[i][0], points[i][1], points[i][2]));
		}
		return this;
	}

	/**
	 * Set number of points.
	 *
	 * @param N
	 *            number of points
	 * @return self
	 */
	public HEMC_VoronoiCells setN(final int N) {
		numberOfPoints = N;
		return this;
	}

	/**
	 * Set voronoi cell offset.
	 *
	 * @param o
	 *            offset
	 * @return self
	 */
	public HEMC_VoronoiCells setOffset(final double o) {
		offset = o;
		return this;
	}

	/**
	 * Set enclosing mesh limiting cells.
	 *
	 * @param container
	 *            enclosing mesh
	 * @return self
	 */
	public HEMC_VoronoiCells setContainer(final HE_Mesh container) {
		this.container = container;
		return this;
	}

	/**
	 * Set optional surface mesh mode.
	 *
	 * @param b
	 *            true, false
	 * @return self
	 */
	public HEMC_VoronoiCells setSurface(final boolean b) {
		surface = b;
		return this;
	}

	/**
	 * Sets the simple cap.
	 *
	 * @param b
	 *            the b
	 * @return the hEM c_ voronoi cells
	 */
	public HEMC_VoronoiCells setSimpleCap(final boolean b) {
		simpleCap = b;
		return this;
	}

	/**
	 * Create skin mesh?.
	 *
	 * @param b
	 *            true, false
	 * @return self
	 */
	public HEMC_VoronoiCells setCreateSkin(final boolean b) {
		createSkin = b;
		return this;
	}

	/**
	 *
	 *
	 * @param b
	 * @return
	 */
	public HEMC_VoronoiCells setBruteForce(final boolean b) {
		bruteForce = b;
		return this;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.hemesh.HE_MultiCreator#create()
	 */
	@Override
	public HE_MeshCollection create() {

		tracker.setStatus(this, "Starting HEMC_VoronoiCells", +1);
		HE_MeshCollection result = new HE_MeshCollection();
		if (container == null) {
			_numberOfMeshes = 0;
			return result;
		}
		if (points == null) {
			result.add(container.copy());
			_numberOfMeshes = 1;
			return result;
		}
		if (numberOfPoints == 0) {
			numberOfPoints = points.size();
		}
		final ArrayList<HE_Selection> linnersel = new ArrayList<HE_Selection>();
		final ArrayList<HE_Selection> loutersel = new ArrayList<HE_Selection>();
		final HEC_VoronoiCell cvc = new HEC_VoronoiCell();
		if (bruteForce || numberOfPoints < 10) {
			cvc.setPoints(points).setN(numberOfPoints).setContainer(container).setSurface(surface).setOffset(offset)
					.setSimpleCap(simpleCap);
			for (int i = 0; i < numberOfPoints; i++) {
				tracker.setStatus(this, "Creating cell " + i + " (" + numberOfPoints + " slices).", 0);
				cvc.setCellIndex(i);
				final HE_Mesh mesh = cvc.createBase();
				linnersel.add(cvc.inner);
				loutersel.add(cvc.outer);
				result.add(mesh);
			}
		} else {
			final int[][] voronoiIndices = WB_Voronoi.getVoronoi3DNeighbors(points);
			cvc.setPoints(points).setN(numberOfPoints).setContainer(container).setSurface(surface).setOffset(offset)
					.setSimpleCap(simpleCap).setLimitPoints(true);
			for (int i = 0; i < numberOfPoints; i++) {
				tracker.setStatus(this, "Creating cell " + i + " (" + voronoiIndices[i].length + " slices).", 0);
				cvc.setCellIndex(i);
				cvc.setPointsToUse(voronoiIndices[i]);

				final HE_Mesh mesh = cvc.createBase();
				linnersel.add(cvc.inner);
				loutersel.add(cvc.outer);
				result.add(mesh);
			}
		}
		inner = new HE_Selection[result.size()];
		outer = new HE_Selection[result.size()];

		for (int i = 0; i < _numberOfMeshes; i++) {

			inner[i] = linnersel.get(i);
			outer[i] = loutersel.get(i);
		}
		if (createSkin) {
			tracker.setStatus(this, "Creating skin.", 0);
			final boolean[] on = new boolean[_numberOfMeshes];
			for (int i = 0; i < _numberOfMeshes; i++) {
				on[i] = true;
			}
			result.add(new HE_Mesh(new HEC_FromVoronoiCells().setActive(on).setCells(result.meshes)));
		}
		_numberOfMeshes = result.size();
		tracker.setStatus(this, "Exiting HEMC_VoronoiCells.", -1);

		return result;
	}

	public static void main(final String[] args) {
		HEC_Torus creator = new HEC_Torus(80, 200, 16, 64);
		HE_Mesh container = new HE_Mesh(creator);
		creator = new HEC_Torus(40, 200, 16, 64);
		HE_Mesh inner = new HE_Mesh(creator);
		HET_MeshOp.flipFaces(inner);
		container.add(inner);

		int numpoints = 50;
		double[][] points = new double[numpoints][3];
		for (int i = 0; i < numpoints; i++) {
			points[i][0] = Math.random() * 200;
			points[i][1] = Math.random() * 200;
			points[i][2] = Math.random() * 200;
		}

		HEMC_VoronoiCells multiCreator = new HEMC_VoronoiCells();
		multiCreator.setPoints(points);
		multiCreator.setN(numpoints);
		multiCreator.setContainer(container);
		multiCreator.setOffset(0);
		multiCreator.setSimpleCap(false);
		multiCreator.create();
	}
}
