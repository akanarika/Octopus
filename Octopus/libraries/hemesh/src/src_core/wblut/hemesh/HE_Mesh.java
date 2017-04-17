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

import java.util.Iterator;
import java.util.List;

import gnu.trove.map.TLongIntMap;
import gnu.trove.map.hash.TLongIntHashMap;
import javolution.util.FastTable;
import wblut.geom.WB_AABB;
import wblut.geom.WB_Coord;
import wblut.geom.WB_Frame;
import wblut.geom.WB_GeometryOp3D;
import wblut.geom.WB_KDTree;
import wblut.geom.WB_KDTree.WB_KDEntry;
import wblut.geom.WB_Mesh;
import wblut.geom.WB_MeshCreator;
import wblut.geom.WB_Point;
import wblut.geom.WB_Polygon;
import wblut.geom.WB_Transform;
import wblut.geom.WB_Triangle;
import wblut.geom.WB_Vector;

/**
 * Half-edge mesh data structure.
 *
 * @author Frederik Vanhoutte (W:Blut)
 *
 */
public class HE_Mesh extends HE_MeshStructure {

	/**
	 * Instantiates a new HE_Mesh.
	 *
	 */
	public HE_Mesh() {
		super();
	}

	/**
	 * Constructor.
	 *
	 * @param creator
	 *            HE_Creator that generates this mesh
	 */
	public HE_Mesh(final HEC_Creator creator) {
		super();
		setNoCopy(creator.create());
	}

	/**
	 *
	 *
	 * @param mesh
	 */
	public HE_Mesh(final WB_Mesh mesh) {
		this(new HEC_FromMesh(mesh));
	}

	/**
	 *
	 *
	 * @param mesh
	 */
	public HE_Mesh(final WB_MeshCreator mesh) {
		this(new HEC_FromMesh(mesh.create()));
	}

	/**
	 *
	 *
	 * @param mesh
	 */
	public HE_Mesh(final HE_Mesh mesh) {
		super();
		set(mesh);
	}

	/**
	 * Deep copy of mesh.
	 *
	 * @return copy as new HE_Mesh
	 */
	public HE_Mesh copy() {
		return new HE_Mesh(new HEC_Copy(this));
	}

	/**
	 * Deep copy of mesh.
	 *
	 * @return copy as new HE_Mesh
	 */
	public HE_Mesh get() {
		return new HE_Mesh(new HEC_Copy(this));
	}

	/**
	 * Replace mesh with deep copy of target.
	 *
	 * @param target
	 *            HE_Mesh to be duplicated
	 */
	public void set(final HE_Mesh target) {
		final HE_Mesh result = target.copy();
		replaceVertices(result);
		replaceFaces(result);
		replaceHalfedges(result);
	}

	/**
	 * Replace mesh with shallow copy of target.
	 *
	 * @param target
	 *            HE_Mesh to be duplicated
	 */
	void setNoCopy(final HE_Mesh target) {
		replaceVertices(target);
		replaceFaces(target);
		replaceHalfedges(target);

	}

	/**
	 * Modify the mesh.
	 *
	 * @param modifier
	 *            HE_Modifier to apply
	 * @return self
	 */
	public HE_Mesh modify(final HEM_Modifier modifier) {
		modifier.apply(this);
		update();
		return this;
	}

	/**
	 * Subdivide the mesh.
	 *
	 * @param subdividor
	 *            HE_Subdividor to apply
	 * @return self
	 */
	public HE_Mesh subdivide(final HES_Subdividor subdividor) {
		subdividor.apply(this);
		update();
		return this;
	}

	/**
	 * Subdivide the mesh a number of times.
	 *
	 * @param subdividor
	 *            HE_Subdividor to apply
	 * @param rep
	 *            subdivision iterations. WARNING: higher values will lead to
	 *            unmanageable number of faces.
	 * @return self
	 */
	public HE_Mesh subdivide(final HES_Subdividor subdividor, final int rep) {

		for (int i = 0; i < rep; i++) {
			subdividor.apply(this);
			update();
		}
		return this;
	}

	/**
	 * Smooth.
	 */
	public void smooth() {

		subdivide(new HES_CatmullClark());
		update();
	}

	/**
	 *
	 *
	 * @param rep
	 */
	public void smooth(final int rep) {

		subdivide(new HES_CatmullClark(), rep);
		update();
	}

	/**
	 * Simplify.
	 *
	 * @param simplifier
	 *            the simplifier
	 * @return the h e_ mesh
	 */
	public HE_Mesh simplify(final HES_Simplifier simplifier) {
		simplifier.apply(this);
		update();
		return this;
	}

	/**
	 * Add all mesh elements to this mesh. No copies are made. Tries to join
	 * geometry.
	 *
	 * @param mesh
	 *            mesh to add
	 */
	public void fuse(final HE_Mesh mesh) {
		addVertices(mesh.getVerticesAsArray());
		addFaces(mesh.getFacesAsArray());
		addHalfedges(mesh.getHalfedgesAsArray());
		setNoCopy(new HE_Mesh(new HEC_FromPolygons().setPolygons(this.getPolygonList())));
	}

	/**
	 *
	 *
	 * @return
	 */
	public WB_Mesh toFacelistMesh() {
		return gf.createMesh(getVerticesAsCoord(), getFacesAsInt());
	}

	/**
	 * Gets the frame.
	 *
	 * @return the frame
	 */
	public WB_Frame getFrame() {
		final WB_Frame frame = new WB_Frame(getVerticesAsCoord());
		final TLongIntMap map = new TLongIntHashMap(10, 0.5f, -1L, -1);
		map.putAll(getVertexKeyToIndexMap());
		final Iterator<HE_Halfedge> eItr = eItr();
		HE_Halfedge e;
		while (eItr.hasNext()) {
			e = eItr.next();
			frame.addStrut(map.get(e.getVertex().key()), map.get(e.getEndVertex().key()));
		}
		return frame;
	}

	public HE_Mesh apply(final WB_Transform T) {

		return new HEC_Transform(this, T).create();
	}

	/**
	 *
	 * @param T
	 * @return
	 */
	public HE_Mesh applySelf(final WB_Transform T) {

		update();
		return modify(new HEM_Transform(T));
	}

	/**
	 * Apply transform to entire mesh.
	 *
	 * @param T
	 *            WB_Transform to apply
	 *
	 * @return self
	 */
	public HE_Mesh transformSelf(final WB_Transform T) {

		update();
		return modify(new HEM_Transform(T));
	}

	/**
	 * Create transformed copy of mesh.
	 *
	 * @param T
	 *            WB_Transform to apply
	 *
	 * @return copy
	 */
	public HE_Mesh transform(final WB_Transform T) {
		return copy().modify(new HEM_Transform(T));
	}

	/**
	 * Translate entire mesh.
	 *
	 * @param x
	 *
	 * @param y
	 *
	 * @param z
	 *
	 * @return self
	 */
	public HE_Mesh moveSelf(final double x, final double y, final double z) {

		update();
		final Iterator<HE_Vertex> vItr = vItr();
		while (vItr.hasNext()) {
			vItr.next().addSelf(x, y, z);
		}
		return this;
	}

	/**
	 * Create translated copy of mesh.
	 *
	 * @param x
	 *
	 * @param y
	 *
	 * @param z
	 *
	 * @return copy
	 */
	public HE_Mesh move(final double x, final double y, final double z) {
		HE_Mesh result = copy();
		final Iterator<HE_Vertex> vItr = result.vItr();
		while (vItr.hasNext()) {
			vItr.next().addSelf(x, y, z);
		}
		return result;
	}

	/**
	 * Translate entire mesh.
	 *
	 * @param v
	 *            the v
	 * @return self
	 */
	public HE_Mesh moveSelf(final WB_Coord v) {

		update();
		return moveSelf(v.xd(), v.yd(), v.zd());
	}

	/**
	 * Created translated copy of mesh.
	 *
	 * @param v
	 *
	 * @return copy
	 */
	public HE_Mesh move(final WB_Coord v) {
		return move(v.xd(), v.yd(), v.zd());
	}

	/**
	 * Translate entire mesh to given position.
	 *
	 * @param x
	 *
	 * @param y
	 *
	 * @param z
	 *
	 * @return self
	 */
	public HE_Mesh moveToSelf(final double x, final double y, final double z) {

		update();
		WB_Point center = getCenter();
		final Iterator<HE_Vertex> vItr = vItr();
		while (vItr.hasNext()) {
			vItr.next().addSelf(x - center.xd(), y - center.yd(), z - center.zd());
		}
		return this;
	}

	/**
	 * Create copy of mesh at given position.
	 *
	 * @param x
	 *
	 * @param y
	 *
	 * @param z
	 *
	 * @return copy
	 */
	public HE_Mesh moveTo(final double x, final double y, final double z) {
		HE_Mesh result = copy();
		WB_Point center = result.getCenter();
		final Iterator<HE_Vertex> vItr = result.vItr();
		while (vItr.hasNext()) {
			vItr.next().addSelf(x - center.xd(), y - center.yd(), z - center.zd());
		}
		return result;
	}

	/**
	 * Translate entire mesh to given position.
	 *
	 * @param v
	 *
	 * @return self
	 */
	public HE_Mesh moveToSelf(final WB_Coord v) {

		update();
		return moveToSelf(v.xd(), v.yd(), v.zd());
	}

	/**
	 * create copy of mesh at given position.
	 *
	 * @param v
	 *
	 * @return copy
	 */
	public HE_Mesh moveTo(final WB_Coord v) {
		return moveTo(v.xd(), v.yd(), v.zd());
	}

	/**
	 * Rotate entire mesh around an arbitrary axis defined by 2 points.
	 *
	 * @param angle
	 *            angle
	 * @param p1x
	 *            x-coordinate of first point on axis
	 * @param p1y
	 *            y-coordinate of first point on axis
	 * @param p1z
	 *            z-coordinate of first point on axis
	 * @param p2x
	 *            x-coordinate of second point on axis
	 * @param p2y
	 *            y-coordinate of second point on axis
	 * @param p2z
	 *            z-coordinate of second point on axis
	 * @return self
	 */
	public HE_Mesh rotateAboutAxis2PSelf(final double angle, final double p1x, final double p1y, final double p1z,
			final double p2x, final double p2y, final double p2z) {

		update();
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = vItr();
		final WB_Transform raa = new WB_Transform();
		raa.addRotateAboutAxis(angle, new WB_Point(p1x, p1y, p1z), new WB_Vector(p2x - p1x, p2y - p1y, p2z - p1z));
		while (vItr.hasNext()) {
			v = vItr.next();
			raa.applyAsPointSelf(v);
		}
		return this;
	}

	/**
	 * Create rotated copy of mesh around an arbitrary axis defined by 2 points.
	 *
	 * @param angle
	 *            angle
	 * @param p1x
	 *            x-coordinate of first point on axis
	 * @param p1y
	 *            y-coordinate of first point on axis
	 * @param p1z
	 *            z-coordinate of first point on axis
	 * @param p2x
	 *            x-coordinate of second point on axis
	 * @param p2y
	 *            y-coordinate of second point on axis
	 * @param p2z
	 *            z-coordinate of second point on axis
	 * @return copy
	 */
	public HE_Mesh rotateAboutAxis2P(final double angle, final double p1x, final double p1y, final double p1z,
			final double p2x, final double p2y, final double p2z) {
		HE_Mesh result = copy();
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = result.vItr();
		final WB_Transform raa = new WB_Transform();
		raa.addRotateAboutAxis(angle, new WB_Point(p1x, p1y, p1z), new WB_Vector(p2x - p1x, p2y - p1y, p2z - p1z));
		while (vItr.hasNext()) {
			v = vItr.next();
			raa.applyAsPointSelf(v);
		}
		return result;
	}

	/**
	 * Rotate entire mesh around an arbitrary axis defined by 2 points..
	 *
	 * @param angle
	 *            angle
	 * @param p1
	 *            first point on axis
	 * @param p2
	 *            second point on axis
	 * @return self
	 */
	public HE_Mesh rotateAboutAxis2PSelf(final double angle, final WB_Coord p1, final WB_Coord p2) {

		update();
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = vItr();
		final WB_Transform raa = new WB_Transform();
		raa.addRotateAboutAxis(angle, p1, new WB_Vector(p1, p2));
		while (vItr.hasNext()) {
			v = vItr.next();
			raa.applyAsPointSelf(v);
		}
		return this;
	}

	/**
	 * Create copy of mesh rotated around an arbitrary axis defined by 2 points.
	 *
	 * @param angle
	 *            angle
	 * @param p1
	 *            first point on axis
	 * @param p2
	 *            second point on axis
	 * @return copy
	 */
	public HE_Mesh rotateAboutAxis2P(final double angle, final WB_Coord p1, final WB_Coord p2) {
		HE_Mesh result = copy();
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = result.vItr();
		final WB_Transform raa = new WB_Transform();
		raa.addRotateAboutAxis(angle, p1, new WB_Vector(p1, p2));
		while (vItr.hasNext()) {
			v = vItr.next();
			raa.applyAsPointSelf(v);
		}
		return result;
	}

	/**
	 * Rotate entire mesh around an arbitrary axis defined by a point and a
	 * direction.
	 *
	 * @param angle
	 *            angle
	 * @param p
	 *            rotation point
	 * @param a
	 *            axis
	 * @return self
	 */
	public HE_Mesh rotateAboutAxisSelf(final double angle, final WB_Coord p, final WB_Coord a) {

		update();
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = vItr();
		final WB_Transform raa = new WB_Transform();
		raa.addRotateAboutAxis(angle, p, a);
		while (vItr.hasNext()) {
			v = vItr.next();
			raa.applyAsPointSelf(v);
		}
		return this;
	}

	/**
	 * Create copy of mesh rotated around an arbitrary axis defined by a point
	 * and a direction.
	 *
	 * @param angle
	 *            angle
	 * @param p
	 *            rotation point
	 * @param a
	 *            axis
	 * @return copy
	 */
	public HE_Mesh rotateAboutAxis(final double angle, final WB_Coord p, final WB_Coord a) {
		HE_Mesh result = copy();
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = result.vItr();
		final WB_Transform raa = new WB_Transform();
		raa.addRotateAboutAxis(angle, p, a);
		while (vItr.hasNext()) {
			v = vItr.next();
			raa.applyAsPointSelf(v);
		}
		return result;
	}

	/**
	 * Rotate entire mesh around an arbitrary axis defined by a point and a
	 * direction.
	 *
	 * @param angle
	 * @param px
	 * @param py
	 * @param pz
	 * @param ax
	 * @param ay
	 * @param az
	 * @return self
	 */
	public HE_Mesh rotateAboutAxisSelf(final double angle, final double px, final double py, final double pz,
			final double ax, final double ay, final double az) {

		update();
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = vItr();
		final WB_Transform raa = new WB_Transform();
		raa.addRotateAboutAxis(angle, new WB_Point(px, py, pz), new WB_Vector(ax, ay, az));
		while (vItr.hasNext()) {
			v = vItr.next();
			raa.applyAsPointSelf(v);
		}
		return this;
	}

	/**
	 * Rotate entire mesh around an arbitrary axis defined by a point and a
	 * direction.
	 *
	 * @param angle
	 * @param px
	 * @param py
	 * @param pz
	 * @param ax
	 * @param ay
	 * @param az
	 * @return copy
	 */
	public HE_Mesh rotateAboutAxis(final double angle, final double px, final double py, final double pz,
			final double ax, final double ay, final double az) {
		HE_Mesh result = copy();
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = result.vItr();
		final WB_Transform raa = new WB_Transform();
		raa.addRotateAboutAxis(angle, new WB_Point(px, py, pz), new WB_Vector(ax, ay, az));
		while (vItr.hasNext()) {
			v = vItr.next();
			raa.applyAsPointSelf(v);
		}
		return result;
	}

	/**
	 * Rotate entire mesh around an arbitrary axis in origin.
	 *
	 * @param angle
	 *            angle
	 * @param a
	 *            axis
	 * @return self
	 */
	public HE_Mesh rotateAboutOriginSelf(final double angle, final WB_Coord a) {

		update();
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = vItr();
		final WB_Transform raa = new WB_Transform();
		raa.addRotateAboutOrigin(angle, a);
		while (vItr.hasNext()) {
			v = vItr.next();
			raa.applyAsPointSelf(v);
		}
		return this;
	}

	/**
	 * Create copy of mesh rotate around an arbitrary axis in origin.
	 *
	 * @param angle
	 *            angle
	 * @param a
	 *            axis
	 * @return copy
	 */
	public HE_Mesh rotateAboutOrigin(final double angle, final WB_Coord a) {
		HE_Mesh result = copy();
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = result.vItr();
		final WB_Transform raa = new WB_Transform();
		raa.addRotateAboutOrigin(angle, a);
		while (vItr.hasNext()) {
			v = vItr.next();
			raa.applyAsPointSelf(v);
		}
		return result;
	}

	/**
	 * Rotate entire mesh around an arbitrary axis in origin.
	 *
	 *
	 * @param angle
	 * @param ax
	 * @param ay
	 * @param az
	 * @return
	 */
	public HE_Mesh rotateAboutOriginSelf(final double angle, final double ax, final double ay, final double az) {

		update();
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = vItr();
		final WB_Transform raa = new WB_Transform();
		raa.addRotateAboutOrigin(angle, new WB_Vector(ax, ay, az));
		while (vItr.hasNext()) {
			v = vItr.next();
			raa.applyAsPointSelf(v);
		}
		return this;
	}

	/**
	 * Create copy of mesh rotated around an arbitrary axis in origin.
	 *
	 *
	 * @param angle
	 * @param ax
	 * @param ay
	 * @param az
	 * @return copy
	 */
	public HE_Mesh rotateAboutOrigin(final double angle, final double ax, final double ay, final double az) {
		HE_Mesh result = copy();
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = result.vItr();
		final WB_Transform raa = new WB_Transform();
		raa.addRotateAboutOrigin(angle, new WB_Vector(ax, ay, az));
		while (vItr.hasNext()) {
			v = vItr.next();
			raa.applyAsPointSelf(v);
		}
		return result;
	}

	/**
	 * Rotate entire mesh around an arbitrary axis in center.
	 *
	 * @param angle
	 *            angle
	 * @param a
	 *            axis
	 * @return self
	 */
	public HE_Mesh rotateAboutCenterSelf(final double angle, final WB_Coord a) {

		update();
		return rotateAboutAxisSelf(angle, getCenter(), a);
	}

	/**
	 * Create copy of mesh rotated around an arbitrary axis in center.
	 *
	 * @param angle
	 *            angle
	 * @param a
	 *            axis
	 * @return self
	 */
	public HE_Mesh rotateAboutCenter(final double angle, final WB_Coord a) {

		return rotateAboutAxis(angle, getCenter(), a);
	}

	/**
	 * Rotate entire mesh around an arbitrary axis in center.
	 *
	 * @param angle
	 * @param ax
	 * @param ay
	 * @param az
	 * @return
	 */
	public HE_Mesh rotateAboutCenterSelf(final double angle, final double ax, final double ay, final double az) {

		update();
		return rotateAboutAxisSelf(angle, getCenter(), new WB_Vector(ax, ay, az));
	}

	/**
	 * Create copy of mesh rotated around an arbitrary axis in center.
	 *
	 * @param angle
	 * @param ax
	 * @param ay
	 * @param az
	 * @return copy
	 */
	public HE_Mesh rotateAboutCenter(final double angle, final double ax, final double ay, final double az) {
		return rotateAboutAxis(angle, getCenter(), new WB_Vector(ax, ay, az));
	}

	/**
	 * Scale entire mesh around center point.
	 *
	 * @param scaleFactorx
	 *            x-coordinate of scale factor
	 * @param scaleFactory
	 *            y-coordinate of scale factor
	 * @param scaleFactorz
	 *            z-coordinate of scale factor
	 * @param c
	 *            center
	 * @return self
	 */
	public HE_Mesh scaleSelf(final double scaleFactorx, final double scaleFactory, final double scaleFactorz,
			final WB_Coord c) {

		update();
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = vItr();
		while (vItr.hasNext()) {
			v = vItr.next();
			v.set(c.xd() + scaleFactorx * (v.xd() - c.xd()), c.yd() + scaleFactory * (v.yd() - c.yd()),
					c.zd() + scaleFactorz * (v.zd() - c.zd()));
		}

		return this;
	}

	/**
	 * Create copy of mesh scaled around center point.
	 *
	 * @param scaleFactorx
	 *            x-coordinate of scale factor
	 * @param scaleFactory
	 *            y-coordinate of scale factor
	 * @param scaleFactorz
	 *            z-coordinate of scale factor
	 * @param c
	 *            center
	 * @return copy
	 */
	public HE_Mesh scale(final double scaleFactorx, final double scaleFactory, final double scaleFactorz,
			final WB_Coord c) {
		HE_Mesh result = copy();
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = result.vItr();
		while (vItr.hasNext()) {
			v = vItr.next();
			v.set(c.xd() + scaleFactorx * (v.xd() - c.xd()), c.yd() + scaleFactory * (v.yd() - c.yd()),
					c.zd() + scaleFactorz * (v.zd() - c.zd()));
		}

		return result;
	}

	/**
	 * Scale entire mesh around center point.
	 *
	 * @param scaleFactor
	 *            scale
	 * @param c
	 *            center
	 * @return self
	 */
	public HE_Mesh scaleSelf(final double scaleFactor, final WB_Coord c) {

		update();
		return scaleSelf(scaleFactor, scaleFactor, scaleFactor, c);
	}

	/**
	 * Create copy of mesh scaled around center point.
	 *
	 * @param scaleFactor
	 *            scale
	 * @param c
	 *            center
	 * @return copy
	 */
	public HE_Mesh scale(final double scaleFactor, final WB_Coord c) {
		return scale(scaleFactor, scaleFactor, scaleFactor, c);
	}

	/**
	 * Scale entire mesh around bodycenter.
	 *
	 * @param scaleFactorx
	 *            x-coordinate of scale factor
	 * @param scaleFactory
	 *            y-coordinate of scale factor
	 * @param scaleFactorz
	 *            z-coordinate of scale factor
	 * @return self
	 */
	public HE_Mesh scaleSelf(final double scaleFactorx, final double scaleFactory, final double scaleFactorz) {

		update();
		WB_Point center = getCenter();
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = vItr();
		while (vItr.hasNext()) {
			v = vItr.next();
			v.set(center.xd() + scaleFactorx * (v.xd() - center.xd()),
					center.yd() + scaleFactory * (v.yd() - center.yd()),
					center.zd() + scaleFactorz * (v.zd() - center.zd()));
		}
		return this;
	}

	/**
	 * Create copy of mesh scaled around bodycenter.
	 *
	 * @param scaleFactorx
	 *            x-coordinate of scale factor
	 * @param scaleFactory
	 *            y-coordinate of scale factor
	 * @param scaleFactorz
	 *            z-coordinate of scale factor
	 * @return copy
	 */
	public HE_Mesh scale(final double scaleFactorx, final double scaleFactory, final double scaleFactorz) {
		HE_Mesh result = copy();
		WB_Point center = result.getCenter();
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = result.vItr();
		while (vItr.hasNext()) {
			v = vItr.next();
			v.set(center.xd() + scaleFactorx * (v.xd() - center.xd()),
					center.yd() + scaleFactory * (v.yd() - center.yd()),
					center.zd() + scaleFactorz * (v.zd() - center.zd()));
		}
		return result;
	}

	/**
	 * Scale entire mesh around bodycenter.
	 *
	 * @param scaleFactor
	 *            scale
	 * @return self
	 */
	public HE_Mesh scaleSelf(final double scaleFactor) {

		update();
		return scaleSelf(scaleFactor, scaleFactor, scaleFactor);
	}

	/**
	 * Create copy of mesh scaled around bodycenter.
	 *
	 * @param scaleFactor
	 *            scale
	 * @return copy
	 */
	public HE_Mesh scale(final double scaleFactor) {
		return scale(scaleFactor, scaleFactor, scaleFactor);
	}

	/**
	 * Fit in aabb.
	 *
	 * @param AABB
	 *
	 */
	public void fitInAABB(final WB_AABB AABB) {
		final WB_AABB self = getAABB();
		moveSelf(new WB_Vector(self.getMin(), AABB.getMin()));
		scaleSelf(AABB.getWidth() / self.getWidth(), AABB.getHeight() / self.getHeight(),
				AABB.getDepth() / self.getDepth(), new WB_Point(AABB.getMin()));
	}

	public void fitInAABB(final WB_AABB from, final WB_AABB to) {

		moveSelf(new WB_Vector(from.getMin(), to.getMin()));
		scaleSelf(to.getWidth() / from.getWidth(), to.getHeight() / from.getHeight(), to.getDepth() / from.getDepth(),
				new WB_Point(to.getMin()));
	}

	/**
	 * Fit in aabb constrained.
	 *
	 * @param AABB
	 *
	 * @return
	 */
	public double fitInAABBConstrained(final WB_AABB AABB) {
		final WB_AABB self = getAABB();
		moveSelf(new WB_Vector(self.getCenter(), AABB.getCenter()));
		double f = Math.min(AABB.getWidth() / self.getWidth(), AABB.getHeight() / self.getHeight());
		f = Math.min(f, AABB.getDepth() / self.getDepth());
		scaleSelf(f, new WB_Point(AABB.getCenter()));
		return f;
	}

	public double fitInAABBConstrained(final WB_AABB from, final WB_AABB to) {

		moveSelf(new WB_Vector(from.getCenter(), to.getCenter()));
		double f = Math.min(to.getWidth() / from.getWidth(), to.getHeight() / from.getHeight());
		f = Math.min(f, to.getDepth() / from.getDepth());
		scaleSelf(f, new WB_Point(to.getCenter()));
		return f;
	}

	/**
	 * Get the center (average of all vertex positions).
	 *
	 * @return the center
	 */

	public WB_Point getCenter() {
		final WB_Point c = new WB_Point(0, 0, 0);
		final Iterator<HE_Vertex> vItr = vItr();
		while (vItr.hasNext()) {
			c.addSelf(vItr.next());
		}
		c.divSelf(getNumberOfVertices());
		return c;
	}

	/**
	 * Delete face and remove all references.
	 *
	 * @param faces
	 *            faces to delete
	 */
	public void deleteFaces(final HE_Selection faces) {
		HE_Face f;
		final Iterator<HE_Face> fItr = faces.fItr();
		while (fItr.hasNext()) {
			f = fItr.next();
			remove(f);
		}
		cleanUnusedElementsByFace();
		capHalfedges();
	}

	/**
	 * Delete face and remove all references. Its halfedges remain and form a
	 * valid boundary loop.
	 *
	 * @param f
	 *            face to delete
	 */
	public void deleteFace(final HE_Face f) {
		HE_Halfedge he = f.getHalfedge();
		do {
			clearFace(he);
			he = he.getNextInFace();
		} while (he != f.getHalfedge());
		remove(f);
	}

	/**
	 * Delete face and remove all references. Its halfedges are removed, the
	 * boundary loop is unpaired.
	 *
	 * @param f
	 */

	public void cutFace(final HE_Face f) {
		HE_Halfedge he = f.getHalfedge();
		do {
			setHalfedge(he.getVertex(), he.getNextInVertex());

			he = he.getNextInFace();
		} while (he != f.getHalfedge());

		do {
			clearFace(he);
			clearPair(he);
			remove(he);
			he = he.getNextInFace();
		} while (he != f.getHalfedge());
		remove(f);
	}

	/**
	 * Delete edge. Adjacent faces are fused.
	 *
	 * @param e
	 *            edge to delete
	 * @return fused face (or null)
	 */
	public HE_Face deleteEdge(final HE_Halfedge e) {
		HE_Face f = null;
		final HE_Halfedge he1 = e.isEdge() ? e : e.getPair();
		final HE_Halfedge he2 = he1.getPair();
		final HE_Halfedge he1n = he1.getNextInFace();
		final HE_Halfedge he2n = he2.getNextInFace();
		final HE_Halfedge he1p = he1.getPrevInFace();
		final HE_Halfedge he2p = he2.getPrevInFace();
		HE_Vertex v = he1.getVertex();
		if (v.getHalfedge() == he1) {
			setHalfedge(v, he1.getNextInVertex());
		}
		v = he2.getVertex();
		if (v.getHalfedge() == he2) {
			setHalfedge(v, he2.getNextInVertex());
		}
		setNext(he1p, he2n);
		setNext(he2p, he1n);
		if (he1.getFace() != null && he2.getFace() != null) {
			f = new HE_Face();
			f.copyProperties(e.getPair().getFace());
			add(f);
			setHalfedge(f, he1p);
			HE_Halfedge he = he1p;
			do {
				setFace(he, f);
				he = he.getNextInFace();
			} while (he != he1p);
		}
		if (he1.getFace() != null) {
			remove(he1.getFace());
		}
		if (he2.getFace() != null) {
			remove(he2.getFace());
		}
		remove(he1);
		remove(he2);
		return f;
	}

	/**
	 *
	 *
	 * @return
	 */
	public List<WB_Triangle> getTriangles() {
		final List<WB_Triangle> result = new FastTable<WB_Triangle>();
		final HE_Mesh trimesh = this.copy();
		trimesh.triangulate();
		final Iterator<HE_Face> fItr = trimesh.fItr();
		HE_Face f;
		while (fItr.hasNext()) {
			f = fItr.next();
			result.add(gf.createTriangle(f.getHalfedge().getVertex(), f.getHalfedge().getNextInFace().getVertex(),
					f.getHalfedge().getPrevInFace().getVertex()));
		}
		return result;
	}

	/**
	 * Triangulate all concave faces.
	 *
	 * @return
	 */
	public HE_Selection triangulateConcaveFaces() {
		return HET_MeshOp.triangulateConcaveFaces(this);
	}

	/**
	 *
	 *
	 * @param sel
	 * @return
	 */
	public HE_Selection triangulateConcaveFaces(final List<HE_Face> sel) {
		return HET_MeshOp.triangulateConcaveFaces(this, sel);
	}

	/**
	 * Triangulate face if concave.
	 *
	 * @param key
	 *            key of face
	 * @return
	 */
	public HE_Selection triangulateConcaveFace(final long key) {
		return HET_MeshOp.triangulateConcaveFace(this, key);
	}

	/**
	 * Triangulate face if concave.
	 *
	 * @param face
	 *            key of face
	 * @return
	 */
	public HE_Selection triangulateConcaveFace(final HE_Face face) {
		return HET_MeshOp.triangulateConcaveFace(this, face);
	}

	/**
	 * Check consistency of datastructure.
	 *
	 * @return true or false
	 */
	public boolean validate() {
		return HET_Diagnosis.validate(this);
	}

	/**
	 * Fuse all coplanar faces connected to face. New face can be concave.
	 *
	 * @param face
	 *            starting face
	 * @param a
	 *            the a
	 * @return new face
	 */
	public HE_Face fuseCoplanarFace(final HE_Face face, final double a) {

		List<HE_Face> neighbors;
		FastTable<HE_Face> facesToCheck = new FastTable<HE_Face>();
		final FastTable<HE_Face> newFacesToCheck = new FastTable<HE_Face>();
		facesToCheck.add(face);
		final HE_Selection sel = new HE_Selection(this);
		sel.add(face);
		HE_Face f;
		HE_Face fn;
		int ni = -1;
		int nf = 0;
		double sa = Math.sin(a);
		sa *= sa;
		while (ni < nf) {
			newFacesToCheck.clear();
			for (int i = 0; i < facesToCheck.size(); i++) {
				f = facesToCheck.get(i);
				neighbors = f.getNeighborFaces();
				for (int j = 0; j < neighbors.size(); j++) {
					fn = neighbors.get(j);
					if (!sel.contains(fn)) {
						if (WB_Vector.isParallel(f.getFaceNormal(), fn.getFaceNormal(), sa)) {
							sel.add(fn);
							newFacesToCheck.add(fn);
						}
					}
				}
			}
			facesToCheck = newFacesToCheck;
			ni = nf;
			nf = sel.getNumberOfFaces();
		}
		if (sel.getNumberOfFaces() == 1) {
			return face;
		}
		final List<HE_Halfedge> halfedges = sel.getOuterHalfedgesInside();
		final HE_Face newFace = new HE_Face();
		add(newFace);
		newFace.copyProperties(sel.getFaceWithIndex(0));
		setHalfedge(newFace, halfedges.get(0));
		for (int i = 0; i < halfedges.size(); i++) {
			final HE_Halfedge hei = halfedges.get(i);
			final HE_Halfedge hep = halfedges.get(i).getPair();
			for (int j = 0; j < halfedges.size(); j++) {
				final HE_Halfedge hej = halfedges.get(j);
				if (i != j && hep.getVertex() == hej.getVertex()) {
					setNext(hei, hej);
				}
			}
			setFace(hei, newFace);
			setHalfedge(hei.getVertex(), hei);
		}
		removeFaces(sel.getFacesAsArray());
		cleanUnusedElementsByFace();
		return newFace;
	}

	/**
	 * Fuse all planar faces. Can lead to concave faces.
	 *
	 */
	public void fuseCoplanarFaces() {
		fuseCoplanarFaces(0);
	}

	/**
	 * Fuse all planar faces. Can lead to concave faces.
	 *
	 * @param a
	 *            the a
	 */
	public void fuseCoplanarFaces(final double a) {
		int ni;
		int no;
		do {
			ni = getNumberOfFaces();
			final List<HE_Face> faces = this.getFaces();
			for (int i = 0; i < faces.size(); i++) {
				final HE_Face f = faces.get(i);
				if (contains(f)) {
					fuseCoplanarFace(f, a);
				}
			}
			no = getNumberOfFaces();
		} while (no < ni);
	}

	/**
	 * Return a KD-tree containing all face centers.
	 *
	 * @return WB_KDTree
	 */
	public WB_KDTree<WB_Coord, Long> getFaceTree() {
		final WB_KDTree<WB_Coord, Long> tree = new WB_KDTree<WB_Coord, Long>();
		HE_Face f;
		final Iterator<HE_Face> fItr = fItr();
		while (fItr.hasNext()) {
			f = fItr.next();
			tree.add(f.getFaceCenter(), f.key());
		}
		return tree;
	}

	/**
	 * Return a KD-tree containing all vertices.
	 *
	 * @return WB_KDTree
	 */
	public WB_KDTree<WB_Coord, Long> getVertexTree() {
		final WB_KDTree<WB_Coord, Long> tree = new WB_KDTree<WB_Coord, Long>();
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = vItr();
		while (vItr.hasNext()) {
			v = vItr.next();
			tree.add(v, v.key());
		}
		return tree;
	}

	/**
	 * Return the closest vertex on the mesh.
	 *
	 * @param p
	 *            query point
	 * @param vertexTree
	 *            KD-tree from mesh (from vertexTree())
	 * @return HE_Vertex closest vertex
	 */
	public HE_Vertex getClosestVertex(final WB_Coord p, final WB_KDTree<WB_Coord, Long> vertexTree) {
		final WB_KDEntry<WB_Coord, Long>[] closestVertex = vertexTree.getNearestNeighbors(p, 1);
		if (closestVertex.length == 0) {
			return null;
		}
		return getVertexWithKey(closestVertex[0].value);
	}

	/**
	 * Return the closest point on the mesh.
	 *
	 * @param p
	 *            query point
	 * @param vertexTree
	 *            KD-tree from mesh (from vertexTree())
	 * @return WB_Coordinate closest point
	 */
	public WB_Coord getClosestPoint(final WB_Coord p, final WB_KDTree<WB_Coord, Long> vertexTree) {
		final WB_KDEntry<WB_Coord, Long>[] closestVertex = vertexTree.getNearestNeighbors(p, 1);
		final HE_Vertex v = getVertexWithKey(closestVertex[0].value);
		if (v == null) {
			return null;
		}
		final List<HE_Face> faces = v.getFaceStar();
		double d;
		double dmin = Double.POSITIVE_INFINITY;
		WB_Coord result = new WB_Point();
		for (int i = 0; i < faces.size(); i++) {
			final WB_Polygon poly = faces.get(i).toPolygon();
			final WB_Coord tmp = WB_GeometryOp3D.getClosestPoint3D(p, poly);
			d = WB_GeometryOp3D.getSqDistance3D(tmp, p);
			if (d < dmin) {
				dmin = d;
				result = tmp;
			}
		}
		return result;
	}

	/**
	 * Split the closest face in the query point.
	 *
	 * @param p
	 *            query point
	 * @param vertexTree
	 *            KD-tree from mesh (from vertexTree())
	 */
	public void addPointInClosestFace(final WB_Coord p, final WB_KDTree<WB_Coord, Long> vertexTree) {
		final WB_KDEntry<WB_Coord, Long>[] closestVertex = vertexTree.getNearestNeighbors(p, 1);
		final HE_Vertex v = getVertexWithKey(closestVertex[0].value);
		final List<HE_Face> faces = v.getFaceStar();
		double d;
		double dmin = Double.POSITIVE_INFINITY;
		HE_Face face = new HE_Face();
		for (int i = 0; i < faces.size(); i++) {
			final WB_Polygon poly = faces.get(i).toPolygon();
			final WB_Coord tmp = WB_GeometryOp3D.getClosestPoint3D(p, poly);
			d = WB_GeometryOp3D.getSqDistance3D(tmp, p);
			if (d < dmin) {
				dmin = d;
				face = faces.get(i);
				;
			}
		}
		final HE_Vertex nv = HEM_TriSplit.splitFaceTri(this, face, p).vItr().next();
		vertexTree.add(nv, nv.key());
	}

	/**
	 *
	 *
	 * @return
	 */
	public double getArea() {
		final Iterator<HE_Face> fItr = fItr();
		double A = 0.0;
		while (fItr.hasNext()) {
			A += fItr.next().getFaceArea();
		}
		return A;
	}

	/**
	 * Triangulate face.
	 *
	 * @param key
	 *            key of face
	 * @return
	 */
	public HE_Selection triangulate(final long key) {
		return triangulate(getFaceWithKey(key));
	}

	/**
	 *
	 *
	 * @param v
	 * @return
	 */
	public HE_Selection triangulateFaceStar(final HE_Vertex v) {
		return HET_MeshOp.triangulateFaceStar(this, v);
	}

	/**
	 *
	 *
	 * @param vertexkey
	 * @return
	 */
	public HE_Selection triangulateFaceStar(final long vertexkey) {
		return HET_MeshOp.triangulateFaceStar(this, vertexkey);
	}

	/**
	 *
	 *
	 * @param face
	 * @return
	 */
	public HE_Selection triangulate(final HE_Face face) {
		return HET_MeshOp.triangulate(this, face);
	}

	/**
	 * Triangulate all faces.
	 *
	 * @return
	 */
	public HE_Selection triangulate() {
		return HET_MeshOp.triangulate(this);
	}

	/**
	 * Triangulate.
	 *
	 * @param sel
	 *            the sel
	 * @return
	 */
	public HE_Selection triangulate(final HE_Selection sel) {
		return HET_MeshOp.triangulate(sel);
	}

	/**
	 * Clean.
	 */
	public void clean() {
		modify(new HEM_Clean());
	}

	/**
	 * Set vertex positions to values in a 2D array. If length of array is not
	 * the same as number of vertices, nothing happens.
	 *
	 * @param values
	 *            2D array of float. First index is number of vertices, second
	 *            index is 3 (x-,y- and z-coordinate)
	 */
	public void setVerticesFromFloat(final float[][] values) {
		if (values.length != getNumberOfVertices()) {
			return;
		}
		int i = 0;
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = vItr();
		while (vItr.hasNext()) {
			v = vItr.next();
			v.set(values[i][0], values[i][1], values[i][2]);
			i++;
		}
		update();
	}

	/**
	 * Set vertex positions to values in a 1D array. If length of array is not
	 * 3* number of vertices, nothing happens.
	 *
	 * @param values
	 *            1D array of float. 3 values, x,y, and z, per point
	 */
	public void setVerticesFromFloat(final float[] values) {
		if (values.length != 3 * getNumberOfVertices()) {
			return;
		}
		int i = 0;
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = vItr();
		while (vItr.hasNext()) {
			v = vItr.next();
			v.set(values[i], values[i + 1], values[i + 2]);
			i += 3;
		}
		update();
	}

	/**
	 * Set vertex positions to array of WB_Coord. If length of array is not the
	 * same as number of vertices, nothing happens.
	 *
	 * @param values
	 *            array of WB_Coord.
	 */
	public void setVerticesFromPoint(final WB_Coord[] values) {
		if (values.length != getNumberOfVertices()) {
			return;
		}
		int i = 0;
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = vItr();
		while (vItr.hasNext()) {
			v = vItr.next();
			v.set(values[i]);
			i++;
		}
		update();
	}

	/**
	 * Set vertex positions to List of WB_Coord. If the size of the List is not
	 * the same as number of vertices, nothing happens.
	 *
	 * @param values
	 *            array of WB_Coord.
	 */
	public void setVerticesFromPoint(final List<? extends WB_Coord> values) {
		if (values.size() != getNumberOfVertices()) {
			return;
		}
		int i = 0;
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = vItr();
		while (vItr.hasNext()) {
			v = vItr.next();
			v.set(values.get(i));
			i++;
		}
		update();
	}

	/**
	 * Set vertex positions to values in a 2D array. If length of array is not
	 * the same as number of vertices, nothing happens.
	 *
	 * @param values
	 *            2D array of double. First index is number of vertices, second
	 *            index is 3 (x-,y- and z-coordinate)
	 */
	public void setVerticesFromDouble(final double[][] values) {
		if (values.length != getNumberOfVertices()) {
			return;
		}
		int i = 0;
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = vItr();
		while (vItr.hasNext()) {
			v = vItr.next();
			v.set(values[i][0], values[i][1], values[i][2]);
			i++;
		}
		update();
	}

	/**
	 * Set vertex positions to values in a 1D array. If length of array is not
	 * 3* number of vertices, nothing happens.
	 *
	 * @param values
	 *            1D array of float. 3 values, x,y, and z, per point
	 */
	public void setVerticesFromFloat(final double[] values) {
		if (values.length != 3 * getNumberOfVertices()) {
			return;
		}
		int i = 0;
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = vItr();
		while (vItr.hasNext()) {
			v = vItr.next();
			v.set(values[i], values[i + 1], values[i + 2]);
			i += 3;
		}
		update();
	}

	/**
	 * Set vertex positions to values in a 2D array. If length of array is not
	 * the same as number of vertices, nothing happens.
	 *
	 * @param values
	 *            2D array of int. First index is number of vertices, second
	 *            index is 3 (x-,y- and z-coordinate)
	 */
	public void setVerticesFromInt(final int[][] values) {
		if (values.length != getNumberOfVertices()) {
			return;
		}
		int i = 0;
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = vItr();
		while (vItr.hasNext()) {
			v = vItr.next();
			v.set(values[i][0], values[i][1], values[i][2]);
			i++;
		}
		update();
	}

	/**
	 * Set vertex positions to values in a 1D array. If length of array is not
	 * 3* number of vertices, nothing happens.
	 *
	 * @param values
	 *            1D array of float. 3 values, x,y, and z, per point
	 */
	public void setVerticesFromFloat(final int[] values) {
		if (values.length != 3 * getNumberOfVertices()) {
			return;
		}
		int i = 0;
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = vItr();
		while (vItr.hasNext()) {
			v = vItr.next();
			v.set(values[i], values[i + 1], values[i + 2]);
			i += 3;
		}
		update();
	}

	public WB_Coord getFaceNormal(final int id) {
		return getFaceWithIndex(id).getFaceNormal();
	}

	public WB_Coord getFaceCenter(final int id) {
		return getFaceWithIndex(id).getFaceCenter();
	}

	public WB_Coord getVertexNormal(final int i) {
		return getVertexWithIndex(i).getVertexNormal();
	}

	public WB_Coord getVertex(final int i) {
		return getVertexWithIndex(i);
	}

	public List<WB_Coord> getPoints() {
		final List<WB_Coord> result = new FastTable<WB_Coord>();
		result.addAll(vertices.getObjects());
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = vItr();
		while (vItr.hasNext()) {
			v = vItr.next();
			result.add(v);

		}
		return result;
	}

}