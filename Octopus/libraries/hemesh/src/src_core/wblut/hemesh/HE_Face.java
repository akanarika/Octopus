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

import static wblut.geom.WB_GeometryOp3D.projectOnPlane;

import java.util.List;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Polygon;
import com.vividsolutions.jts.operation.valid.IsValidOp;

import javolution.util.FastTable;
import wblut.geom.WB_AABB;
import wblut.geom.WB_Classification;
import wblut.geom.WB_Coord;
import wblut.geom.WB_GeometryOp3D;
import wblut.geom.WB_Map2D;
import wblut.geom.WB_Plane;
import wblut.geom.WB_Point;
import wblut.geom.WB_Polygon;
import wblut.geom.WB_PolygonTriangulator;
import wblut.geom.WB_Triangle;
import wblut.geom.WB_Vector;
import wblut.math.WB_Epsilon;

/**
 * Face element of half-edge data structure.
 *
 * @author Frederik Vanhoutte (W:Blut)
 *
 */
public class HE_Face extends HE_MeshElement implements Comparable<HE_Face> {
	/** Halfedge associated with this face. */
	private HE_Halfedge _halfedge;
	private int textureId;
	int[] triangles;
	WB_Coord normal;
	WB_Point center;

	/**
	 * Instantiates a new HE_Face.
	 */
	public HE_Face() {
		super();
		triangles = null;
		normal = null;
		center = null;
	}

	/**
	 *
	 *
	 * @return
	 */
	public HE_FaceEdgeCirculator feCrc() {
		return new HE_FaceEdgeCirculator(this);
	}

	/**
	 *
	 *
	 * @return
	 */
	public HE_FaceFaceCirculator ffCrc() {
		return new HE_FaceFaceCirculator(this);
	}

	/**
	 *
	 *
	 * @return
	 */
	public HE_FaceVertexCirculator fvCrc() {
		return new HE_FaceVertexCirculator(this);
	}

	/**
	 *
	 *
	 * @return
	 */
	public HE_FaceHalfedgeInnerCirculator fheiCrc() {
		return new HE_FaceHalfedgeInnerCirculator(this);
	}

	/**
	 *
	 *
	 * @return
	 */
	public HE_FaceHalfedgeOuterCirculator fheoCrc() {
		return new HE_FaceHalfedgeOuterCirculator(this);
	}

	/**
	 *
	 *
	 * @return
	 */
	public HE_FaceEdgeRevCirculator feRevCrc() {
		return new HE_FaceEdgeRevCirculator(this);
	}

	/**
	 *
	 *
	 * @return
	 */
	public HE_FaceFaceRevCirculator ffRevCrc() {
		return new HE_FaceFaceRevCirculator(this);
	}

	/**
	 *
	 *
	 * @return
	 */
	public HE_FaceHalfedgeInnerRevCirculator fheiRevCrc() {
		return new HE_FaceHalfedgeInnerRevCirculator(this);
	}

	/**
	 *
	 *
	 * @return
	 */
	public HE_FaceHalfedgeOuterRevCirculator fheoRevCrc() {
		return new HE_FaceHalfedgeOuterRevCirculator(this);
	}

	/**
	 *
	 *
	 * @return
	 */
	public HE_FaceVertexRevCirculator fvRevCrc() {
		return new HE_FaceVertexRevCirculator(this);
	}

	/**
	 *
	 *
	 * @return
	 */
	public long key() {
		return super.getKey();
	}

	/**
	 *
	 *
	 * @return
	 */
	public WB_Coord getFaceCenter() {
		if (center != null) {
			return center;
		}
		if (_halfedge == null) {
			return null;
		}
		HE_Halfedge he = _halfedge;
		center = new WB_Point();
		int c = 0;
		do {
			center.addSelf(he.getVertex());
			c++;
			he = he.getNextInFace();
		} while (he != _halfedge);
		center.divSelf(c);
		return center;
	}

	/**
	 *
	 *
	 * @param d
	 * @return
	 */
	public WB_Coord getFaceCenter(final double d) {
		if (center != null) {
			return center.addMul(d, getFaceNormal());
		}
		if (_halfedge == null) {
			return null;
		}
		HE_Halfedge he = _halfedge;
		center = new WB_Point();
		int c = 0;
		do {
			center.addSelf(he.getVertex());
			c++;
			he = he.getNextInFace();
		} while (he != _halfedge);
		center.divSelf(c);
		return center.addMul(d, getFaceNormal());
	}

	/**
	 *
	 *
	 * @return
	 */
	public WB_Coord getFaceNormal() {
		if (normal == null) {
			normal = HET_MeshOp.getFaceNormal(this);
		}
		return normal;
	}

	/**
	 *
	 *
	 * @return
	 */
	public WB_Coord getNonNormFaceNormal() {
		return HET_MeshOp.getNonNormFaceNormal(this);
	}

	/**
	 *
	 *
	 * @return
	 */
	public double getFaceArea() {
		return HET_MeshOp.getFaceArea(this);
	}

	/**
	 *
	 *
	 * @return
	 */
	public WB_Classification getFaceType() {
		return HET_MeshOp.getFaceType(this);
	}

	/**
	 *
	 *
	 * @return
	 */
	public List<HE_Vertex> getUniqueFaceVertices() {
		final FastTable<HE_Vertex> fv = new FastTable<HE_Vertex>();
		if (_halfedge == null) {
			return fv;
		}
		HE_Halfedge he = _halfedge;
		do {
			if (!fv.contains(he.getVertex())) {
				fv.add(he.getVertex());
			}
			he = he.getNextInFace();
		} while (he != _halfedge);
		return fv.unmodifiable();
	}

	public List<HE_Vertex> getFaceVertices() {
		final FastTable<HE_Vertex> fv = new FastTable<HE_Vertex>();
		if (_halfedge == null) {
			return fv;
		}
		HE_Halfedge he = _halfedge;
		do {

			fv.add(he.getVertex());

			he = he.getNextInFace();
		} while (he != _halfedge);
		return fv.unmodifiable();
	}

	/**
	 *
	 *
	 * @return
	 */
	public List<HE_TextureCoordinate> getFaceUVWs() {
		final FastTable<HE_TextureCoordinate> fv = new FastTable<HE_TextureCoordinate>();
		if (_halfedge == null) {
			return fv;
		}
		HE_Halfedge he = _halfedge;
		do {
			if (!fv.contains(he.getVertex())) {
				fv.add(he.getUVW());
			}
			he = he.getNextInFace();
		} while (he != _halfedge);
		return fv.unmodifiable();
	}

	/**
	 *
	 *
	 * @return
	 */
	public int getFaceOrder() {
		int result = 0;
		if (_halfedge == null) {
			return 0;
		}
		HE_Halfedge he = _halfedge;
		do {
			result++;
			he = he.getNextInFace();
		} while (he != _halfedge);
		return result;
	}

	/**
	 *
	 *
	 * @return
	 */
	public List<HE_Halfedge> getFaceHalfedges() {
		final FastTable<HE_Halfedge> fhe = new FastTable<HE_Halfedge>();
		if (_halfedge == null) {
			return fhe;
		}
		HE_Halfedge he = _halfedge;
		do {
			if (!fhe.contains(he)) {
				fhe.add(he);
			}
			he = he.getNextInFace();
		} while (he != _halfedge);
		return fhe.unmodifiable();
	}

	/**
	 *
	 *
	 * @return
	 */
	public List<HE_Halfedge> getFaceHalfedgesTwoSided() {
		final FastTable<HE_Halfedge> fhe = new FastTable<HE_Halfedge>();
		if (_halfedge == null) {
			return fhe;
		}
		HE_Halfedge he = _halfedge;
		do {
			if (!fhe.contains(he)) {
				fhe.add(he);
				if (he.getPair() != null) {
					if (!fhe.contains(he.getPair())) {
						fhe.add(he.getPair());
					}
				}
			}
			he = he.getNextInFace();
		} while (he != _halfedge);
		return fhe.unmodifiable();
	}

	/**
	 *
	 *
	 * @return
	 */
	public List<HE_Halfedge> getFaceEdges() {
		final FastTable<HE_Halfedge> fe = new FastTable<HE_Halfedge>();
		if (_halfedge == null) {
			return fe;
		}
		HE_Halfedge he = _halfedge;
		do {
			if (he.isEdge()) {
				if (!fe.contains(he)) {
					fe.add(he);
				}
			} else {
				if (!fe.contains(he.getPair())) {
					fe.add(he.getPair());
				}
			}
			he = he.getNextInFace();
		} while (he != _halfedge);
		return fe.unmodifiable();
	}

	/**
	 *
	 *
	 * @return
	 */
	public HE_Halfedge getHalfedge() {
		return _halfedge;
	}

	/**
	 *
	 *
	 * @param v
	 * @return
	 */
	public HE_Halfedge getHalfedge(final HE_Vertex v) {
		HE_Halfedge he = _halfedge;
		if (he == null) {
			return null;
		}
		do {
			if (he.getVertex() == v) {
				return he;
			}
			he = he.getNextInFace();
		} while (he != _halfedge);
		return null;
	}

	/**
	 *
	 *
	 * @param halfedge
	 */
	protected void _setHalfedge(final HE_Halfedge halfedge) {
		_halfedge = halfedge;
	}

	/**
	 *
	 *
	 * @param c
	 */
	public void push(final WB_Coord c) {
		HE_Halfedge he = _halfedge;
		do {
			he.getVertex().addSelf(c);
			he = he.getNextInFace();
		} while (he != _halfedge);
	}

	/**
	 *
	 */
	protected void _clearHalfedge() {
		_halfedge = null;
	}

	/**
	 *
	 *
	 * @return
	 * @deprecated Use {@link #getPlane()} instead
	 */
	@Deprecated
	public WB_Plane toPlane() {
		return getPlane();
	}

	/**
	 *
	 *
	 * @return
	 */
	public WB_Plane getPlane() {
		WB_Coord fn = getFaceNormal();
		if (WB_Vector.getSqLength3D(fn) < 0.5) {
			if (WB_Epsilon.isEqualAbs(_halfedge.getVertex().xd(), _halfedge.getEndVertex().xd())) {
				fn = new WB_Vector(1, 0, 0);
			} else {
				fn = new WB_Vector(0, 0, 1);
			}
		}
		return new WB_Plane(getFaceCenter(), fn);
	}

	/**
	 *
	 *
	 * @param d
	 * @return
	 */
	public WB_Plane getPlane(final double d) {
		final WB_Coord fn = getFaceNormal();
		return new WB_Plane(WB_Point.addMul(getFaceCenter(), d, fn), fn);
	}

	/**
	 *
	 */
	public void sort() {
		if (_halfedge != null) {
			HE_Halfedge he = _halfedge;
			HE_Halfedge leftmost = he;
			do {
				he = he.getNextInFace();
				if (he.getVertex().compareTo(leftmost.getVertex()) < 0) {
					leftmost = he;
				}
			} while (he != _halfedge);
			_halfedge = leftmost;
		}
	}

	/**
	 *
	 *
	 * @param f
	 * @return
	 */
	@Override
	public int compareTo(final HE_Face f) {

		if (f.getHalfedge() == null) {
			if (getHalfedge() == null) {

				return 0;
			} else {
				return 1;
			}
		} else if (getHalfedge() == null) {
			return -1;
		}

		return getHalfedge().compareTo(f.getHalfedge());
	}

	/**
	 *
	 *
	 * @return
	 */
	public int[] getTriangles() {
		return getTriangles(true);
	}

	/**
	 *
	 *
	 * @param optimize
	 * @return
	 */
	public int[] getTriangles(final boolean optimize) {
		if (triangles != null) {
			return triangles;
		}
		final int fo = getFaceOrder();
		if (fo < 3) {
			return new int[] { 0, 0, 0 };
		} else if (fo == 3) {
			return new int[] { 0, 1, 2 };
		} else if (isDegenerate()) {
			triangles = new int[3 * (fo - 2)];
			for (int i = 0; i < fo - 2; i++) {
				triangles[3 * i] = 0;
				triangles[3 * i + 1] = i + 1;
				triangles[3 * i + 2] = i + 2;
			}
		} else if (fo == 4) {
			final WB_Point[] points = new WB_Point[4];
			int i = 0;
			HE_Halfedge he = _halfedge;
			do {
				points[i] = new WB_Point(he.getVertex().xd(), he.getVertex().yd(), he.getVertex().zd());
				he = he.getNextInFace();
				i++;
			} while (he != _halfedge);

			return triangles = WB_PolygonTriangulator.triangulateQuad(points[0], points[1], points[2], points[3]);
		}
		return triangles = new WB_PolygonTriangulator().triangulateFace(this, optimize);
	}

	/**
	 *
	 *
	 * @return
	 */
	public WB_AABB toAABB() {
		final WB_AABB aabb = new WB_AABB();
		HE_Halfedge he = getHalfedge();
		do {
			aabb.expandToInclude(he.getVertex());
			he = he.getNextInFace();
		} while (he != getHalfedge());
		return aabb;
	}

	/**
	 *
	 *
	 * @return
	 */
	public WB_Triangle toTriangle() {
		if (getFaceOrder() != 3) {
			return null;
		}
		return new WB_Triangle(_halfedge.getVertex(), _halfedge.getEndVertex(),
				_halfedge.getNextInFace().getEndVertex());
	}

	/**
	 *
	 *
	 * @return
	 */
	public WB_Polygon toPolygon() {
		final int n = getFaceOrder();
		if (n == 0) {
			return null;
		}
		final WB_Point[] points = new WB_Point[n];
		int i = 0;
		HE_Halfedge he = _halfedge;
		do {
			points[i++] = new WB_Point(he.getVertex().xd(), he.getVertex().yd(), he.getVertex().zd());
			he = he.getNextInFace();
		} while (he != _halfedge);
		return gf.createSimplePolygon(points);
	}

	/**
	 *
	 *
	 * @return
	 */
	public WB_Polygon toOrthoPolygon() {
		final int n = getFaceOrder();
		if (n == 0) {
			return null;
		}
		final WB_Point[] points = new WB_Point[n];
		int i = 0;
		HE_Halfedge he = _halfedge;
		do {
			points[i++] = new WB_Point(he.getVertex().xd(), he.getVertex().yd(), he.getVertex().zd());
			he = he.getNextInFace();
		} while (he != _halfedge);
		return gf.createSimplePolygon(points);
	}

	/**
	 *
	 *
	 * @return
	 */
	public WB_Polygon toPlanarPolygon() {
		final int n = getFaceOrder();
		if (n == 0) {
			return null;
		}
		final WB_Point[] points = new WB_Point[n];
		final WB_Plane P = getPlane();
		int i = 0;
		HE_Halfedge he = _halfedge;
		do {
			points[i] = projectOnPlane(he.getVertex(), P);
			he = he.getNextInFace();
			i++;
		} while (he != _halfedge);
		return gf.createSimplePolygon(points);
	}

	/**
	 *
	 *
	 * @return
	 */
	public List<HE_Face> getNeighborFaces() {
		final FastTable<HE_Face> ff = new FastTable<HE_Face>();
		if (getHalfedge() == null) {
			return ff;
		}
		HE_Halfedge he = getHalfedge();
		do {
			final HE_Halfedge hep = he.getPair();
			if (hep != null && hep.getFace() != null) {
				if (hep.getFace() != this) {
					if (!ff.contains(hep.getFace())) {
						ff.add(hep.getFace());
					}
				}
			}
			he = he.getNextInFace();
		} while (he != getHalfedge());
		return ff.unmodifiable();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.geom.Point3D#toString()
	 */
	@Override
	public String toString() {
		String s = "HE_Face key: " + key() + ". Connects " + getFaceOrder() + " vertices: ";
		HE_Halfedge he = getHalfedge();
		for (int i = 0; i < getFaceOrder() - 1; i++) {
			s += he.getVertex().key + "-";
			he = he.getNextInFace();
		}
		s += he.getVertex().key + "." + " (" + getLabel() + "," + getInternalLabel() + ")";
		return s;
	}

	/**
	 *
	 *
	 * @return
	 */
	public boolean isPlanar() {
		final WB_Plane P = getPlane();
		HE_Halfedge he = getHalfedge();
		do {
			if (!WB_Epsilon.isZero(WB_GeometryOp3D.getDistance3D(he.getVertex(), P))) {
				return false;
			}
			he = he.getNextInFace();
		} while (he != getHalfedge());
		return true;
	}

	/**
	 * Checks if is boundary.
	 *
	 * @return true, if is boundary
	 */
	public boolean isBoundary() {
		HE_Halfedge he = _halfedge;
		do {
			if (he.getPair().getFace() == null) {
				return true;
			}
			he = he.getNextInFace();
		} while (he != _halfedge);
		return false;
	}

	/**
	 *
	 *
	 * @return
	 */
	public boolean isDegenerate() {
		return WB_Vector.getLength3D(getFaceNormal()) < 0.5;
	}

	/**
	 *
	 *
	 * @param el
	 */
	public void copyProperties(final HE_Face el) {
		super.copyProperties(el);
		textureId = el.textureId;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.hemesh.HE_Element#clear()
	 */
	@Override
	public void clear() {
		_halfedge = null;
	}

	/**
	 *
	 */
	public void checkValidity() {
		final Coordinate[] coords = new Coordinate[getFaceOrder() + 1];
		final WB_Point point = gf.createPoint();
		final WB_Map2D context = gf.createEmbeddedPlane(getPlane());
		HE_Halfedge he = _halfedge;
		int i = 0;
		do {
			context.mapPoint3D(he.getVertex(), point);
			coords[i] = new Coordinate(point.xd(), point.yd(), i);
			he = he.getNextInFace();
			i++;
		} while (he != _halfedge);
		context.mapPoint3D(he.getVertex(), point);
		coords[i] = new Coordinate(point.xd(), point.yd(), i);
		he = he.getNextInFace();
		final Polygon inputPolygon = new GeometryFactory().createPolygon(coords);
		final IsValidOp isValidOp = new IsValidOp(inputPolygon);
		if (!IsValidOp.isValid(inputPolygon)) {
			System.out.println(this);
			System.out.println(this.getFaceArea() + " " + this.getFaceNormal());
			he = _halfedge;
			i = 0;
			do {
				System.out.println("  " + i + ": " + he.getVertex());
				he = he.getNextInFace();
				i++;
			} while (he != _halfedge);
			System.out.println(isValidOp.getValidationError());
		}
	}

	/**
	 *
	 *
	 * @return
	 */
	public int getTextureId() {
		return textureId;
	}

	/**
	 *
	 *
	 * @param i
	 */
	public void setTextureId(final int i) {
		textureId = i;
	}

	public void update() {
		triangles = null;
		normal = null;
		center = null;

	}
}
