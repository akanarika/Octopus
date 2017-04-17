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

/**
 *
 * WB_EmbeddedPlane maps coordinates from world space into the coordinate system
 * associated with a plane. The plane can be X-, Y- or Z-plane with offset,
 * reverse X-, Y- or Z-plane with offset, or an arbitrary
 * {@link wblut.geom.WB_Plane}
 *
 */
public class WB_PlanarMap extends WB_CoordinateSystem3D implements WB_Map2D {
	/**
	 *
	 */
	private double offset;
	/**
	 *
	 */
	int id;
	/**
	 *
	 */
	private final WB_Transform T2D3D;
	/**
	 *
	 */
	private int mode;
	/**
	 *
	 */
	public static final int YZ = 0;
	/**
	 *
	 */
	public static final int XZ = 1;
	/**
	 *
	 */
	public static final int XY = 2;
	/**
	 *
	 */
	public static final int YZrev = 3;
	/**
	 *
	 */
	public static final int XZrev = 4;
	/**
	 *
	 */
	public static final int XYrev = 5;
	/**
	 *
	 */
	public static final int PLANE = 6;
	/**
	 *
	 */
	private WB_GeometryFactory geometryfactory = new WB_GeometryFactory();

	/**
	 *
	 */
	public WB_PlanarMap() {
		this(XY, 0);
	}

	public WB_PlanarMap(final WB_Coord c) {
		if (Math.abs(c.xd()) > Math.abs(c.yd())) {
			mode = Math.abs(c.xd()) > Math.abs(c.zd()) ? YZ : XY;
		} else {
			mode = Math.abs(c.yd()) > Math.abs(c.zd()) ? XZ : XY;
		}

		if (mode < 0 || mode > 5) {
			throw new IndexOutOfBoundsException();
		}
		if (mode == YZ) {
			set(geometryfactory.createPoint(offset, 0, 0), geometryfactory.Y(), geometryfactory.Z(),
					geometryfactory.X(), geometryfactory.WORLD());
			this.mode = YZ;
		} else if (mode == XZ) {
			set(geometryfactory.createPoint(0, offset, 0), geometryfactory.Z(), geometryfactory.X(),
					geometryfactory.Y(), geometryfactory.WORLD());
			this.mode = XZ;
		} else if (mode == YZrev) {
			set(geometryfactory.createPoint(0, offset, 0), geometryfactory.Z(), geometryfactory.Y(),
					geometryfactory.minX(), geometryfactory.WORLD());
			this.mode = YZrev;
		} else if (mode == XZrev) {
			set(geometryfactory.createPoint(0, offset, 0), geometryfactory.X(), geometryfactory.Z(),
					geometryfactory.minY(), geometryfactory.WORLD());
			this.mode = XZrev;
		} else if (mode == XYrev) {
			set(geometryfactory.createPoint(0, offset, 0), geometryfactory.Y(), geometryfactory.X(),
					geometryfactory.minZ(), geometryfactory.WORLD());
			this.mode = XYrev;
		} else {// XY
			set(geometryfactory.createPoint(0, 0, offset), geometryfactory.X(), geometryfactory.Y(),
					geometryfactory.Z(), geometryfactory.WORLD());
			this.mode = XY;
		}
		T2D3D = getTransformToWorld();

	}

	/**
	 *
	 *
	 * @param mode
	 * @param offset
	 */
	public WB_PlanarMap(final int mode, final double offset) {
		super();
		this.mode = mode;
		this.offset = offset;
		if (mode < 0 || mode > 5) {
			throw new IndexOutOfBoundsException();
		}
		if (mode == YZ) {
			set(geometryfactory.createPoint(offset, 0, 0), geometryfactory.Y(), geometryfactory.Z(),
					geometryfactory.X(), geometryfactory.WORLD());
			this.mode = YZ;
		} else if (mode == XZ) {
			set(geometryfactory.createPoint(0, offset, 0), geometryfactory.Z(), geometryfactory.X(),
					geometryfactory.Y(), geometryfactory.WORLD());
			this.mode = XZ;
		} else if (mode == YZrev) {
			set(geometryfactory.createPoint(0, offset, 0), geometryfactory.Z(), geometryfactory.Y(),
					geometryfactory.minX(), geometryfactory.WORLD());
			this.mode = YZrev;
		} else if (mode == XZrev) {
			set(geometryfactory.createPoint(0, offset, 0), geometryfactory.X(), geometryfactory.Z(),
					geometryfactory.minY(), geometryfactory.WORLD());
			this.mode = XZrev;
		} else if (mode == XYrev) {
			set(geometryfactory.createPoint(0, offset, 0), geometryfactory.Y(), geometryfactory.X(),
					geometryfactory.minZ(), geometryfactory.WORLD());
			this.mode = XYrev;
		} else {// XY
			set(geometryfactory.createPoint(0, 0, offset), geometryfactory.X(), geometryfactory.Y(),
					geometryfactory.Z(), geometryfactory.WORLD());
			this.mode = XY;
		}
		T2D3D = getTransformToWorld();
	}

	/**
	 *
	 *
	 * @param mode
	 */
	public WB_PlanarMap(final int mode) {
		this(mode, 0);
	}

	/**
	 *
	 *
	 * @param P
	 */
	public WB_PlanarMap(final WB_Plane P) {
		super(P.getOrigin(), P.getU(), P.getV(), P.getW(), new WB_CoordinateSystem3D());
		mode = PLANE;
		T2D3D = getTransformToWorld();
	}

	/**
	 *
	 *
	 * @param P
	 * @param offset
	 */
	public WB_PlanarMap(final WB_Plane P, final double offset) {
		super(P.getOrigin().addMul(offset, P.getNormal()), P.getU(), P.getV(), P.getW(), new WB_CoordinateSystem3D());
		mode = PLANE;
		T2D3D = getTransformToWorld();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.geom.WB_Context2D#pointTo2D(wblut.geom.WB_Coordinate,
	 * wblut.geom.WB_MutableCoordinate)
	 */
	@Override
	public void mapPoint3D(final WB_Coord p, final WB_MutableCoord result) {
		switch (mode) {
		case XY:
			result.set(p.xd(), p.yd(), p.zd() - offset);
			break;
		case YZ:
			result.set(p.yd(), p.zd(), p.xd() - offset);
			break;
		case XZ:
			result.set(p.zd(), p.xd(), p.yd() - offset);
			break;
		case XYrev:
			result.set(p.yd(), p.xd(), offset - p.zd());
			break;
		case YZrev:
			result.set(p.zd(), p.yd(), offset - p.xd());
			break;
		case XZrev:
			result.set(p.xd(), p.zd(), offset - p.yd());
			break;
		default:
			T2D3D.applyInvAsPointInto(result, p);
		}
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.geom.WB_Context2D#pointTo2D(double, double, double,
	 * wblut.geom.WB_MutableCoordinate)
	 */
	@Override
	public void mapPoint3D(final double x, final double y, final double z, final WB_MutableCoord result) {
		switch (mode) {
		case XY:
			result.set(x, y, z - offset);
			break;
		case YZ:
			result.set(y, z, x - offset);
			break;
		case XZ:
			result.set(z, x, y - offset);
			break;
		case XYrev:
			result.set(y, x, offset - z);
			break;
		case YZrev:
			result.set(z, y, offset - x);
			break;
		case XZrev:
			result.set(x, z, offset - y);
			break;
		default:
			T2D3D.applyInvAsPointInto(result, x, y, z);
		}
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.geom.WB_Context2D#pointTo3D(wblut.geom.WB_Coordinate,
	 * wblut.geom.WB_MutableCoordinate)
	 */
	@Override
	public void unmapPoint3D(final WB_Coord p, final WB_MutableCoord result) {
		switch (mode) {
		case XY:
			result.set(p.xd(), p.yd(), p.zd() + offset);
			break;
		case YZ:
			result.set(p.zd() + offset, p.xd(), p.yd());
			break;
		case XZ:
			result.set(p.yd(), p.zd() + offset, p.xd());
			break;
		case XYrev:
			result.set(p.yd(), p.xd(), offset - p.zd());
			break;
		case YZrev:
			result.set(offset - p.zd(), p.yd(), p.xd());
			break;
		case XZrev:
			result.set(p.xd(), offset - p.zd(), p.yd());
			break;
		default:
			T2D3D.applyAsPointInto(p, result);
		}
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.geom.WB_Context2D#pointTo3D(double, double, double,
	 * wblut.geom.WB_MutableCoordinate)
	 */
	@Override
	public void unmapPoint3D(final double u, final double v, final double w, final WB_MutableCoord result) {
		switch (mode) {
		case XY:
			result.set(u, v, w + offset);
			break;
		case YZ:
			result.set(w + offset, u, v);
			break;
		case XZ:
			result.set(v, w + offset, u);
			break;
		case XYrev:
			result.set(v, u, offset - w);
			break;
		case YZrev:
			result.set(offset - w, v, u);
			break;
		case XZrev:
			result.set(u, offset - w, v);
			break;
		default:
			T2D3D.applyAsPointInto(u, v, w, result);
		}
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.geom.WB_Context2D#pointTo3D(double, double,
	 * wblut.geom.WB_MutableCoordinate)
	 */
	@Override
	public void unmapPoint2D(final double u, final double v, final WB_MutableCoord result) {
		switch (mode) {
		case XY:
			result.set(u, v, offset);
			break;
		case YZ:
			result.set(offset, u, v);
			break;
		case XZ:
			result.set(v, offset, u);
			break;
		case XYrev:
			result.set(v, u, offset);
			break;
		case YZrev:
			result.set(offset, v, u);
			break;
		case XZrev:
			result.set(u, offset, v);
			break;
		default:
			T2D3D.applyAsPointInto(u, v, 0, result);
		}
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.geom.WB_Map2D#unmapPoint2D(wblut.geom.WB_Coord,
	 * wblut.geom.WB_MutableCoord)
	 */
	@Override
	public void unmapPoint2D(final WB_Coord p, final WB_MutableCoord result) {
		switch (mode) {
		case XY:
			result.set(p.xf(), p.yf(), offset);
			break;
		case YZ:
			result.set(offset, p.xf(), p.yf());
			break;
		case XZ:
			result.set(p.yf(), offset, p.xf());
			break;
		case XYrev:
			result.set(p.yf(), p.xf(), offset);
			break;
		case YZrev:
			result.set(offset, p.yf(), p.xf());
			break;
		case XZrev:
			result.set(p.xf(), offset, p.yf());
			break;
		default:
			T2D3D.applyAsPointInto(p.xf(), p.yf(), 0, result);
		}
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.geom.WB_Context2D#vectorTo2D(wblut.geom.WB_Coordinate,
	 * wblut.geom.WB_MutableCoordinate)
	 */
	@Override
	public void mapVector3D(final WB_Coord v, final WB_MutableCoord result) {
		switch (mode) {
		case XY:
			result.set(v.xd(), v.yd(), v.zd() - offset);
			break;
		case YZ:
			result.set(v.yd(), v.zd(), v.xd() - offset);
			break;
		case XZ:
			result.set(v.zd(), v.xd(), v.yd() - offset);
			break;
		case XYrev:
			result.set(v.yd(), v.xd(), offset - v.zd());
			break;
		case YZrev:
			result.set(v.zd(), v.yd(), offset - v.xd());
			break;
		case XZrev:
			result.set(v.xd(), v.zd(), offset - v.yd());
			break;
		default:
			T2D3D.applyInvAsVectorInto(result, v);
		}
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.geom.WB_Context2D#vectorTo2D(double, double, double,
	 * wblut.geom.WB_MutableCoordinate)
	 */
	@Override
	public void mapVector3D(final double x, final double y, final double z, final WB_MutableCoord result) {
		switch (mode) {
		case XY:
			result.set(x, y, z - offset);
			break;
		case YZ:
			result.set(y, z, x - offset);
			break;
		case XZ:
			result.set(z, x, y - offset);
			break;
		case XYrev:
			result.set(y, x, offset - z);
			break;
		case YZrev:
			result.set(z, y, offset - x);
			break;
		case XZrev:
			result.set(x, z, offset - y);
			break;
		default:
			T2D3D.applyInvAsVectorInto(result, x, y, z);
		}
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.geom.WB_Context2D#vectorTo3D(wblut.geom.WB_Coordinate,
	 * wblut.geom.WB_MutableCoordinate)
	 */
	@Override
	public void unmapVector3D(final WB_Coord v, final WB_MutableCoord result) {
		switch (mode) {
		case XY:
			result.set(v.xd(), v.yd(), v.zd() + offset);
			break;
		case YZ:
			result.set(v.zd() + offset, v.xd(), v.yd());
			break;
		case XZ:
			result.set(v.yd(), v.zd() + offset, v.xd());
			break;
		case XYrev:
			result.set(v.yd(), v.xd(), offset - v.zd());
			break;
		case YZrev:
			result.set(offset - v.zd(), v.yd(), v.xd());
			break;
		case XZrev:
			result.set(v.xd(), offset - v.zd(), v.yd());
			break;
		default:
			T2D3D.applyAsVectorInto(result, v);
		}
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.geom.WB_Context2D#vectorTo3D(double, double, double,
	 * wblut.geom.WB_MutableCoordinate)
	 */
	@Override
	public void unmapVector3D(final double u, final double v, final double w, final WB_MutableCoord result) {
		switch (mode) {
		case XY:
			result.set(u, v, w + offset);
			break;
		case YZ:
			result.set(w + offset, u, v);
			break;
		case XZ:
			result.set(v, w + offset, u);
			break;
		case XYrev:
			result.set(v, u, offset - w);
			break;
		case YZrev:
			result.set(offset - w, v, u);
			break;
		case XZrev:
			result.set(u, offset - w, v);
			break;
		default:
			T2D3D.applyAsVectorInto(result, u, v, w);
		}
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.geom.WB_Context2D#vectorTo3D(double, double,
	 * wblut.geom.WB_MutableCoordinate)
	 */
	@Override
	public void unmapVector2D(final double u, final double v, final WB_MutableCoord result) {
		switch (mode) {
		case XY:
			result.set(u, v, offset);
			break;
		case YZ:
			result.set(0, u, v + offset);
			break;
		case XZ:
			result.set(v, 0, u + offset);
			break;
		case XYrev:
			result.set(v, u, offset);
			break;
		case YZrev:
			result.set(offset, v, u);
			break;
		case XZrev:
			result.set(u, offset, v);
			break;
		default:
			T2D3D.applyAsVectorInto(result, u, v, 0);
		}
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.geom.WB_Map2D#unmapVector2D(wblut.geom.WB_Coord,
	 * wblut.geom.WB_MutableCoord)
	 */
	@Override
	public void unmapVector2D(final WB_Coord v, final WB_MutableCoord result) {
		switch (mode) {
		case XY:
			result.set(v.xf(), v.yf(), offset);
			break;
		case YZ:
			result.set(offset, v.xf(), v.yf());
			break;
		case XZ:
			result.set(v.yf(), offset, v.xf());
			break;
		case XYrev:
			result.set(v.yf(), v.xf(), offset);
			break;
		case YZrev:
			result.set(offset, v.yf(), v.xf());
			break;
		case XZrev:
			result.set(v.xf(), offset, v.yf());
			break;
		default:
			T2D3D.applyAsVectorInto(result, v.xf(), v.yf(), 0);
		}
	}
}
