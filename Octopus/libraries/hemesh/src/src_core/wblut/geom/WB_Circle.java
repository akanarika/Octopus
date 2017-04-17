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

import wblut.math.WB_Epsilon;
import wblut.math.WB_Math;

public class WB_Circle implements WB_Curve {

	private WB_Point center;

	private WB_Vector normal;

	private double radius, r2;

	private WB_GeometryFactory geometryfactory = new WB_GeometryFactory();

	/**
	 *
	 */
	public WB_Circle() {
		center = geometryfactory.createPoint();
		normal = geometryfactory.createVector(0, 0, 1);
		radius = 1;
		r2 = 1;
	}

	/**
	 *
	 *
	 * @param center
	 * @param radius
	 */
	public WB_Circle(final WB_Coord center, final double radius) {
		this.center = geometryfactory.createPoint(center);
		this.radius = WB_Math.fastAbs(radius);
		r2 = this.radius * this.radius;
		normal = geometryfactory.createVector(0, 0, 1);
	}

	/**
	 *
	 *
	 * @param center
	 * @param normal
	 * @param radius
	 */
	public WB_Circle(final WB_Coord center, final WB_Coord normal, final double radius) {
		this.center = geometryfactory.createPoint(center);
		this.radius = WB_Math.fastAbs(radius);
		r2 = this.radius * this.radius;
		this.normal = geometryfactory.createNormalizedVector(normal);
	}

	/**
	 *
	 *
	 * @param x
	 * @param y
	 * @param r
	 */
	public WB_Circle(final double x, final double y, final double r) {
		center = geometryfactory.createPoint(x, y);
		radius = WB_Math.fastAbs(r);
		normal = geometryfactory.createVector(0, 0, 1);
	}

	/**
	 *
	 *
	 * @return
	 */
	public double getRadius() {
		return radius;
	}

	/**
	 *
	 *
	 * @return
	 */
	public WB_Coord getCenter() {
		return center;
	}

	/**
	 *
	 *
	 * @return
	 */
	public WB_Coord getNormal() {
		return normal;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.geom.WB_Geometry#apply(wblut.geom.WB_Transform)
	 */

	public WB_Circle apply(final WB_Transform T) {
		WB_Point p = geometryfactory.createPoint(center).applyAsPointSelf(T);
		WB_Point q = geometryfactory.createPoint(center).addSelf(radius, 0, 0).applyAsPointSelf(T);
		double newradius = p.getDistance2D(q);
		return geometryfactory.createCircleWithRadius(p, geometryfactory.createVector(normal).applyAsNormalSelf(T),
				newradius);
	}

	/**
	 *
	 *
	 * @param T
	 * @return
	 */
	public WB_Circle applySelf(final WB_Transform T) {
		WB_Point p = geometryfactory.createPoint(center).applyAsPointSelf(T);
		WB_Point q = geometryfactory.createPoint(center).addSelf(radius, 0, 0).applyAsPointSelf(T);
		double newradius = p.getDistance2D(q);
		center.set(p);
		normal.applyAsNormalSelf(T);
		radius = newradius;
		return this;
	}

	/**
	 *
	 *
	 * @param c
	 */
	public void set(final WB_Circle c) {
		center.set(c.getCenter());
		normal.set(c.getNormal());
		radius = c.getRadius();

	}

	/**
	 *
	 *
	 * @param x
	 * @param y
	 */
	public void setCenter(final double x, final double y) {
		center.set(x, y);
	}

	/**
	 *
	 *
	 * @param x
	 * @param y
	 * @param z
	 */
	public void setCenter(final double x, final double y, final double z) {
		center.set(x, y, z);
	}

	/**
	 *
	 *
	 * @param c
	 */
	public void setCenter(final WB_Coord c) {
		center.set(c);
	}

	/**
	 *
	 *
	 * @param x
	 * @param y
	 * @param z
	 */
	public void setNormal(final double x, final double y, final double z) {
		normal.set(x, y, z);
	}

	/**
	 *
	 *
	 * @param c
	 */
	public void setNormal(final WB_Coord c) {
		normal.set(c);
	}

	/**
	 *
	 *
	 * @param radius
	 */
	public void setRadius(final double radius) {
		this.radius = radius;
		r2 = this.radius * this.radius;
	}

	/**
	 *
	 *
	 * @param diameter
	 */
	public void setDiameter(final double diameter) {
		this.radius = diameter * 0.5;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(final Object o) {
		if (o == this) {
			return true;
		}
		if (!(o instanceof WB_Circle)) {
			return false;
		}
		return WB_Epsilon.isEqualAbs(radius, ((WB_Circle) o).getRadius()) && center.equals(((WB_Circle) o).getCenter())
				&& normal.equals(((WB_Circle) o).getNormal());
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode() {
		return 31 * (31 * center.hashCode() + hashCode(radius)) + normal.hashCode();
	}

	/**
	 *
	 *
	 * @param v
	 * @return
	 */
	private int hashCode(final double v) {
		final long tmp = Double.doubleToLongBits(v);
		return (int) (tmp ^ tmp >>> 32);
	}

	public boolean contains(final WB_Coord p) {
		return center.getSqDistance2D(p) <= r2;

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see wblut.geom.WB_Curve#curvePoint(double)
	 */
	@Override
	public WB_Point curvePoint(final double u) {
		// TODO Auto-generated method stub
		return null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see wblut.geom.WB_Curve#curveDirection(double)
	 */
	@Override
	public WB_Vector curveDirection(final double u) {
		// TODO Auto-generated method stub
		return null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see wblut.geom.WB_Curve#curveDerivative(double)
	 */
	@Override
	public WB_Vector curveDerivative(final double u) {
		// TODO Auto-generated method stub
		return null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see wblut.geom.WB_Curve#getLowerU()
	 */
	@Override
	public double getLowerU() {

		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see wblut.geom.WB_Curve#getUpperU()
	 */
	@Override
	public double getUpperU() {

		return 1;
	}

}
