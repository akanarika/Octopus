/*
 * This file is part of HE_Mesh, a library for creating and manipulating meshes.
 * It is dedicated to the public domain. To the extent possible under law,
 * I , Frederik Vanhoutte, have waived all copyright and related or neighboring
 * rights.
 * 
 * This work is published from Belgium. (http://creativecommons.org/publicdomain/zero/1.0/)
 * 
 */
package wblut.math;

public interface WB_Noise {
	
	/**
	 * 
	 *
	 * @param seed 
	 */
	public void setSeed(long seed);

	/**
	 * 
	 *
	 * @param x 
	 * @return 
	 */
	public double value1D(double x);

	/**
	 * 
	 *
	 * @param x 
	 * @param y 
	 * @return 
	 */
	public double value2D(double x, double y);

	/**
	 * 
	 *
	 * @param x 
	 * @param y 
	 * @param z 
	 * @return 
	 */
	public double value3D(double x, double y, double z);

	/**
	 * 
	 *
	 * @param x 
	 * @param y 
	 * @param z 
	 * @param w 
	 * @return 
	 */
	public double value4D(double x, double y, double z, double w);

	/**
	 * 
	 *
	 * @param sx 
	 */
	public void setScale(double sx);

	/**
	 * 
	 *
	 * @param sx 
	 * @param sy 
	 */
	public void setScale(double sx, double sy);

	/**
	 * 
	 *
	 * @param sx 
	 * @param sy 
	 * @param sz 
	 */
	public void setScale(double sx, double sy, double sz);

	/**
	 * 
	 *
	 * @param sx 
	 * @param sy 
	 * @param sz 
	 * @param sw 
	 */
	public void setScale(double sx, double sy, double sz, double sw);
}
