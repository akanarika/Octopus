/*
 * This file is part of HE_Mesh, a library for creating and manipulating meshes.
 * It is dedicated to the public domain. To the extent possible under law,
 * I , Frederik Vanhoutte, have waived all copyright and related or neighboring
 * rights.
 *
 * This work is published from Belgium. (http://creativecommons.org/publicdomain/zero/1.0/)
 *
 */

package wblut.core;

public class WB_Disclaimer {

	public static final WB_Disclaimer CURRENT_DISCLAIMER = new WB_Disclaimer();

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		final String dis = "License: https://github.com/wblut/HE_Mesh#license";
		return dis;
	}

	/**
	 *
	 *
	 * @return
	 */
	public static String disclaimer() {
		return CURRENT_DISCLAIMER.toString();
	}
}