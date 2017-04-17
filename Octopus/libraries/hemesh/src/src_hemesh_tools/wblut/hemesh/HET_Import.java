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


/**
 * The Class HET_Import.
 */
public class HET_Import {


	/**
	 * 
	 *
	 * @param path 
	 * @return 
	 */
	public static HE_Mesh readFromHemeshFile(final String path) {
		HEC_FromHemeshFile ff=new HEC_FromHemeshFile();
		ff.setPath(path);
		return new HE_Mesh(ff);
	}

	/**
	 * 
	 *
	 * @param path 
	 * @return 
	 */
	public static HE_Mesh readFromBinaryHemeshFile(final String path) {
		HEC_FromBinaryHemeshFile ff=new HEC_FromBinaryHemeshFile();
		ff.setPath(path);
		return new HE_Mesh(ff);
	}

	/**
	 * 
	 *
	 * @param path 
	 * @return 
	 */
	public static HE_Mesh readFromBinarySTLFile(final String path) {
		HEC_FromBinarySTLFile ff=new HEC_FromBinarySTLFile();
		ff.setPath(path);
		return new HE_Mesh(ff);
	}

	/**
	 * 
	 *
	 * @param path 
	 * @return 
	 */
	public static HE_Mesh readFromOBJFile(final String path) {
		HEC_FromOBJFile ff=new HEC_FromOBJFile();
		ff.setPath(path);
		return new HE_Mesh(ff);
	}

	/**
	 * 
	 *
	 * @param path 
	 * @return 
	 */
	public static HE_Mesh readFromSimpleMeshFile(final String path) {
		HEC_FromSimpleMeshFile ff=new HEC_FromSimpleMeshFile();
		ff.setPath(path);
		return new HE_Mesh(ff);
	}

}
