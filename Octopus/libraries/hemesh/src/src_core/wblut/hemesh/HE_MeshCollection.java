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

import java.util.List;

import javolution.util.FastTable;

/** OBject consisting of several meshes */
public class HE_MeshCollection {
	final List<HE_Mesh> meshes;

	/**
	 * 
	 */
	public HE_MeshCollection() {
		meshes = new FastTable<HE_Mesh>();
	}

	/**
	 * 
	 *
	 * @return 
	 */
	public HE_MeshIterator mItr() {
		return new HE_MeshIterator(this);
	}

	/**
	 * 
	 *
	 * @param i 
	 * @return 
	 */
	public HE_Mesh getMesh(final int i) {
		return meshes.get(i);
	}

	/**
	 * 
	 *
	 * @param i 
	 * @param mesh 
	 * @return 
	 */
	public HE_Mesh set(final int i, final HE_Mesh mesh) {
		return meshes.set(i, mesh);
	}

	/**
	 * 
	 *
	 * @param mesh 
	 * @return 
	 */
	public boolean add(final HE_Mesh mesh) {
		return meshes.add(mesh);
	}

	/**
	 * 
	 *
	 * @param mesh 
	 * @return 
	 */
	public boolean remove(final HE_Mesh mesh) {
		return meshes.remove(mesh);
	}

	/**
	 * 
	 *
	 * @return 
	 */
	public int getNumberOfMeshes() {
		return meshes.size();
	}

	/**
	 * 
	 *
	 * @return 
	 */
	public int size() {
		return meshes.size();
	}
}
