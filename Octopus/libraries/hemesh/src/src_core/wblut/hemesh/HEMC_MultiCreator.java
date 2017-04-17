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
 * 
 */
public abstract class HEMC_MultiCreator extends HE_Machine {

	/**
	 * 
	 */
	protected int _numberOfMeshes;

	/**
	 * 
	 */
	public HEMC_MultiCreator() {
		super();
		_numberOfMeshes = 0;
	}

	/**
	 * 
	 *
	 * @return
	 */
	public HE_MeshCollection create() {
		final HE_MeshCollection result = new HE_MeshCollection();
		return result;
	}

	/**
	 * 
	 *
	 * @return
	 */
	public int numberOfMeshes() {
		return _numberOfMeshes;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.hemesh.HE_Machine#apply(wblut.hemesh.HE_Mesh)
	 */
	@Override
	public HE_Mesh apply(final HE_Mesh mesh) {
		mesh.setNoCopy(create().getMesh(0));
		return mesh;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.hemesh.HE_Machine#apply(wblut.hemesh.HE_Selection)
	 */
	@Override
	public HE_Mesh apply(final HE_Selection sel) {
		return create().getMesh(0);
	}
}
