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

import java.util.Collections;
import java.util.Comparator;

/**
 *
 *
 * @author Frederik Vanhoutte (W:Blut)
 *
 */
public class HEM_KeepLargestParts extends HEM_Modifier {

	/**
	 *
	 */
	private int n;

	/**
	 *
	 */
	public HEM_KeepLargestParts() {
		super();
		n = 1;
	}

	public HEM_KeepLargestParts(final int number) {
		super();
		n = number;
	}

	/**
	 *
	 *
	 * @param n
	 * @return
	 */
	public HEM_KeepLargestParts setNumberOfParts(final int n) {
		this.n = n;
		return this;
	}


	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.hemesh.HE_Modifier#apply(wblut.hemesh.HE_Mesh)
	 */
	@Override
	public HE_Mesh apply(final HE_Mesh mesh) {
		final HEMC_Explode explode = new HEMC_Explode().setMesh(mesh);
		final HE_MeshCollection fragments = explode.create();
		Collections.sort(fragments.meshes, new MeshSizeComparator());
		mesh.clear();


		for (int i=0;i<Math.min(n,fragments.size());i++) {

			mesh.add(fragments.getMesh(i));
		}



		return mesh;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.hemesh.HE_Modifier#apply(wblut.hemesh.HE_Mesh)
	 */
	@Override
	public HE_Mesh apply(final HE_Selection selection) {
		return apply(selection.parent);
	}


	static class MeshSizeComparator implements Comparator<HE_Mesh>
	{

		/* (non-Javadoc)
		 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
		 */
		@Override
		public int compare(final HE_Mesh mesh1, final HE_Mesh mesh2)
		{
			return 0-Integer.compare(mesh1.getNumberOfFaces(), mesh2.getNumberOfFaces());
		}
	}
}
