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

import wblut.core.WB_ProgressCounter;
import wblut.geom.WB_Vector;
import wblut.math.WB_Epsilon;

/**
 *
 */
public class HEM_Triangulate extends HEM_Modifier {
	/**
	 *
	 */
	public HE_Selection triangles;

	/**
	 *
	 */
	public HEM_Triangulate() {
		super();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.hemesh.HE_Modifier#apply(wblut.hemesh.HE_Mesh)
	 */
	@Override
	public HE_Mesh apply(final HE_Mesh mesh) {
		triangles = new HE_Selection(mesh);
		tracker.setStatus(this, "Starting HEM_Triangulate.", +1);
		final HE_Face[] f = mesh.getFacesAsArray();
		final int n = mesh.getNumberOfFaces();
		WB_ProgressCounter counter = new WB_ProgressCounter(n, 10);
		tracker.setStatus(this, "Triangulating faces.", counter);
		for (int i = 0; i < n; i++) {
			if (!WB_Epsilon.isZero(WB_Vector.getLength3D(f[i].getFaceNormal()))) {
				triangulateNoPairing(f[i], mesh);
			} else {
				final HE_Halfedge he = f[i].getHalfedge();
				do {

					mesh.clearPair(he);
					mesh.setHalfedge(he.getVertex(),he);
				} while (he != f[i].getHalfedge());
			}
			counter.increment();
		}
		mesh.pairHalfedges();
		mesh.capHalfedges();
		tracker.setStatus(this, "Exiting HEM_Triangulate.", -1);
		return mesh;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see wblut.hemesh.HE_Modifier#apply(wblut.hemesh.HE_Mesh)
	 */
	@Override
	public HE_Mesh apply(final HE_Selection selection) {
		triangles = new HE_Selection(selection.parent);
		tracker.setStatus(this, "Starting HEM_Triangulate.", +1);
		final HE_Face[] f = selection.getFacesAsArray();
		final int n = selection.getNumberOfFaces();
		WB_ProgressCounter counter = new WB_ProgressCounter(n, 10);
		tracker.setStatus(this, "Triangulating faces.", counter);
		for (int i = 0; i < n; i++) {
			if (!WB_Epsilon.isZero(WB_Vector.getLength3D(f[i].getFaceNormal()))) {
				triangulateNoPairing(f[i], selection.parent);
			} else {
				final HE_Halfedge he = f[i].getHalfedge();
				do {
					selection.parent.clearPair(he);
					selection.parent.setHalfedge(he.getVertex(),he);
				} while (he != f[i].getHalfedge());
			}
			counter.increment();
		}
		selection.parent.pairHalfedges();
		selection.parent.capHalfedges();
		selection.clearFaces();
		selection.add(triangles);
		tracker.setStatus(this, "Exiting HEM_Triangulate.", -1);
		return selection.parent;
	}

	/**
	 *
	 *
	 * @param face
	 * @param mesh
	 */
	private void triangulateNoPairing(final HE_Face face, final HE_Mesh mesh) {
		if (face.getFaceOrder() == 3) {
			triangles.add(face);
		} else if (face.getFaceOrder() > 3) {
			final int[] tris = face.getTriangles(false);
			final List<HE_Vertex> vertices = face.getFaceVertices();
			final List<HE_TextureCoordinate> UVWs = face.getFaceUVWs();
			HE_Halfedge he = face.getHalfedge();
			do {

				mesh.clearPair(he);
				mesh.remove(he);
				he = he.getNextInFace();
			} while (he != face.getHalfedge());
			for (int i = 0; i < tris.length; i += 3) {

				final HE_Face f = new HE_Face();
				mesh.add(f);
				triangles.add(f);
				f.copyProperties(face);
				final HE_Halfedge he1 = new HE_Halfedge();
				final HE_Halfedge he2 = new HE_Halfedge();
				final HE_Halfedge he3 = new HE_Halfedge();
				he1.setUVW(UVWs.get(tris[i]));
				he2.setUVW(UVWs.get(tris[i + 1]));
				he3.setUVW(UVWs.get(tris[i + 2]));
				mesh.setVertex(he1,vertices.get(tris[i]));
				mesh.setVertex(he2,vertices.get(tris[i + 1]));
				mesh.setVertex(he3,vertices.get(tris[i + 2]));
				mesh.setHalfedge(he1.getVertex(),he1);
				mesh.setHalfedge(he2.getVertex(),he2);
				mesh.setHalfedge(he3.getVertex(),he3);
				mesh.setFace(he1,f);
				mesh.setFace(he2,f);
				mesh.setFace(he3,f);
				mesh.setNext(he1,he2);
				mesh.setNext(he2,he3);
				mesh.setNext(he3,he1);
				mesh.setHalfedge(f,he1);
				mesh.add(he1);
				mesh.add(he2);
				mesh.add(he3);
			}
			mesh.remove(face);
		}
	}
}
