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

import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import gnu.trove.iterator.TLongObjectIterator;
import gnu.trove.map.TLongIntMap;
import gnu.trove.map.TLongObjectMap;
import gnu.trove.map.hash.TLongIntHashMap;
import gnu.trove.map.hash.TLongObjectHashMap;
import javolution.util.FastMap;
import javolution.util.FastTable;
import wblut.core.WB_ProgressCounter;
import wblut.geom.WB_AABB;
import wblut.geom.WB_Coord;
import wblut.geom.WB_GeometryFactory;
import wblut.geom.WB_GeometryOp3D;
import wblut.geom.WB_Point;
import wblut.geom.WB_Polygon;
import wblut.geom.WB_Segment;
import wblut.geom.WB_Sphere;
import wblut.hemesh.HE_RAS.HE_RASTrove;

/**
*
*/

/**
 * Collection of mesh elements. Contains methods to manipulate the data
 * structures.
 *
 * @author Frederik Vanhoutte (W:Blut)
 *
 */
public class HE_MeshStructure extends HE_MeshElement {
	protected WB_GeometryFactory gf = new WB_GeometryFactory();
	protected HE_RAS<HE_Vertex> vertices;
	protected HE_RAS<HE_Halfedge> halfedges;
	protected HE_RAS<HE_Halfedge> edges;
	protected HE_RAS<HE_Halfedge> unpairedHalfedges;
	protected HE_RAS<HE_Face> faces;

	/**
	 * Instantiates a new HE_MeshStructure.
	 */
	public HE_MeshStructure() {
		super();
		vertices = new HE_RAS.HE_RASTrove<HE_Vertex>();
		halfedges = new HE_RAS.HE_RASTrove<HE_Halfedge>();
		edges = new HE_RAS.HE_RASTrove<HE_Halfedge>();
		unpairedHalfedges = new HE_RAS.HE_RASTrove<HE_Halfedge>();
		faces = new HE_RAS.HE_RASTrove<HE_Face>();

	}

	/**
	 * Instantiates a new HE_MeshStructure and populate it with a shallow copy
	 * of mesh structure as parameter.
	 *
	 * @param ms
	 *            mesh structure to copy
	 */
	public HE_MeshStructure(final HE_MeshStructure ms) {
		this();
		for (final HE_Face f : ms.faces) {
			add(f);
		}
		for (final HE_Halfedge he : ms.halfedges) {
			halfedges.add(he);
		}

		for (final HE_Halfedge e : ms.edges) {
			edges.add(e);
		}

		for (final HE_Halfedge uhe : ms.unpairedHalfedges) {
			unpairedHalfedges.add(uhe);
		}
		for (final HE_Vertex v : ms.vertices) {
			add(v);
		}
	}

	/**
	 * Number of faces.
	 *
	 * @return the number of faces
	 */
	public final int getNumberOfFaces() {
		return faces.size();
	}

	/**
	 * Number of halfedges.
	 *
	 * @return the number of halfedges
	 */
	public final int getNumberOfHalfedges() {
		return halfedges.size() + edges.size() + unpairedHalfedges.size();
	}

	/**
	 * Number of edges.
	 *
	 * @return the number of edges
	 */
	public int getNumberOfEdges() {
		return edges.size();
	}

	/**
	 * Number of vertices.
	 *
	 * @return the number of vertices
	 */
	public final int getNumberOfVertices() {
		return vertices.size();
	}

	/**
	 * Get face with key. The key of a mesh element is unique and never changes.
	 *
	 * @param key
	 *            face key
	 * @return face
	 */
	public final HE_Face getFaceWithKey(final long key) {
		return faces.getWithKey(key);
	}

	/**
	 * Get halfedge with key. The key of a mesh element is unique and never
	 * changes.
	 *
	 * @param key
	 *            halfedge key
	 * @return halfedge
	 */
	public final HE_Halfedge getHalfedgeWithKey(final long key) {
		HE_Halfedge he = edges.getWithKey(key);
		if (he != null) {
			return he;
		}
		he = halfedges.getWithKey(key);
		if (he != null) {
			return he;
		}
		return unpairedHalfedges.getWithKey(key);
	}

	/**
	 * Get edge with key. The key of a mesh element is unique and never changes.
	 *
	 * @param key
	 *            halfedge key
	 * @return halfedge
	 */
	public final HE_Halfedge getEdgeWithKey(final long key) {
		HE_Halfedge he = edges.getWithKey(key);
		if (he != null) {
			return he;
		}
		he = halfedges.getWithKey(key);
		if (he != null) {
			return he;
		}
		return unpairedHalfedges.getWithKey(key);
	}

	/**
	 * Get vertex with key. The key of a mesh element is unique and never
	 * changes.
	 *
	 * @param key
	 *            vertex key
	 * @return vertex
	 */
	public final HE_Vertex getVertexWithKey(final long key) {
		return vertices.getWithKey(key);
	}

	/**
	 * Get face with index. Indices of mesh elements are not fixed and will
	 * change when the mesh is modified.
	 *
	 * @param i
	 *            face index
	 * @return
	 */
	public final HE_Face getFaceWithIndex(final int i) {
		return faces.get(i);
	}

	/**
	 * Get halfedge with index. Indices of mesh elements are not fixed and will
	 * change when the mesh is modified.
	 *
	 * @param i
	 *            halfedge index
	 * @return
	 */
	public final HE_Halfedge getHalfedgeWithIndex(final int i) {
		if (i >= edges.size() + halfedges.size()) {
			return unpairedHalfedges.get(i - edges.size() - halfedges.size());
		} else if (i >= edges.size()) {
			return halfedges.get(i - edges.size());
		}
		return edges.get(i);
	}

	/**
	 * Get edge with index. Indices of mesh elements are not fixed and will
	 * change when the mesh is modified.
	 *
	 * @param i
	 *            edge index
	 * @return
	 */
	public final HE_Halfedge getEdgeWithIndex(final int i) {
		return edges.get(i);
	}

	/**
	 * Get vertex with index. Indices of mesh elements are not fixed and will
	 * change when the mesh is modified.
	 *
	 * @param i
	 *            vertex index
	 * @return
	 */
	public final HE_Vertex getVertexWithIndex(final int i) {
		return vertices.get(i);
	}

	/**
	 * Add face.
	 *
	 * @param f
	 *            face to add
	 */
	public final void add(final HE_Face f) {
		faces.add(f);
	}

	/**
	 * Adds halfedge.
	 *
	 * @param he
	 *            halfedge to add
	 */
	public final void add(final HE_Halfedge he) {
		if (he.getPair() == null) {
			unpairedHalfedges.add(he);
		} else if (he.isEdge()) {
			edges.add(he);
		} else {
			halfedges.add(he);
		}
	}

	/**
	 * Add vertex.
	 *
	 * @param v
	 *            vertex to add
	 */
	public final void add(final HE_Vertex v) {
		vertices.add(v);
	}

	/**
	 * Add all mesh elements to this mesh. No copies are made.
	 *
	 * @param mesh
	 *            mesh to add
	 */
	public void add(final HE_Mesh mesh) {
		addVertices(mesh);
		addFaces(mesh);
		addHalfedges(mesh);
	}

	/**
	 * Adds faces.
	 *
	 * @param faces
	 *            faces to add as HE_Face[]
	 */
	public final void addFaces(final HE_Face[] faces) {
		for (final HE_Face face : faces) {
			add(face);
		}
	}

	/**
	 * Adds faces.
	 *
	 * @param faces
	 *            faces to add as Collection<? extends HE_Face>
	 */
	public final void addFaces(final Collection<? extends HE_Face> faces) {
		for (HE_Face f : faces) {
			add(f);
		}
	}

	/**
	 *
	 *
	 * @param source
	 */
	public final void addFaces(final HE_MeshStructure source) {
		faces.addAll(source.faces);
	}

	/**
	 * Adds halfedges.
	 *
	 * @param halfedges
	 *            halfedges to add as HE_Halfedge[]
	 */
	public final void addHalfedges(final HE_Halfedge[] halfedges) {
		for (final HE_Halfedge halfedge : halfedges) {
			add(halfedge);
		}
	}

	/**
	 * Adds halfedges.
	 *
	 * @param halfedges
	 *            halfedges to add as Collection<? extends HE_Halfedge>
	 */
	public final void addHalfedges(final Collection<? extends HE_Halfedge> halfedges) {
		for (HE_Halfedge he : halfedges) {
			add(he);
		}
	}

	/**
	 *
	 *
	 * @param source
	 */
	public final void addHalfedges(final HE_MeshStructure source) {
		halfedges.addAll(source.halfedges);
		edges.addAll(source.edges);
		unpairedHalfedges.addAll(source.unpairedHalfedges);

	}

	/**
	 *
	 *
	 * @param source
	 */
	public final void addEdges(final HE_MeshStructure source) {
		edges.addAll(source.edges);
	}

	/**
	 * Adds vertices.
	 *
	 * @param vertices
	 *            vertices to add as HE_Vertex[]
	 */
	public final void addVertices(final HE_Vertex[] vertices) {
		for (final HE_Vertex vertex : vertices) {
			add(vertex);
		}
	}

	/**
	 *
	 *
	 * @param source
	 */
	public final void addVertices(final HE_MeshStructure source) {
		vertices.addAll(source.vertices);
	}

	/**
	 * Adds vertices.
	 *
	 * @param vertices
	 *            vertices to add as Collection<? extends HE_Vertex>
	 */
	public final void addVertices(final Collection<? extends HE_Vertex> vertices) {
		for (HE_Vertex v : vertices) {
			add(v);
		}
	}

	/**
	 * Removes face.
	 *
	 * @param f
	 *            face to remove
	 */
	public final void remove(final HE_Face f) {

		faces.remove(f);
	}

	/**
	 * Removes halfedge.
	 *
	 * @param he
	 *            halfedge to remove
	 */
	public final void remove(final HE_Halfedge he) {
		edges.remove(he);
		halfedges.remove(he);
		unpairedHalfedges.remove(he);
	}

	/**
	 * Removes vertex.
	 *
	 * @param v
	 *            vertex to remove
	 */
	public final void remove(final HE_Vertex v) {
		vertices.remove(v);
	}

	/**
	 * Removes faces.
	 *
	 * @param faces
	 *            faces to remove as HE_Face[]
	 */
	public final void removeFaces(final HE_Face[] faces) {
		for (final HE_Face face : faces) {
			remove(face);
		}
	}

	/**
	 * Removes faces.
	 *
	 * @param faces
	 *            faces to remove as Collection<? extends HE_Face>
	 */
	public final void removeFaces(final Collection<? extends HE_Face> faces) {
		for (final HE_Face f : faces) {
			remove(f);
		}
	}

	/**
	 * Removes halfedges.
	 *
	 * @param halfedges
	 *            halfedges to remove as HE_Halfedge[]
	 */
	public final void removeHalfedges(final HE_Halfedge[] halfedges) {
		for (final HE_Halfedge halfedge : halfedges) {
			remove(halfedge);
		}
	}

	/**
	 * Removes halfedges.
	 *
	 * @param halfedges
	 *            halfedges to remove as Collection<? extends HE_Halfedge>
	 */
	public final void removeHalfedges(final Collection<? extends HE_Halfedge> halfedges) {
		for (final HE_Halfedge he : halfedges) {
			remove(he);
		}
	}

	/**
	 * Removes vertices.
	 *
	 * @param vertices
	 *            vertices to remove as HE_Vertex[]
	 */
	public final void removeVertices(final HE_Vertex[] vertices) {
		for (final HE_Vertex vertice : vertices) {
			remove(vertice);
		}
	}

	/**
	 * Removes vertices.
	 *
	 * @param vertices
	 *            vertices to remove as Collection<? extends HE_Vertex>
	 */
	public final void removeVertices(final Collection<? extends HE_Vertex> vertices) {
		for (final HE_Vertex v : vertices) {
			remove(v);
		}
	}

	/**
	 * Clear entire structure.
	 */
	@Override
	public void clear() {
		clearVertices();
		clearHalfedges();
		clearFaces();
	}

	/**
	 * Clear faces.
	 */
	public final void clearFaces() {
		faces = new HE_RASTrove<HE_Face>();
	}

	/**
	 * Clear halfedges.
	 */
	public final void clearHalfedges() {
		halfedges = new HE_RASTrove<HE_Halfedge>();
		edges = new HE_RASTrove<HE_Halfedge>();
		unpairedHalfedges = new HE_RASTrove<HE_Halfedge>();

	}

	/**
	 * Clear edges.
	 */
	protected final void clearEdges() {

		edges = new HE_RASTrove<HE_Halfedge>();

	}

	/**
	 * Clear vertices.
	 */
	public final void clearVertices() {
		vertices = new HE_RASTrove<HE_Vertex>();
	}

	/**
	 * Check if structure contains face.
	 *
	 * @param f
	 *            face
	 * @return true, if successful
	 */
	public final boolean contains(final HE_Face f) {
		return faces.contains(f);
	}

	/**
	 * Check if structure contains halfedge.
	 *
	 * @param he
	 *            halfedge
	 * @return true, if successful
	 */
	public final boolean contains(final HE_Halfedge he) {
		return edges.contains(he) || halfedges.contains(he) || unpairedHalfedges.contains(he);

	}

	/**
	 * Check if structure contains vertex.
	 *
	 * @param v
	 *            vertex
	 * @return true, if successful
	 */
	public final boolean contains(final HE_Vertex v) {
		return vertices.contains(v);
	}

	/**
	 * Get axis-aligned bounding box surrounding mesh.
	 *
	 * @return WB_AABB axis-aligned bounding box
	 */
	public final WB_AABB getAABB() {
		final double[] result = getLimits();
		final WB_Point min = gf.createPoint(result[0], result[1], result[2]);
		final WB_Point max = gf.createPoint(result[3], result[4], result[5]);
		return new WB_AABB(min, max);
	}

	/**
	 *
	 *
	 * @return
	 */
	public final WB_Sphere getBoundingSphere() {

		return WB_GeometryOp3D.getBoundingSphere(vertices);
	}

	/**
	 * Get range of vertex coordinates.
	 *
	 * @return array of limit values: min x, min y, min z, max x, max y, max z
	 */
	public final double[] getLimits() {
		final double[] result = new double[6];
		for (int i = 0; i < 3; i++) {
			result[i] = Double.POSITIVE_INFINITY;
		}
		for (int i = 3; i < 6; i++) {
			result[i] = Double.NEGATIVE_INFINITY;
		}
		HE_Vertex v;
		for (int i = 0; i < vertices.size(); i++) {
			v = getVertexWithIndex(i);
			result[0] = Math.min(result[0], v.xd());
			result[1] = Math.min(result[1], v.yd());
			result[2] = Math.min(result[2], v.zd());
			result[3] = Math.max(result[3], v.xd());
			result[4] = Math.max(result[4], v.yd());
			result[5] = Math.max(result[5], v.zd());
		}
		return result;
	}

	/**
	 *
	 *
	 * @return
	 */
	public final List<HE_Vertex> getVertices() {
		return vertices.getObjects();
	}

	/**
	 * Vertices as array.
	 *
	 * @return all vertices as HE_Vertex[]
	 */
	public final HE_Vertex[] getVerticesAsArray() {
		final HE_Vertex[] vertices = new HE_Vertex[getNumberOfVertices()];
		final Collection<HE_Vertex> _vertices = this.vertices;
		final Iterator<HE_Vertex> vitr = _vertices.iterator();
		int i = 0;
		while (vitr.hasNext()) {
			vertices[i] = vitr.next();
			i++;
		}
		return vertices;
	}

	/**
	 *
	 *
	 * @return
	 */
	public final List<HE_Halfedge> getHalfedges() {
		final List<HE_Halfedge> halfedges = new FastTable<HE_Halfedge>();
		halfedges.addAll(this.halfedges);
		halfedges.addAll(this.edges);
		halfedges.addAll(this.unpairedHalfedges);
		Collections.sort(halfedges);
		return halfedges;
	}

	/**
	 * Halfedges as array.
	 *
	 * @return all halfedges as HE_Halfedge[]
	 */
	public final HE_Halfedge[] getHalfedgesAsArray() {
		List<HE_Halfedge> hes = getHalfedges();
		final HE_Halfedge[] halfedges = new HE_Halfedge[hes.size()];
		int i = 0;
		for (HE_Halfedge he : hes) {
			halfedges[i] = he;
			i++;
		}

		return halfedges;
	}

	/**
	 *
	 *
	 * @return
	 */
	public final List<HE_Halfedge> getEdges() {
		return edges.getObjects();
	}

	/**
	 * Edges as array.
	 *
	 * @return all edges as HE_Halfedge[]
	 */
	public final HE_Halfedge[] getEdgesAsArray() {
		final HE_Halfedge[] edges = new HE_Halfedge[getNumberOfEdges()];
		final Iterator<HE_Halfedge> eItr = eItr();
		int i = 0;
		while (eItr.hasNext()) {
			edges[i] = eItr.next();
			i++;
		}
		return edges;
	}

	/**
	 *
	 *
	 * @return
	 */
	public final List<HE_Face> getFaces() {
		return faces.getObjects();
	}

	/**
	 * Faces as array.
	 *
	 * @return all faces as HE_Face[]
	 */
	public final HE_Face[] getFacesAsArray() {
		final HE_Face[] faces = new HE_Face[getNumberOfFaces()];
		final Iterator<HE_Face> fItr = this.faces.iterator();
		int i = 0;
		while (fItr.hasNext()) {
			faces[i] = fItr.next();
			i++;
		}
		return faces;
	}

	/**
	 *
	 *
	 * @param key
	 * @return
	 */
	public final boolean containsFace(final long key) {
		return faces.containsKey(key);
	}

	/**
	 *
	 *
	 * @param key
	 * @return
	 */
	public final boolean containsHalfedge(final long key) {
		return halfedges.containsKey(key) || edges.containsKey(key) || unpairedHalfedges.containsKey(key);
	}

	/**
	 *
	 *
	 * @param key
	 * @return
	 */
	public final boolean containsEdge(final long key) {
		return edges.containsKey(key);
	}

	/**
	 *
	 *
	 * @param key
	 * @return
	 */
	public final boolean containsVertex(final long key) {
		return vertices.containsKey(key);
	}

	/**
	 *
	 *
	 * @param f
	 * @return
	 */
	public final int getIndex(final HE_Face f) {
		return faces.indexOf(f);
	}

	/**
	 *
	 *
	 * @param v
	 * @return
	 */
	public final int getIndex(final HE_Vertex v) {
		return vertices.indexOf(v);
	}

	/**
	 *
	 *
	 * @param mesh
	 */
	protected final void replaceFaces(final HE_Mesh mesh) {
		clearFaces();
		addFaces(mesh);
	}

	/**
	 * Replace faces.
	 *
	 * @param faces
	 *            faces to replace with as HE_Face[]
	 */
	protected final void replaceFaces(final HE_Face[] faces) {
		clearFaces();
		for (final HE_Face face : faces) {
			add(face);
		}
	}

	/**
	 * Replace faces.
	 *
	 * @param faces
	 *            faces to replace with as Collection<? extends HE_Face>
	 */
	protected final void replaceFaces(final Collection<? extends HE_Face> faces) {
		clearFaces();
		addFaces(faces);
	}

	/**
	 *
	 *
	 * @param mesh
	 */
	protected final void replaceVertices(final HE_Mesh mesh) {
		clearVertices();
		addVertices(mesh);
	}

	/**
	 * Replace vertices.
	 *
	 * @param vertices
	 *            vertices to replace with as HE_Vertex[]
	 */
	protected final void replaceVertices(final HE_Vertex[] vertices) {
		clearVertices();
		for (final HE_Vertex vertice : vertices) {
			add(vertice);
		}
	}

	/**
	 * Replace vertices.
	 *
	 * @param vertices
	 *            vertices to replace with as Collection<? extends HE_Vertex>
	 */
	protected final void replaceVertices(final Collection<? extends HE_Vertex> vertices) {
		clearVertices();
		addVertices(vertices);
	}

	/**
	 *
	 *
	 * @param mesh
	 */
	protected final void replaceHalfedges(final HE_Mesh mesh) {
		clearHalfedges();
		HE_HalfedgeIterator heItr = mesh.heItr();
		while (heItr.hasNext()) {
			add(heItr.next());

		}

	}

	/**
	 * Replace halfedges.
	 *
	 * @param halfedges
	 *            halfedges to replace with as HE_Halfedge[]
	 */
	protected final void replaceHalfedges(final HE_Halfedge[] halfedges) {
		clearHalfedges();

		for (final HE_Halfedge halfedge : halfedges) {
			add(halfedge);
		}
	}

	/**
	 * Replace halfedges.
	 *
	 * @param halfedges
	 *            halfedges to replace with as Collection<? extends HE_Halfedge>
	 */
	protected final void replaceHalfedges(final Collection<? extends HE_Halfedge> halfedges) {
		clearHalfedges();

		addHalfedges(halfedges);

	}

	/**
	 * Vertex iterator.
	 *
	 * @return vertex iterator
	 */
	public HE_VertexIterator vItr() {
		return new HE_VertexIterator(vertices);
	}

	/**
	 * Edge iterator.
	 *
	 * @return edge iterator
	 */
	public HE_EdgeIterator eItr() {
		return new HE_EdgeIterator(edges);
	}

	/**
	 * Halfedge iterator.
	 *
	 * @return halfedge iterator
	 */
	public HE_HalfedgeIterator heItr() {
		return new HE_HalfedgeIterator(edges, halfedges, unpairedHalfedges);
	}

	/**
	 * Face iterator.
	 *
	 * @return face iterator
	 */
	public HE_FaceIterator fItr() {
		return new HE_FaceIterator(faces);
	}

	/**
	 * Link face to halfedge
	 *
	 * @param he
	 * @param f
	 */
	public void setFace(final HE_Halfedge he, final HE_Face f) {
		he._setFace(f);
		if (he.getPair() != null) {
			setPair(he, he.getPair());
		}

	}

	/**
	 * Unlink face from halfedge
	 *
	 * @param he
	 */
	public void clearFace(final HE_Halfedge he) {
		he._clearFace();
		if (he.getPair() != null) {
			setPair(he, he.getPair());
		}

	}

	/**
	 * Pair two halfedges
	 *
	 * @param he1
	 * @param he2
	 */
	public void setPair(final HE_Halfedge he1, final HE_Halfedge he2) {
		remove(he1);
		remove(he2);
		he1._setPair(he2);
		he2._setPair(he1);
		add(he1);
		add(he2);
	}

	/**
	 * Unpair halfedge. If the halfedge was paired, its pair is unpaired as
	 * well.
	 *
	 * @param he
	 */
	public void clearPair(final HE_Halfedge he) {
		if (he.getPair() == null) {
			return;
		}
		HE_Halfedge hep = he.getPair();
		remove(he);
		remove(hep);
		he._clearPair();
		hep._clearPair();
		add(he);
		add(hep);
	}

	/**
	 * Link hen to he as next halfedge, he is linked as previous halfedge to he.
	 *
	 * @param he
	 * @param hen
	 */
	public void setNext(final HE_Halfedge he, final HE_Halfedge hen) {
		he._setNext(hen);
		hen._setPrev(he);

	}

	/**
	 * Link vertex to halfedge
	 *
	 * @param he
	 * @param v
	 */
	public void setVertex(final HE_Halfedge he, final HE_Vertex v) {
		he._setVertex(v);

	}

	/**
	 * Link halfedge to face
	 *
	 * @param f
	 * @param he
	 */
	public void setHalfedge(final HE_Face f, final HE_Halfedge he) {
		f._setHalfedge(he);

	}

	/**
	 * Unlink halfedge from face
	 *
	 * @param f
	 */
	public void clearHalfedge(final HE_Face f) {
		f._clearHalfedge();
	}

	/**
	 * Link halfedge to vertex
	 *
	 * @param v
	 * @param he
	 */
	public void setHalfedge(final HE_Vertex v, final HE_Halfedge he) {
		v._setHalfedge(he);

	}

	/**
	 * Unlink halfedge from vertex
	 *
	 * @param v
	 */
	public void clearHalfedge(final HE_Vertex v) {
		v._clearHalfedge();
	}

	/**
	 * Unlink next halfedge from halfedge, unlinks the corresponding "previous"
	 * relationship.
	 *
	 * @param he
	 */
	public void clearNext(final HE_Halfedge he) {
		if (he.getNextInFace() != null) {
			he.getNextInFace()._clearPrev();
		}
		he._clearNext();
	}

	/**
	 * Unlink previous halfedge from halfedge, unlinks the corresponding "next"
	 * relationship.
	 *
	 * @param he
	 */
	public void clearPrev(final HE_Halfedge he) {
		if (he.getPrevInFace() != null) {
			he.getPrevInFace()._clearNext();
		}
		he._clearPrev();
	}

	/**
	 * Unlink vertex from halfedge
	 *
	 * @param he
	 */
	public void clearVertex(final HE_Halfedge he) {
		he._clearVertex();
	}

	/**
	 * Sort all faces and vertices in lexographical order
	 */
	public void sort() {
		HE_FaceIterator fItr = fItr();
		while (fItr.hasNext()) {
			fItr.next().sort();
		}

		List<HE_Face> sortedFaces = new FastTable<HE_Face>();
		sortedFaces.addAll(getFaces());
		Collections.sort(sortedFaces);
		clearFaces();
		addFaces(sortedFaces);

		List<HE_Vertex> sortedVertices = new FastTable<HE_Vertex>();
		sortedVertices.addAll(getVertices());
		Collections.sort(sortedVertices);
		clearVertices();
		addVertices(sortedVertices);
	}

	/**
	 *
	 *
	 * @param color
	 */
	public void setFaceColor(final int color) {
		final HE_FaceIterator fitr = fItr();
		while (fitr.hasNext()) {
			fitr.next().setColor(color);
		}
	}

	/**
	 *
	 *
	 * @param color
	 */
	public void setVertexColor(final int color) {
		final HE_VertexIterator vitr = vItr();
		while (vitr.hasNext()) {
			vitr.next().setColor(color);
		}
	}

	/**
	 *
	 *
	 * @param color
	 */
	public void setHalfedgeColor(final int color) {
		final HE_HalfedgeIterator heitr = heItr();
		while (heitr.hasNext()) {
			heitr.next().setColor(color);
		}
	}

	/**
	 *
	 *
	 * @param color
	 * @param i
	 */
	public void setFaceColorWithLabel(final int color, final int i) {
		final HE_FaceIterator fitr = fItr();
		HE_Face f;
		while (fitr.hasNext()) {
			f = fitr.next();
			if (f.getLabel() == i) {
				f.setColor(color);
			}
		}
	}

	/**
	 *
	 *
	 * @param color
	 * @param i
	 */
	public void setFaceColorWithOtherLabel(final int color, final int i) {
		final HE_FaceIterator fitr = fItr();
		HE_Face f;
		while (fitr.hasNext()) {
			f = fitr.next();
			if (f.getLabel() != i) {
				f.setColor(color);
			}
		}
	}

	/**
	 *
	 *
	 * @param color
	 * @param i
	 */
	public void setFaceColorWithInternalLabel(final int color, final int i) {
		final HE_FaceIterator fitr = fItr();
		HE_Face f;
		while (fitr.hasNext()) {
			f = fitr.next();
			if (f.getInternalLabel() == i) {
				f.setColor(color);
			}
		}
	}

	/**
	 *
	 *
	 * @param color
	 * @param i
	 */
	public void setFaceColorWithOtherInternalLabel(final int color, final int i) {
		final HE_FaceIterator fitr = fItr();
		HE_Face f;
		while (fitr.hasNext()) {
			f = fitr.next();
			if (f.getInternalLabel() != i) {
				f.setColor(color);
			}
		}
	}

	/**
	 *
	 *
	 * @param color
	 * @param i
	 */
	public void setVertexColorWithLabel(final int color, final int i) {
		final HE_VertexIterator fitr = vItr();
		HE_Vertex f;
		while (fitr.hasNext()) {
			f = fitr.next();
			if (f.getLabel() == i) {
				f.setColor(color);
			}
		}
	}

	/**
	 *
	 *
	 * @param color
	 * @param i
	 */
	public void setVertexColorWithOtherLabel(final int color, final int i) {
		final HE_VertexIterator fitr = vItr();
		HE_Vertex f;
		while (fitr.hasNext()) {
			f = fitr.next();
			if (f.getLabel() != i) {
				f.setColor(color);
			}
		}
	}

	/**
	 *
	 *
	 * @param color
	 * @param i
	 */
	public void setVertexColorWithInternalLabel(final int color, final int i) {
		final HE_VertexIterator fitr = vItr();
		HE_Vertex f;
		while (fitr.hasNext()) {
			f = fitr.next();
			if (f.getInternalLabel() == i) {
				f.setColor(color);
			}
		}
	}

	/**
	 *
	 *
	 * @param color
	 * @param i
	 */
	public void setVertexColorWithOtherInternalLabel(final int color, final int i) {
		final HE_VertexIterator vitr = vItr();
		HE_Vertex v;
		while (vitr.hasNext()) {
			v = vitr.next();
			if (v.getInternalLabel() != i) {
				v.setColor(color);
			}
		}
	}

	/**
	 *
	 *
	 * @param color
	 * @param i
	 */
	public void setHalfedgeColorWithLabel(final int color, final int i) {
		final HE_HalfedgeIterator fitr = heItr();
		HE_Halfedge f;
		while (fitr.hasNext()) {
			f = fitr.next();
			if (f.getLabel() == i) {
				f.setColor(color);
			}
		}
	}

	/**
	 *
	 *
	 * @param color
	 * @param i
	 */
	public void setHalfedgeColorWithOtherLabel(final int color, final int i) {
		final HE_HalfedgeIterator heitr = heItr();
		HE_Halfedge he;
		while (heitr.hasNext()) {
			he = heitr.next();
			if (he.getLabel() != i) {
				he.setColor(color);
			}
		}
	}

	/**
	 *
	 *
	 * @param color
	 * @param i
	 */
	public void setHalfedgeColorWithInternalLabel(final int color, final int i) {
		final HE_HalfedgeIterator fitr = heItr();
		HE_Halfedge f;
		while (fitr.hasNext()) {
			f = fitr.next();
			if (f.getInternalLabel() == i) {
				f.setColor(color);
			}
		}
	}

	/**
	 *
	 *
	 * @param color
	 * @param i
	 */
	public void setHalfedgeColorWithOtherInternalLabel(final int color, final int i) {
		final HE_HalfedgeIterator heitr = heItr();
		;
		HE_Halfedge f;
		while (heitr.hasNext()) {
			f = heitr.next();
			if (f.getInternalLabel() != i) {
				f.setColor(color);
			}
		}
	}

	/**
	 * Reset all labels.
	 */
	public void resetLabels() {
		resetVertexLabels();
		resetFaceLabels();
		resetHalfedgeLabels();
	}

	/**
	 * Reset all vertex labels to -1.
	 */
	public void resetVertexLabels() {
		final Iterator<HE_Vertex> vItr = vItr();
		while (vItr.hasNext()) {
			vItr.next().setLabel(-1);
		}
	}

	/**
	 * Reset all face labels to -1.
	 */
	public void resetFaceLabels() {
		final Iterator<HE_Face> fItr = fItr();
		while (fItr.hasNext()) {
			fItr.next().setLabel(-1);
		}
	}

	/**
	 * Reset all edge labels to -1.
	 */
	public void resetEdgeLabels() {
		final Iterator<HE_Halfedge> eItr = eItr();
		while (eItr.hasNext()) {
			eItr.next().setLabel(-1);
		}
	}

	/**
	 * Reset all halfedge labels to -1.
	 */
	public void resetHalfedgeLabels() {
		final Iterator<HE_Halfedge> heItr = heItr();
		while (heItr.hasNext()) {
			heItr.next().setLabel(-1);
		}
	}

	/**
	 * Set all vertex labels to value.
	 *
	 * @param label
	 */
	public void setVertexLabels(final int label) {
		final Iterator<HE_Vertex> vItr = vItr();
		while (vItr.hasNext()) {
			vItr.next().setLabel(label);
		}
	}

	/**
	 * Reset all face labels to value.
	 *
	 * @param label
	 */
	public void setFaceLabels(final int label) {
		final Iterator<HE_Face> fItr = fItr();
		while (fItr.hasNext()) {
			fItr.next().setLabel(label);
		}
	}

	/**
	 * Reset all edge labels to value.
	 *
	 * @param label
	 */
	public void setEdgeLabels(final int label) {
		final Iterator<HE_Halfedge> eItr = eItr();
		while (eItr.hasNext()) {
			eItr.next().setLabel(label);
		}
	}

	/**
	 * Reset all halfedge labels to value.
	 *
	 * @param label
	 */
	public void setHalfedgeLabels(final int label) {
		final Iterator<HE_Halfedge> heItr = heItr();
		while (heItr.hasNext()) {
			heItr.next().setLabel(label);
		}
	}

	/**
	 * Reset all internal labels to -1.
	 *
	 */

	protected void resetInternalLabels() {
		resetVertexInternalLabels();
		resetFaceInternalLabels();
		resetHalfedgeInternalLabels();
	}

	/**
	 * Reset face labels to -1.
	 *
	 *
	 */

	protected void resetFaceInternalLabels() {
		final Iterator<HE_Face> fItr = fItr();
		while (fItr.hasNext()) {
			fItr.next().setInternalLabel(-1);
		}
	}

	/**
	 * Reset all halfedge labels to -1.
	 *
	 */

	protected void resetHalfedgeInternalLabels() {
		final Iterator<HE_Halfedge> heItr = heItr();
		while (heItr.hasNext()) {
			heItr.next().setInternalLabel(-1);
		}
	}

	/**
	 * Reset all edge labels to -1.
	 *
	 */

	protected void resetEdgeInternalLabels() {
		final Iterator<HE_Halfedge> eItr = eItr();
		while (eItr.hasNext()) {
			eItr.next().setInternalLabel(-1);
		}
	}

	/**
	 * Reset all vertex labels to -1.
	 *
	 *
	 */

	protected void resetVertexInternalLabels() {
		final Iterator<HE_Vertex> vItr = vItr();
		while (vItr.hasNext()) {
			vItr.next().setInternalLabel(-1);
		}
	}

	/**
	 * Set face labels to value.
	 *
	 * @param label
	 */

	protected void setFaceInternalLabels(final int label) {
		final Iterator<HE_Face> fItr = fItr();
		while (fItr.hasNext()) {
			fItr.next().setInternalLabel(label);
		}
	}

	/**
	 * Det edge labels to value.
	 *
	 * @param label
	 *
	 */
	protected void setHalfedgeInternalLabels(final int label) {
		final Iterator<HE_Halfedge> heItr = heItr();
		while (heItr.hasNext()) {
			heItr.next().setInternalLabel(label);
		}
	}

	/**
	 * Set edge labels to value.
	 *
	 * @param label
	 *
	 */
	protected void setEdgeInternalLabels(final int label) {
		final Iterator<HE_Halfedge> eItr = eItr();
		while (eItr.hasNext()) {
			eItr.next().setInternalLabel(label);
		}
	}

	/**
	 * Set all vertex labels to value.
	 *
	 * @param label
	 */

	protected void setVertexInternalLabels(final int label) {
		final Iterator<HE_Vertex> vItr = vItr();
		while (vItr.hasNext()) {
			vItr.next().setInternalLabel(label);
		}
	}

	/**
	 * Return all vertex positions as an array .
	 *
	 * @return 2D array of float. First index gives vertex. Second index gives
	 *         x-,y- or z-coordinate.
	 */
	public float[][] getVerticesAsFloat() {
		final float[][] result = new float[getNumberOfVertices()][3];
		int i = 0;
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = vItr();
		while (vItr.hasNext()) {
			v = vItr.next();
			result[i][0] = v.xf();
			result[i][1] = v.yf();
			result[i][2] = v.zf();
			i++;
		}
		return result;
	}

	/**
	 * Return all vertex positions as an array .
	 *
	 * @return 2D array of double. First index gives vertex. Second index gives
	 *         x-,y- or z-coordinate.
	 */
	public double[][] getVerticesAsDouble() {
		final double[][] result = new double[getNumberOfVertices()][3];
		int i = 0;
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = vItr();
		while (vItr.hasNext()) {
			v = vItr.next();
			result[i][0] = v.xd();
			result[i][1] = v.yd();
			result[i][2] = v.zd();
			i++;
		}
		return result;
	}

	/**
	 * Return all vertex positions as an immutable List of immutable WB_Coord.
	 *
	 * @return array of WB_Coord.
	 */
	public List<WB_Coord> getVerticesAsCoord() {
		final FastTable<WB_Coord> result = new FastTable<WB_Coord>();
		;
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = vItr();
		while (vItr.hasNext()) {
			v = vItr.next();
			result.add(v);
		}
		return result.unmodifiable();
	}

	/**
	 * Vertex key to index.
	 *
	 * @return the map
	 */
	public Map<Long, Integer> getVertexKeyToIndexMap() {
		final Map<Long, Integer> map = new FastMap<Long, Integer>();
		int i = 0;
		final Iterator<HE_Vertex> vItr = vItr();
		while (vItr.hasNext()) {
			map.put(vItr.next().key(), i);
			i++;
		}
		return map;
	}

	/**
	 * Return all vertex normals.
	 *
	 * @return array of WB_Coordinate.
	 */
	public WB_Coord[] getVertexNormals() {
		final WB_Coord[] result = new WB_Coord[getNumberOfVertices()];
		int i = 0;
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = vItr();
		while (vItr.hasNext()) {
			v = vItr.next();
			result[i] = v.getVertexNormal();
			i++;
		}
		return result;
	}

	/**
	 * Get vertex normals mapped on vertex key
	 *
	 *
	 * @return
	 */
	public Map<Long, WB_Coord> getKeyedVertexNormals() {
		final Map<Long, WB_Coord> result = new FastMap<Long, WB_Coord>();
		HE_Vertex v;
		final Iterator<HE_Vertex> vItr = vItr();
		while (vItr.hasNext()) {
			v = vItr.next();
			result.put(v.key(), v.getVertexNormal());
		}
		return result;
	}

	/**
	 * Return the faces as array of vertex indices.
	 *
	 * @return 2D array of int. First index gives face. Second index gives
	 *         vertices.
	 */

	public int[][] getFacesAsInt() {
		final int[][] result = new int[getNumberOfFaces()][];
		final TLongIntMap vertexKeys = new TLongIntHashMap(10, 0.5f, -1L, -1);
		final Iterator<HE_Vertex> vItr = vItr();
		int i = 0;
		while (vItr.hasNext()) {
			vertexKeys.put(vItr.next().key(), i);
			i++;
		}
		final Iterator<HE_Face> fItr = fItr();
		HE_Halfedge he;
		HE_Face f;
		i = 0;
		while (fItr.hasNext()) {
			f = fItr.next();
			result[i] = new int[f.getFaceOrder()];
			he = f.getHalfedge();
			int j = 0;
			do {
				result[i][j] = vertexKeys.get(he.getVertex().key());
				he = he.getNextInFace();
				j++;
			} while (he != f.getHalfedge());
			i++;
		}
		return result;
	}

	/**
	 * Return all face normals.
	 *
	 * @return array of WB_Coordinate.
	 */
	public WB_Coord[] getFaceNormals() {
		final WB_Coord[] result = new WB_Coord[getNumberOfFaces()];
		int i = 0;
		HE_Face f;
		final Iterator<HE_Face> fItr = fItr();
		while (fItr.hasNext()) {
			f = fItr.next();
			result[i] = f.getFaceNormal();
			i++;
		}
		return result;
	}

	/**
	 * Return all face normals.
	 *
	 * @return FastMap of WB_Coordinate.
	 */
	public Map<Long, WB_Coord> getKeyedFaceNormals() {
		final Map<Long, WB_Coord> result = new FastMap<Long, WB_Coord>();
		HE_Face f;
		final Iterator<HE_Face> fItr = fItr();
		while (fItr.hasNext()) {
			f = fItr.next();
			result.put(f.key(), f.getFaceNormal());
		}
		return result;
	}

	/**
	 * Return all face centers.
	 *
	 * @return array of WB_Coordinate.
	 */
	public WB_Coord[] getFaceCenters() {
		final WB_Coord[] result = new WB_Coord[getNumberOfFaces()];
		int i = 0;
		HE_Face f;
		final Iterator<HE_Face> fItr = fItr();
		while (fItr.hasNext()) {
			f = fItr.next();
			result[i] = f.getFaceCenter();
			i++;
		}
		return result;
	}

	/**
	 * Return all face centers.
	 *
	 * @return FastMap of WB_Coordinate.
	 */
	public Map<Long, WB_Coord> getKeyedFaceCenters() {
		final Map<Long, WB_Coord> result = new FastMap<Long, WB_Coord>();
		HE_Face f;
		final Iterator<HE_Face> fItr = fItr();
		while (fItr.hasNext()) {
			f = fItr.next();
			result.put(f.key(), f.getFaceCenter());
		}
		return result;
	}

	/**
	 * Return the mesh as polygon soup.
	 *
	 * @return array of WB_polygon
	 *
	 */
	public WB_Polygon[] getPolygons() {
		final WB_Polygon[] result = new WB_Polygon[getNumberOfFaces()];
		final Iterator<HE_Face> fItr = fItr();
		HE_Face f;
		int i = 0;
		while (fItr.hasNext()) {
			f = fItr.next();
			result[i] = f.toPolygon();
			i++;
		}
		return result;
	}

	/**
	 * Gets the polygon list.
	 *
	 * @return the polygon list
	 */
	public List<WB_Polygon> getPolygonList() {
		final List<WB_Polygon> result = new FastTable<WB_Polygon>();
		final Iterator<HE_Face> fItr = fItr();
		HE_Face f;
		while (fItr.hasNext()) {
			f = fItr.next();
			result.add(f.toPolygon());
		}
		return result;
	}

	public int[][] getEdgesAsInt() {
		final int[][] result = new int[getNumberOfEdges()][2];
		final TLongIntMap vertexKeys = new TLongIntHashMap(10, 0.5f, -1L, -1);
		final Iterator<HE_Vertex> vItr = vItr();
		int i = 0;
		while (vItr.hasNext()) {
			vertexKeys.put(vItr.next().key(), i);
			i++;
		}
		final Iterator<HE_Halfedge> eItr = eItr();
		HE_Halfedge he;
		i = 0;
		while (eItr.hasNext()) {
			he = eItr.next();
			result[i][0] = vertexKeys.get(he.getVertex().key());
			he = he.getPair();
			result[i][1] = vertexKeys.get(he.getVertex().key());
			i++;
		}
		return result;
	}

	/**
	 * Return all edge normals.
	 *
	 * @return array of WB_Coordinate.
	 */
	public WB_Coord[] getEdgeNormals() {
		final WB_Coord[] result = new WB_Coord[getNumberOfEdges()];
		int i = 0;
		HE_Halfedge e;
		final Iterator<HE_Halfedge> eItr = eItr();
		while (eItr.hasNext()) {
			e = eItr.next();
			result[i] = e.getEdgeNormal();
			i++;
		}
		return result;
	}

	/**
	 * Return all edge normals.
	 *
	 * @return FastMap of WB_Coordinate.
	 */
	public Map<Long, WB_Coord> getKeyedEdgeNormals() {
		final Map<Long, WB_Coord> result = new FastMap<Long, WB_Coord>();
		HE_Halfedge e;
		final Iterator<HE_Halfedge> eItr = eItr();
		while (eItr.hasNext()) {
			e = eItr.next();
			result.put(e.key(), e.getEdgeNormal());
		}
		return result;
	}

	/**
	 * Return all edge centers.
	 *
	 * @return array of WB_Coordinate.
	 */
	public WB_Coord[] getEdgeCenters() {
		final WB_Coord[] result = new WB_Coord[getNumberOfEdges()];
		int i = 0;
		HE_Halfedge e;
		final Iterator<HE_Halfedge> eItr = eItr();
		while (eItr.hasNext()) {
			e = eItr.next();
			result[i] = e.getHalfedgeCenter();
			i++;
		}
		return result;
	}

	/**
	 * Return all edge centers.
	 *
	 * @return FastMap of WB_Coordinate.
	 */
	public Map<Long, WB_Coord> getKeyedEdgeCenters() {
		final Map<Long, WB_Coord> result = new FastMap<Long, WB_Coord>();
		HE_Halfedge e;
		final Iterator<HE_Halfedge> eItr = eItr();
		while (eItr.hasNext()) {
			e = eItr.next();
			result.put(e.key(), e.getHalfedgeCenter());
		}
		return result;
	}

	/**
	 * Gets the segments.
	 *
	 * @return the segments
	 */
	public WB_Segment[] getSegments() {
		final WB_Segment[] result = new WB_Segment[getNumberOfEdges()];
		final Iterator<HE_Halfedge> eItr = eItr();
		HE_Halfedge e;
		int i = 0;
		while (eItr.hasNext()) {
			e = eItr.next();
			result[i] = new WB_Segment(e.getVertex(), e.getEndVertex());
			i++;
		}
		return result;
	}

	/**
	 *
	 *
	 * @return
	 */
	public HE_Path[] getBoundaryAsPath() {
		final List<HE_Halfedge> boundaryhes = getBoundaryLoopHalfedges();
		final HE_Path[] result = new HE_Path[boundaryhes.size()];
		for (int i = 0; i < boundaryhes.size(); i++) {
			result[i] = new HE_Path(boundaryhes.get(i));
		}
		return result;
	}

	/**
	 *
	 *
	 * @return
	 */
	public List<WB_Polygon> getBoundaryAsPolygons() {
		final List<WB_Polygon> polygons = new FastTable<WB_Polygon>();
		final List<HE_Halfedge> halfedges = getBoundaryHalfedges();
		final List<HE_Halfedge> loop = new FastTable<HE_Halfedge>();
		final List<WB_Coord> points = new FastTable<WB_Coord>();
		while (halfedges.size() > 0) {
			points.clear();
			loop.clear();
			HE_Halfedge he = halfedges.get(0);
			do {
				loop.add(he);
				points.add(he.getVertex());
				he = he.getNextInFace();
				if (loop.contains(he)) {
					break;
				}
			} while (he != halfedges.get(0));
			polygons.add(gf.createSimplePolygon(points));
			halfedges.removeAll(loop);
		}
		return polygons;
	}

	/**
	 * Collect all boundary halfedges.
	 *
	 * @return boundary halfedges
	 */
	public List<HE_Halfedge> getBoundaryHalfedges() {
		final List<HE_Halfedge> boundaryHalfedges = new FastTable<HE_Halfedge>();
		HE_Halfedge he;
		final Iterator<HE_Halfedge> heItr = heItr();
		while (heItr.hasNext()) {
			he = heItr.next();
			if (he.getFace() == null) {
				boundaryHalfedges.add(he);
			}
		}
		return boundaryHalfedges;
	}

	/**
	 * Collect all boundary vertices.
	 *
	 * @return boundary vertices
	 */
	public List<HE_Vertex> getBoundaryVertices() {
		final List<HE_Vertex> boundaryVertices = new FastTable<HE_Vertex>();
		HE_Halfedge he;
		final Iterator<HE_Halfedge> heItr = heItr();
		while (heItr.hasNext()) {
			he = heItr.next();
			if (he.getFace() == null) {
				boundaryVertices.add(he.getVertex());
			}
		}
		return boundaryVertices;
	}

	/**
	 *
	 *
	 * @return
	 */
	public List<HE_Halfedge> getBoundaryLoopHalfedges() {
		final List<HE_Halfedge> hes = new FastTable<HE_Halfedge>();
		final List<HE_Halfedge> halfedges = getBoundaryHalfedges();
		final List<HE_Halfedge> loop = new FastTable<HE_Halfedge>();
		while (halfedges.size() > 0) {
			loop.clear();
			HE_Halfedge he = halfedges.get(0);
			hes.add(he);
			do {
				loop.add(he);
				he = he.getNextInFace();
				if (loop.contains(he)) {
					break;
				}
			} while (he != halfedges.get(0));
			halfedges.removeAll(loop);
		}
		return hes;
	}

	/**
	 * Clean all mesh elements not used by any faces.
	 *
	 * @return self
	 */
	public HE_MeshStructure cleanUnusedElementsByFace() {

		return HET_MeshOp.cleanUnusedElementsByFace(this);
	}

	/**
	 *
	 *
	 * @param vertices
	 * @param loop
	 * @return
	 */
	public HE_Path createPathFromIndices(final int[] vertices, final boolean loop) {
		final List<HE_Halfedge> halfedges = new FastTable<HE_Halfedge>();
		if (vertices.length > 1) {
			HE_Halfedge he;
			for (int i = 0; i < vertices.length - 1; i++) {
				he = getHalfedgeFromTo(getVertexWithIndex(vertices[i]), getVertexWithIndex(vertices[i + 1]));
				if (he == null) {
					// throw new IllegalArgumentException(
					// "Two vertices " + vertices[i] + " and " + vertices[i + 1]
					// + " in path are not connected.");
				} else {

					halfedges.add(he);
				}
			}
			if (loop) {
				he = getHalfedgeFromTo(getVertexWithIndex(vertices[vertices.length - 1]),
						getVertexWithIndex(vertices[0]));
				if (he == null) {
					throw new IllegalArgumentException("Vertices " + vertices[vertices.length - 1] + " and "
							+ vertices[0] + " in path are not connected: path is not a loop.");
				}
			}
		}
		final HE_Path path = new HE_Path(halfedges, loop);
		return path;
	}

	/**
	 *
	 *
	 * @param partition1
	 * @param partition2
	 * @return
	 */
	boolean isNeighbor(final HE_RAS<HE_Face> partition1, final HE_RAS<HE_Face> partition2) {
		HE_Halfedge he1, he2;
		HE_Vertex v1;
		HE_FaceHalfedgeInnerCirculator heitr1, heitr2;
		for (final HE_Face f1 : partition1) {
			heitr1 = new HE_FaceHalfedgeInnerCirculator(f1);
			while (heitr1.hasNext()) {
				he1 = heitr1.next();
				if (he1.getPair() == null) {
					v1 = he1.getNextInFace().getVertex();
					for (final HE_Face f2 : partition2) {
						heitr2 = new HE_FaceHalfedgeInnerCirculator(f2);
						while (heitr2.hasNext()) {
							he2 = heitr2.next();
							if (he2.getPair() == null && he2.getNextInFace().getVertex() == v1) {
								return true;
							}
						}
					}
				}
			}
		}
		return false;
	}

	/**
	 * Cycle halfedges.
	 *
	 * @param halfedges
	 *            halfedges to cycle
	 */
	public void cycleHalfedges(final List<HE_Halfedge> halfedges) {
		HE_Halfedge he;
		final int n = halfedges.size();
		if (n > 0) {
			for (int j = 0; j < n - 1; j++) {
				he = halfedges.get(j);
				setNext(he, halfedges.get(j + 1));
			}
			he = halfedges.get(n - 1);
			setNext(he, halfedges.get(0));
		}
	}

	/**
	 * Cycle halfedges.
	 *
	 * @param halfedges
	 *            halfedges to cycle
	 */
	public void cycleHalfedgesReverse(final List<HE_Halfedge> halfedges) {
		HE_Halfedge he;
		final int n = halfedges.size();
		if (n > 0) {
			he = halfedges.get(0);
			setNext(he, halfedges.get(n - 1));

			for (int j = 1; j < n; j++) {
				he = halfedges.get(j);
				setNext(he, halfedges.get(j - 1));

			}
		}
	}

	public void orderHalfedges(final List<HE_Halfedge> halfedges) {
		HE_Halfedge he;
		final int n = halfedges.size();
		if (n > 0) {
			for (int j = 0; j < n - 1; j++) {
				he = halfedges.get(j);
				setNext(he, halfedges.get(j + 1));
			}
		}
	}

	/**
	 * Cycle halfedges.
	 *
	 * @param halfedges
	 *            halfedges to cycle
	 */
	public void orderHalfedgesReverse(final List<HE_Halfedge> halfedges) {
		HE_Halfedge he;
		final int n = halfedges.size();
		if (n > 0) {

			for (int j = 1; j < n; j++) {
				he = halfedges.get(j);
				setNext(he, halfedges.get(j - 1));

			}
		}
	}

	/**
	 * Uncap halfedges.
	 */
	public void uncapBoundaryHalfedges() {
		tracker.setStatus(this, "Uncapping boundary halfedges.", +1);
		WB_ProgressCounter counter = new WB_ProgressCounter(getNumberOfHalfedges(), 10);
		tracker.setStatus(this, "Detecting and uncapping boundary edges.", 0);

		List<HE_Halfedge> halfedges = getHalfedges();
		final HE_RAS<HE_Halfedge> keep = new HE_RASTrove<HE_Halfedge>();
		for (HE_Halfedge he : halfedges) {

			if (he.getFace() == null) {
				setHalfedge(he.getVertex(), he.getNextInVertex());
				clearPair(he);
			} else {
				keep.add(he);
			}
			counter.increment();
		}
		clearHalfedges();
		addHalfedges(keep);
		tracker.setStatus(this, "Removing outer boundary halfedges.", -1);

	}

	/**
	 * Cap all remaining unpaired halfedges. Only use after pairHalfedges();
	 */
	public void capHalfedges() {

		tracker.setStatus(this, "Capping unpaired halfedges.", +1);
		final List<HE_Halfedge> unpairedHalfedges = getUnpairedHalfedges();
		final int nuh = unpairedHalfedges.size();
		final HE_Halfedge[] newHalfedges = new HE_Halfedge[nuh];
		HE_Halfedge he1, he2;
		WB_ProgressCounter counter = new WB_ProgressCounter(nuh, 10);
		tracker.setStatus(this, "Capping unpaired halfedges.", counter);
		for (int i = 0; i < nuh; i++) {
			he1 = unpairedHalfedges.get(i);
			he2 = new HE_Halfedge();
			setVertex(he2, he1.getNextInFace().getVertex());
			setPair(he1, he2);
			newHalfedges[i] = he2;
			add(he2);
			counter.increment();
		}
		counter = new WB_ProgressCounter(nuh, 10);
		tracker.setStatus(this, "Cycling new halfedges.", counter);
		for (int i = 0; i < nuh; i++) {
			he1 = newHalfedges[i];
			if (he1.getNextInFace() == null) {
				for (int j = 0; j < nuh; j++) {
					he2 = newHalfedges[j];
					if (!he2.isVisited()) {
						if (he2.getVertex() == he1.getPair().getVertex()) {
							setNext(he1, he2);
							he2.setVisited();
							break;
						}
					}
				}
			}
			counter.increment();
		}
		tracker.setStatus(this, "Processed unpaired halfedges.", -1);
	}

	/**
	 * Iterate through all halfedges and reset the halfedge link to its vertex
	 * to itself. v=he.getVertex() v.setHalfedge(he)
	 */
	public void fixHalfedgeVertexAssignment() {
		final Iterator<HE_Halfedge> heItr = heItr();
		HE_Halfedge he;
		while (heItr.hasNext()) {
			he = heItr.next();
			setHalfedge(he.getVertex(), he);
		}
	}

	/**
	 * Iterate through all halfedges and reset the halfedge link to its face to
	 * itself. f=he.getFace() f.setHalfedge(he)
	 */
	public void fixHalfedgeFaceAssignment() {
		final Iterator<HE_Halfedge> heItr = heItr();
		HE_Halfedge he;
		while (heItr.hasNext()) {
			he = heItr.next();
			if (he.getFace() != null) {
				setHalfedge(he.getFace(), he);
			}
		}
	}

	/**
	 * Collect all unpaired halfedges.
	 *
	 * @return the unpaired halfedges
	 */
	public List<HE_Halfedge> getUnpairedHalfedges() {
		final List<HE_Halfedge> halfedges = new FastTable<HE_Halfedge>();
		halfedges.addAll(this.unpairedHalfedges);
		return halfedges;
	}

	/**
	 * Try to pair all unpaired halfedges.
	 *
	 * @return
	 */
	public List<HE_Halfedge> pairHalfedgesOnePass() {
		tracker.setStatus(this, "Pairing halfedges.", +1);
		class VertexInfo {
			FastTable<HE_Halfedge> out;
			FastTable<HE_Halfedge> in;

			VertexInfo() {
				out = new FastTable<HE_Halfedge>();
				in = new FastTable<HE_Halfedge>();
			}
		}
		final TLongObjectMap<VertexInfo> vertexLists = new TLongObjectHashMap<VertexInfo>(1024, 0.5f, -1L);
		tracker.setStatus(this, "Collecting unpaired halfedges.", 0);
		final List<HE_Halfedge> unpairedHalfedges = getUnpairedHalfedges();
		HE_Vertex v;
		VertexInfo vi;
		WB_ProgressCounter counter = new WB_ProgressCounter(unpairedHalfedges.size(), 10);
		tracker.setStatus(this, "Classifying unpaired halfedges.", counter);
		for (final HE_Halfedge he : unpairedHalfedges) {
			v = he.getVertex();
			vi = vertexLists.get(v.key());
			if (vi == null) {
				vi = new VertexInfo();
				vertexLists.put(v.key(), vi);
			}
			vi.out.add(he);
			v = he.getNextInFace().getVertex();
			vi = vertexLists.get(v.key());
			if (vi == null) {
				vi = new VertexInfo();
				vertexLists.put(v.key(), vi);
			}
			vi.in.add(he);
			counter.increment();
		}
		HE_Halfedge he;
		HE_Halfedge he2;
		counter = new WB_ProgressCounter(vertexLists.size(), 10);
		tracker.setStatus(this, "Pairing unpaired halfedges per vertex.", counter);
		final TLongObjectIterator<VertexInfo> vitr = vertexLists.iterator();
		VertexInfo vInfo;
		final List<HE_Halfedge> mismatchedHalfedges = new FastTable<HE_Halfedge>();
		while (vitr.hasNext()) {
			vitr.advance();
			vInfo = vitr.value();
			for (int i = 0; i < vInfo.out.size(); i++) {
				he = vInfo.out.get(i);
				if (he.getPair() == null) {
					for (int j = 0; j < vInfo.in.size(); j++) {
						he2 = vInfo.in.get(j);
						if (he2.getPair() == null) {
							if (he.getVertex() == he2.getNextInFace().getVertex()
									&& he2.getVertex() == he.getNextInFace().getVertex()) {
								setPair(he, he2);
								break;
							}
						}
					}
					for (int j = 0; j < vInfo.out.size(); j++) {
						he2 = vInfo.out.get(j);
						if (he2 != he && he2.getPair() == null) {
							if (he.getNextInFace().getVertex() == he2.getNextInFace().getVertex()) {
								mismatchedHalfedges.add(he);
								mismatchedHalfedges.add(he2);
								break;
							}
						}
					}
				}
			}
			counter.increment();
		}
		tracker.setStatus(this, "Processed unpaired halfedges.", -1);
		return mismatchedHalfedges;
	}

	/**
	 *
	 */
	public void pairHalfedges() {
		pairHalfedgesOnePass();
	}

	/**
	 * Pair halfedges.
	 *
	 * @param unpairedHalfedges
	 *            the unpaired halfedges
	 */
	public void pairHalfedges(final List<HE_Halfedge> unpairedHalfedges) {
		tracker.setStatus(this, "Pairing halfedges.", +1);
		class VertexInfo {
			FastTable<HE_Halfedge> out;
			FastTable<HE_Halfedge> in;

			VertexInfo() {
				out = new FastTable<HE_Halfedge>();
				in = new FastTable<HE_Halfedge>();
			}
		}
		final TLongObjectMap<VertexInfo> vertexLists = new TLongObjectHashMap<VertexInfo>(1024, 0.5f, -1L);
		HE_Vertex v;
		VertexInfo vi;
		WB_ProgressCounter counter = new WB_ProgressCounter(unpairedHalfedges.size(), 10);
		tracker.setStatus(this, "Classifying unpaired halfedges.", counter);
		for (final HE_Halfedge he : unpairedHalfedges) {
			v = he.getVertex();
			vi = vertexLists.get(v.key());
			if (vi == null) {
				vi = new VertexInfo();
				vertexLists.put(v.key(), vi);
			}
			vi.out.add(he);
			v = he.getNextInFace().getVertex();
			vi = vertexLists.get(v.key());
			if (vi == null) {
				vi = new VertexInfo();
				vertexLists.put(v.key(), vi);
			}
			vi.in.add(he);
			counter.increment();
		}
		HE_Halfedge he;
		HE_Halfedge he2;
		counter = new WB_ProgressCounter(vertexLists.size(), 10);
		tracker.setStatus(this, "Pairing unpaired halfedges per vertex.", counter);
		final TLongObjectIterator<VertexInfo> vitr = vertexLists.iterator();
		VertexInfo vInfo;
		while (vitr.hasNext()) {
			vitr.advance();
			vInfo = vitr.value();
			for (int i = 0; i < vInfo.out.size(); i++) {
				he = vInfo.out.get(i);
				if (he.getPair() == null) {
					for (int j = 0; j < vInfo.in.size(); j++) {
						he2 = vInfo.in.get(j);
						if (he2.getPair() == null) {
							if (he.getVertex() == he2.getNextInFace().getVertex()
									&& he2.getVertex() == he.getNextInFace().getVertex()) {
								setPair(he, he2);
								break;
							}
						}
					}
					for (int j = 0; j < vInfo.out.size(); j++) {
						he2 = vInfo.out.get(j);
						if (he2 != he && he2.getPair() == null) {
							if (he.getNextInFace().getVertex() == he2.getNextInFace().getVertex()) {
								System.out.println("Two identical halfedges found!");
								break;
							}
						}
					}
				}
			}
			counter.increment();
		}
		tracker.setStatus(this, "Processed unpaired halfedges.", -1);
	}

	/**
	 * Return a halfedge from vertex v0 to vertex v1. If no such halfedge exists
	 * return null.
	 *
	 * @param v0
	 * @param v1
	 * @return
	 */
	public HE_Halfedge getHalfedgeFromTo(final HE_Vertex v0, final HE_Vertex v1) {
		final List<HE_Halfedge> hes = v0.getHalfedgeStar();
		for (final HE_Halfedge he : hes) {
			if (he.getEndVertex() == v1) {
				return he;
			}
		}
		return null;
	}

	/**
	 *
	 *
	 * @param v1
	 * @param v2
	 * @return
	 */
	public List<HE_Face> getSharedFaces(final HE_Vertex v1, final HE_Vertex v2) {
		final List<HE_Face> result = v1.getFaceStar();
		final List<HE_Face> compare = v2.getFaceStar();
		final Iterator<HE_Face> it = result.iterator();
		while (it.hasNext()) {
			if (!compare.contains(it.next())) {
				it.remove();
			}
		}
		return result;
	}

	/**
	 *
	 *
	 * @return
	 */
	public int getNumberOfBoundaryComponents() {
		return getBoundaryLoopHalfedges().size();

	}

	/**
	 *
	 *
	 * @return
	 */
	public int getEulerCharacteristic() {
		return getNumberOfVertices() - getNumberOfEdges() + getNumberOfFaces();
	}

	/**
	 *
	 *
	 * @return
	 */
	public int getGenus() {
		return (2 - getEulerCharacteristic() - getNumberOfBoundaryComponents()) / 2;
	}

	/**
		 *
		 */
	public void clearVisitedElements() {
		final HE_FaceIterator fitr = fItr();
		while (fitr.hasNext()) {
			fitr.next().clearVisited();
		}
		final HE_VertexIterator vitr = vItr();
		while (vitr.hasNext()) {
			vitr.next().clearVisited();
		}
		final HE_HalfedgeIterator heitr = heItr();
		while (heitr.hasNext()) {
			heitr.next().clearVisited();
		}
	}

	public void update() {
		updateFaces();
	}

	public void updateFaces() {
		HE_FaceIterator fItr = fItr();
		while (fItr.hasNext()) {
			fItr.next().update();

		}
	}

}
