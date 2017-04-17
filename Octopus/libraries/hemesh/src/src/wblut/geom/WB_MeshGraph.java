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

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.PriorityQueue;

import javolution.util.FastTable;
import wblut.hemesh.HE_Mesh;

/**
 *
 */
public class WB_MeshGraph {
	/**
	 *
	 */
	private final WB_GraphVertex[] vertices;
	/**
	 *
	 */
	private int lastSource;

	/**
	 *
	 *
	 * @param mesh
	 */
	public WB_MeshGraph(final WB_Mesh mesh) {
		vertices = new WB_GraphVertex[mesh.getNumberOfVertices()];
		for (int i = 0; i < mesh.getNumberOfVertices(); i++) {
			vertices[i] = new WB_GraphVertex(i, mesh.getVertex(i));
		}
		final int[][] meshedges = mesh.getEdgesAsInt();
		WB_Coord p0;
		WB_Coord p1;
		WB_GraphVertex v0;
		WB_GraphVertex v1;
		double d;
		for (int i = 0; i < meshedges.length; i++) {
			if (meshedges[i][0] != meshedges[i][1]) {
				p0 = mesh.getVertex(meshedges[i][0]);
				p1 = mesh.getVertex(meshedges[i][1]);
				d = WB_GeometryOp3D.getDistance3D(p0, p1);
				v0 = vertices[meshedges[i][0]];
				v1 = vertices[meshedges[i][1]];
				v0.neighbors.add(new WB_GraphEdge(v1, d));
				v1.neighbors.add(new WB_GraphEdge(v0, d));
			}
		}
		lastSource = -1;
	}

	/**
	 *
	 *
	 * @param mesh
	 * @param offset
	 */
	public WB_MeshGraph(final WB_Mesh mesh, final double offset) {
		vertices = new WB_GraphVertex[mesh.getNumberOfVertices()];
		for (int i = 0; i < mesh.getNumberOfVertices(); i++) {
			vertices[i] = new WB_GraphVertex(i,
					new WB_Point(mesh.getVertex(i)).addMulSelf(offset, mesh.getVertexNormal(i)));
		}
		final int[][] meshedges = mesh.getEdgesAsInt();
		WB_Coord p0;
		WB_Coord p1;
		WB_GraphVertex v0;
		WB_GraphVertex v1;
		double d;
		for (int i = 0; i < meshedges.length; i++) {
			if (meshedges[i][0] != meshedges[i][1]) {
				p0 = mesh.getVertex(meshedges[i][0]);
				p1 = mesh.getVertex(meshedges[i][1]);
				d = WB_GeometryOp3D.getDistance3D(p0, p1);
				v0 = vertices[meshedges[i][0]];
				v1 = vertices[meshedges[i][1]];
				v0.neighbors.add(new WB_GraphEdge(v1, d));
				v1.neighbors.add(new WB_GraphEdge(v0, d));
			}
		}
		lastSource = -1;
	}

	/**
	 *
	 *
	 * @param mesh
	 */
	public WB_MeshGraph(final HE_Mesh mesh) {
		vertices = new WB_GraphVertex[mesh.getNumberOfVertices()];
		for (int i = 0; i < mesh.getNumberOfVertices(); i++) {
			vertices[i] = new WB_GraphVertex(i, mesh.getVertex(i));
		}
		final int[][] meshedges = mesh.getEdgesAsInt();
		WB_Coord p0;
		WB_Coord p1;
		WB_GraphVertex v0;
		WB_GraphVertex v1;
		double d;
		for (int i = 0; i < meshedges.length; i++) {
			if (meshedges[i][0] != meshedges[i][1]) {
				p0 = mesh.getVertex(meshedges[i][0]);
				p1 = mesh.getVertex(meshedges[i][1]);
				d = WB_GeometryOp3D.getDistance3D(p0, p1);
				v0 = vertices[meshedges[i][0]];
				v1 = vertices[meshedges[i][1]];
				v0.neighbors.add(new WB_GraphEdge(v1, d));
				v1.neighbors.add(new WB_GraphEdge(v0, d));
			}
		}
		lastSource = -1;
	}

	/**
	 *
	 *
	 * @param mesh
	 * @param offset
	 */
	public WB_MeshGraph(final HE_Mesh mesh, final double offset) {
		vertices = new WB_GraphVertex[mesh.getNumberOfVertices()];
		for (int i = 0; i < mesh.getNumberOfVertices(); i++) {
			vertices[i] = new WB_GraphVertex(i,
					new WB_Point(mesh.getVertex(i)).addMulSelf(offset, mesh.getVertexNormal(i)));
		}
		final int[][] meshedges = mesh.getEdgesAsInt();
		WB_Coord p0;
		WB_Coord p1;
		WB_GraphVertex v0;
		WB_GraphVertex v1;
		double d;
		for (int i = 0; i < meshedges.length; i++) {
			if (meshedges[i][0] != meshedges[i][1]) {
				p0 = mesh.getVertex(meshedges[i][0]);
				p1 = mesh.getVertex(meshedges[i][1]);
				d = WB_GeometryOp3D.getDistance3D(p0, p1);
				v0 = vertices[meshedges[i][0]];
				v1 = vertices[meshedges[i][1]];
				v0.neighbors.add(new WB_GraphEdge(v1, d));
				v1.neighbors.add(new WB_GraphEdge(v0, d));
			}
		}
		lastSource = -1;
	}

	public WB_MeshGraph(final List<? extends WB_Coord> points, final WB_Triangulation3D triangulation) {
		vertices = new WB_GraphVertex[points.size()];
		for (int i = 0; i < points.size(); i++) {
			vertices[i] = new WB_GraphVertex(i, points.get(i));
		}
		final int[] meshedges = triangulation.getEdges();
		WB_Coord p0;
		WB_Coord p1;
		WB_GraphVertex v0;
		WB_GraphVertex v1;
		double d;
		for (int i = 0; i < meshedges.length; i += 2) {
			if (meshedges[i] != meshedges[i + 1]) {
				p0 = points.get(meshedges[i]);
				p1 = points.get(meshedges[i + 1]);
				d = WB_GeometryOp3D.getDistance3D(p0, p1);
				v0 = vertices[meshedges[i]];
				v1 = vertices[meshedges[i + 1]];
				v0.neighbors.add(new WB_GraphEdge(v1, d));
				v1.neighbors.add(new WB_GraphEdge(v0, d));
			}
		}
		lastSource = -1;
	}

	public WB_MeshGraph(final WB_Coord[] points, final WB_Triangulation3D triangulation) {
		vertices = new WB_GraphVertex[points.length];
		for (int i = 0; i < points.length; i++) {
			vertices[i] = new WB_GraphVertex(i, points[i]);
		}
		final int[] meshedges = triangulation.getEdges();
		WB_Coord p0;
		WB_Coord p1;
		WB_GraphVertex v0;
		WB_GraphVertex v1;
		double d;
		for (int i = 0; i < meshedges.length; i += 2) {
			if (meshedges[i] != meshedges[i + 1]) {
				p0 = points[meshedges[i]];
				p1 = points[meshedges[i + 1]];
				d = WB_GeometryOp3D.getDistance3D(p0, p1);
				v0 = vertices[meshedges[i]];
				v1 = vertices[meshedges[i + 1]];
				v0.neighbors.add(new WB_GraphEdge(v1, d));
				v1.neighbors.add(new WB_GraphEdge(v0, d));
			}
		}
		lastSource = -1;
	}

	public WB_MeshGraph(final WB_Coord[] points, final WB_Triangulation2D triangulation) {
		vertices = new WB_GraphVertex[points.length];
		for (int i = 0; i < points.length; i++) {
			vertices[i] = new WB_GraphVertex(i, points[i]);
		}
		final int[] meshedges = triangulation.getEdges();
		WB_Coord p0;
		WB_Coord p1;
		WB_GraphVertex v0;
		WB_GraphVertex v1;
		double d;
		for (int i = 0; i < meshedges.length; i += 2) {
			if (meshedges[i] != meshedges[i + 1]) {
				p0 = points[meshedges[i]];
				p1 = points[meshedges[i + 1]];
				d = WB_GeometryOp3D.getDistance3D(p0, p1);
				v0 = vertices[meshedges[i]];
				v1 = vertices[meshedges[i + 1]];
				v0.neighbors.add(new WB_GraphEdge(v1, d));
				v1.neighbors.add(new WB_GraphEdge(v0, d));
			}
		}
		lastSource = -1;
	}

	/**
	 *
	 *
	 * @param i
	 * @return
	 */
	public int getVertexIndex(final int i) {
		return vertices[i].index;
	}

	/**
	 *
	 *
	 * @param i
	 */
	public void computePathsToVertex(final int i) {
		final WB_GraphVertex source = vertices[i];
		for (int j = 0; j < vertices.length; j++) {
			vertices[j].reset();
		}
		source.minWeight = 0.;
		final PriorityQueue<WB_GraphVertex> vertexQueue = new PriorityQueue<WB_GraphVertex>();
		vertexQueue.add(source);
		while (!vertexQueue.isEmpty()) {
			final WB_GraphVertex u = vertexQueue.poll();
			// Visit each edge exiting u
			for (final WB_GraphEdge e : u.neighbors) {
				final WB_GraphVertex v = e.target;
				final double weight = e.weight;
				final double distanceThroughU = u.minWeight + weight;
				if (distanceThroughU < v.minWeight) {
					vertexQueue.remove(v);
					v.minWeight = distanceThroughU;
					v.previous = u;
					vertexQueue.add(v);
				}
			}
		}
		lastSource = i;
	}

	/**
	 *
	 *
	 * @param source
	 * @param target
	 * @return
	 */
	public int[] getShortestPathBetweenVertices(final int source, final int target) {
		if (source != lastSource) {
			computePathsToVertex(source);
		}
		if (source == target) {
			return new int[] { source };
		}
		final List<WB_GraphVertex> path = new ArrayList<WB_GraphVertex>();
		for (WB_GraphVertex vertex = vertices[target]; vertex != null; vertex = vertex.previous) {
			path.add(vertex);
		}
		Collections.reverse(path);
		final int[] result = new int[path.size()];
		for (int i = 0; i < path.size(); i++) {
			result[i] = path.get(i).index;
		}
		return result;
	}

	/**
	 *
	 *
	 * @param i
	 * @return
	 */
	public WB_Frame getFrame(final int i) {
		final WB_Frame frame = new WB_Frame();
		computePathsToVertex(i);
		for (final WB_GraphVertex v : vertices) {
			frame.addNode(v.x, v.y, v.z, 0);
		}
		for (final WB_GraphVertex v : vertices) {
			final int[] path = getShortestPathBetweenVertices(i, v.index);
			for (int j = 0; j < path.length - 1; j++) {
				frame.nodes.get(path[j]).value = Math.max(frame.nodes.get(path[j]).value, 1.0 - j * 1.0 / path.length);
				frame.addStrut(path[j], path[j + 1]);
			}
			frame.nodes.get(path[path.length - 1]).value = Math.max(frame.nodes.get(path[path.length - 1]).value,
					1.0 / path.length);
		}
		return frame;
	}

	/**
	 *
	 *
	 * @param i
	 * @param maxnodes
	 * @return
	 */
	public WB_Frame getFrame(final int i, final int maxnodes) {
		final WB_Frame frame = new WB_Frame();
		computePathsToVertex(i);
		for (final WB_GraphVertex v : vertices) {
			frame.addNode(v.x, v.y, v.z, 0);
		}
		for (final WB_GraphVertex v : vertices) {
			final int[] path = getShortestPathBetweenVertices(i, v.index);
			final int nodes = Math.min(maxnodes, path.length);
			for (int j = 0; j < nodes - 1; j++) {
				frame.nodes.get(path[j]).value = Math.max(frame.nodes.get(path[j]).value, 1.0 - j * 1.0 / nodes);
				frame.addStrut(path[j], path[j + 1]);
			}
			frame.nodes.get(path[nodes - 1]).value = Math.max(frame.nodes.get(path[nodes - 1]).value, 1.0 / nodes);
		}
		return frame;
	}

	public WB_Frame getFrame(final int i, final int maxnodes, final double offset) {
		final WB_Frame frame = new WB_Frame();
		computePathsToVertex(i);
		for (final WB_GraphVertex v : vertices) {
			frame.addNode(v.x, v.y, v.z, 0);
		}
		for (final WB_GraphVertex v : vertices) {
			final int[] path = getShortestPathBetweenVertices(i, v.index);
			final int nodes = Math.min(maxnodes, path.length);
			for (int j = 0; j < nodes - 1; j++) {
				frame.nodes.get(path[j]).value = Math.max(frame.nodes.get(path[j]).value,
						1.0 - j * 1.0 / nodes + offset);
				frame.addStrut(path[j], path[j + 1]);
			}
			frame.nodes.get(path[nodes - 1]).value = Math.max(frame.nodes.get(path[nodes - 1]).value,
					1.0 / nodes + offset);
		}
		return frame;
	}

	/**
	 *
	 *
	 * @param i
	 * @param maxnodes
	 * @param cuttail
	 * @return
	 */
	public WB_Frame getFrame(final int i, final int maxnodes, final int cuttail) {
		final WB_Frame frame = new WB_Frame();
		computePathsToVertex(i);
		for (final WB_GraphVertex v : vertices) {
			frame.addNode(v.x, v.y, v.z, 0);
		}
		for (final WB_GraphVertex v : vertices) {
			final int[] path = getShortestPathBetweenVertices(i, v.index);
			final int nodes = Math.min(maxnodes, path.length - cuttail);
			if (nodes <= 1) {
				continue;
			}
			for (int j = 0; j < nodes - 1; j++) {
				frame.nodes.get(path[j]).value = Math.max(frame.nodes.get(path[j]).value, 1.0 - j * 1.0 / nodes);
				frame.addStrut(path[j], path[j + 1]);
			}
			frame.nodes.get(path[nodes - 1]).value = Math.max(frame.nodes.get(path[nodes - 1]).value, 1.0 / nodes);
		}
		return frame;
	}

	/**
	 *
	 *
	 * @param args
	 */
	public static void main(final String[] args) {
		final WB_Geodesic geo = new WB_Geodesic(1.0, 2, 0, WB_Geodesic.ICOSAHEDRON);
		WB_MeshGraph graph = new WB_MeshGraph(geo.create());
		for (final WB_GraphVertex v : graph.vertices) {
			final int[] path = graph.getShortestPathBetweenVertices(5, v.index);
			System.out.println("Distance to " + v + ": " + v.minWeight);
			System.out.print("Path: ");
			for (int i = 0; i < path.length - 1; i++) {
				System.out.print(path[i] + "->");
			}
			System.out.println(path[path.length - 1] + ".");
		}
		final HE_Mesh mesh = new HE_Mesh(geo.create());
		mesh.smooth();
		graph = new WB_MeshGraph(mesh);
		for (final WB_GraphVertex v : graph.vertices) {
			final int[] path = graph.getShortestPathBetweenVertices(0, v.index);
			System.out.println("Distance to " + v + ": " + v.minWeight);
			System.out.print("Path: ");
			for (int i = 0; i < path.length - 1; i++) {
				System.out.print(path[i] + "->");
			}
			System.out.println(path[path.length - 1] + ".");
		}
		for (final WB_GraphVertex v : graph.vertices) {
			final int[] path = graph.getShortestPathBetweenVertices(5, v.index);
			System.out.println("Distance to " + v + ": " + v.minWeight);
			System.out.print("Path: ");
			for (int i = 0; i < path.length - 1; i++) {
				System.out.print(path[i] + "->");
			}
			System.out.println(path[path.length - 1] + ".");
		}
	}

	/**
	 *
	 */
	public class WB_GraphVertex implements Comparable<WB_GraphVertex> {
		/**
		 *
		 */
		public final int index;
		/**
		 *
		 */
		public List<WB_GraphEdge> neighbors;
		/**
		 *
		 */
		public double minWeight = Double.POSITIVE_INFINITY;
		/**
		 *
		 */
		public WB_GraphVertex previous;
		/**
		 *
		 */
		public double x, y, z;

		/**
		 *
		 *
		 * @param id
		 * @param pos
		 */
		public WB_GraphVertex(final int id, final WB_Coord pos) {
			index = id;
			neighbors = new FastTable<WB_GraphEdge>();
			x = pos.xd();
			y = pos.yd();
			z = pos.zd();
		}

		/*
		 * (non-Javadoc)
		 *
		 * @see java.lang.Object#toString()
		 */
		@Override
		public String toString() {
			return "Vertex " + index;
		}

		/*
		 * (non-Javadoc)
		 *
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(final WB_GraphVertex other) {
			return Double.compare(minWeight, other.minWeight);
		}

		/**
		 *
		 */
		public void reset() {
			minWeight = Double.POSITIVE_INFINITY;
			previous = null;
		}
	}

	/**
	 *
	 */
	public class WB_GraphEdge {
		/**
		 *
		 */
		public final WB_GraphVertex target;
		/**
		 *
		 */
		public final double weight;

		/**
		 *
		 *
		 * @param target
		 * @param weight
		 */
		public WB_GraphEdge(final WB_GraphVertex target, final double weight) {
			this.target = target;
			this.weight = weight;
		}
	}
}
