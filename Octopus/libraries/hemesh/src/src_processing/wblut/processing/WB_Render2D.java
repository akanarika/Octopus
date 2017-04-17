/*
 * This file is part of HE_Mesh, a library for creating and manipulating meshes.
 * It is dedicated to the public domain. To the extent possible under law,
 * I , Frederik Vanhoutte, have waived all copyright and related or neighboring
 * rights.
 *
 * This work is published from Belgium. (http://creativecommons.org/publicdomain/zero/1.0/)
 *
 */
package wblut.processing;

import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import processing.core.PApplet;
import processing.core.PConstants;
import processing.core.PGraphics;
import wblut.geom.WB_Circle;
import wblut.geom.WB_Coord;
import wblut.geom.WB_Line;
import wblut.geom.WB_PolyLine;
import wblut.geom.WB_Polygon;
import wblut.geom.WB_Ray;
import wblut.geom.WB_Ring;
import wblut.geom.WB_Segment;
import wblut.geom.WB_Triangle;
import wblut.geom.WB_Triangulation2D;

/**
 *
 */
public class WB_Render2D extends WB_Processing {
	/**
	 *
	 */

	protected WB_Render2D() {
		super();
	}

	/**
	 *
	 *
	 * @param home
	 */
	public WB_Render2D(final PApplet home) {
		super(home);
		if (home.g == null) {
			throw new IllegalArgumentException("WB_Render3D can only be used after size()");
		}
	}

	/**
	 *
	 *
	 * @param home
	 */
	public WB_Render2D(final PGraphics home) {
		this.home = home;
	}

	/**
	 *
	 *
	 * @param p
	 */

	public void vertex2D(final WB_Coord p) {
		home.vertex(p.xf(), p.yf());

	}

	/**
	 *
	 *
	 * @param p
	 */

	public void drawPoint2D(final WB_Coord p) {
		home.point(p.xf(), p.yf());
	}

	/**
	 *
	 *
	 * @param p
	 * @param r
	 */

	public void drawPoint2D(final WB_Coord p, final double r) {
		home.ellipse(p.xf(), p.yf(), 2 * (float) r, 2 * (float) r);
	}

	/**
	 *
	 *
	 * @param points
	 */

	public void drawPoint2D(final Collection<? extends WB_Coord> points) {
		for (final WB_Coord p : points) {
			drawPoint2D(p);
		}
	}

	/**
	 *
	 *
	 * @param points
	 */

	public void drawPoint2D(final WB_Coord[] points) {
		for (final WB_Coord p : points) {
			drawPoint2D(p);
		}
	}

	/**
	 *
	 *
	 * @param points
	 * @param r
	 */

	public void drawPoint2D(final Collection<? extends WB_Coord> points, final double r) {
		for (final WB_Coord p : points) {
			drawPoint2D(p, r);
		}
	}

	/**
	 *
	 *
	 * @param points
	 * @param r
	 */

	public void drawPoint2D(final WB_Coord[] points, final double r) {
		for (final WB_Coord p : points) {
			drawPoint2D(p, r);
		}
	}

	/**
	 *
	 *
	 * @param p
	 * @param v
	 * @param r
	 */

	public void drawVector2D(final WB_Coord p, final WB_Coord v, final double r) {
		home.pushMatrix();
		home.translate(p.xf(), p.yf());
		home.line(0f, 0f, (float) (r * v.xd()), (float) (r * v.yd()));
		home.popMatrix();
	}

	/**
	 *
	 *
	 * @param L
	 * @param d
	 */

	public void drawLine2D(final WB_Line L, final double d) {
		home.line((float) (L.getOrigin().xd() - d * L.getDirection().xd()),
				(float) (L.getOrigin().yd() - d * L.getDirection().yd()),
				(float) (L.getOrigin().xd() + d * L.getDirection().xd()),
				(float) (L.getOrigin().yd() + d * L.getDirection().yd()));
	}

	/**
	 *
	 *
	 * @param R
	 * @param d
	 */

	public void drawRay2D(final WB_Ray R, final double d) {
		home.line((float) R.getOrigin().xd(), (float) R.getOrigin().yd(),
				(float) (R.getOrigin().xd() + d * R.getDirection().xd()),
				(float) (R.getOrigin().yd() + d * R.getDirection().yd()));
	}

	/**
	 *
	 *
	 * @param p
	 * @param q
	 */

	public void drawSegment2D(final WB_Coord p, final WB_Coord q) {
		home.line((float) p.xd(), (float) p.yd(), (float) q.xd(), (float) q.yd());
	}

	/**
	 *
	 *
	 * @param P
	 */

	public void drawPolyLine2D(final WB_PolyLine P) {
		for (int i = 0; i < P.getNumberOfPoints() - 1; i++) {
			home.line((float) P.getPoint(i).xd(), (float) P.getPoint(i).yd(), (float) P.getPoint(i + 1).xd(),
					(float) P.getPoint(i + 1).yd());
		}
	}

	/**
	 *
	 *
	 * @param P
	 */

	public void drawRing2D(final WB_Ring P) {
		for (int i = 0, j = P.getNumberOfPoints() - 1; i < P.getNumberOfPoints(); j = i++) {
			home.line((float) P.getPoint(j).xd(), (float) P.getPoint(j).yd(), (float) P.getPoint(i).xd(),
					(float) P.getPoint(i).yd());
		}
	}

	/**
	 *
	 *
	 * @param P
	 */

	public void drawPolygon2D(final WB_Polygon P) {
		final int[] tris = P.getTriangles();
		for (int i = 0; i < tris.length; i += 3) {
			drawTriangle2D(P.getPoint(tris[i]), P.getPoint(tris[i + 1]), P.getPoint(tris[i + 2]));
		}
	}

	public void drawPolygon2D(final WB_Polygon[] P) {
		for (WB_Polygon poly : P) {
			drawPolygon2D(poly);
		}

	}

	public void drawPolygon2D(final Collection<? extends WB_Polygon> P) {
		for (WB_Polygon poly : P) {
			drawPolygon2D(poly);
		}

	}

	/**
	 *
	 *
	 * @param P
	 */

	public void drawPolygonEdges2D(final WB_Polygon P) {
		final int[] npc = P.getNumberOfPointsPerContour();
		int index = 0;
		for (int i = 0; i < P.getNumberOfContours(); i++) {
			home.beginShape();
			for (int j = 0; j < npc[i]; j++) {
				vertex2D(P.getPoint(index++));
			}
			home.endShape(PConstants.CLOSE);
		}
	}

	public void drawPolygonEdges2D(final WB_Polygon[] P) {
		for (WB_Polygon poly : P) {
			drawPolygon2D(poly);
		}

	}

	public void drawPolygonEdges2D(final Collection<? extends WB_Polygon> P) {
		for (WB_Polygon poly : P) {
			drawPolygon2D(poly);
		}

	}

	/**
	 *
	 *
	 * @param C
	 */

	public void drawCircle2D(final WB_Circle C) {
		home.ellipse((float) C.getCenter().xd(), (float) C.getCenter().yd(), 2 * (float) C.getRadius(),
				2 * (float) C.getRadius());
	}

	/**
	 *
	 *
	 * @param T
	 */

	public void drawTriangle2D(final WB_Triangle T) {
		home.beginShape(PConstants.TRIANGLES);
		vertex2D(T.p1());
		vertex2D(T.p2());
		vertex2D(T.p3());
		home.endShape();
	}

	/**
	 *
	 *
	 * @param p1
	 * @param p2
	 * @param p3
	 */

	public void drawTriangle2D(final WB_Coord p1, final WB_Coord p2, final WB_Coord p3) {
		home.beginShape(PConstants.TRIANGLES);
		vertex2D(p1);
		vertex2D(p2);
		vertex2D(p3);
		home.endShape();
	}

	/**
	 *
	 *
	 * @param tri
	 * @param points
	 */

	public void drawTriangle2D(final int[] tri, final List<? extends WB_Coord> points) {
		for (int i = 0; i < tri.length; i += 3) {
			home.beginShape(PConstants.TRIANGLES);
			vertex2D(points.get(tri[i]));
			vertex2D(points.get(tri[i + 1]));
			vertex2D(points.get(tri[i + 2]));
			home.endShape();
		}
	}

	/**
	 *
	 *
	 * @param tri
	 * @param points
	 */

	public void drawTriangle2D(final int[] tri, final WB_Coord[] points) {
		for (int i = 0; i < tri.length; i += 3) {
			home.beginShape(PConstants.TRIANGLES);
			vertex2D(points[tri[i]]);
			vertex2D(points[tri[i + 1]]);
			vertex2D(points[tri[i + 2]]);
			home.endShape();
		}
	}

	/**
	 *
	 *
	 * @param tri
	 * @param points
	 */

	public void drawTriangulation2D(final WB_Triangulation2D tri, final List<? extends WB_Coord> points) {
		final int[] triangles = tri.getTriangles();
		home.beginShape(PConstants.TRIANGLES);
		for (int i = 0; i < triangles.length; i += 3) {
			vertex2D(points.get(triangles[i]));
			vertex2D(points.get(triangles[i + 1]));
			vertex2D(points.get(triangles[i + 2]));
		}
		home.endShape();
	}

	/**
	 *
	 *
	 * @param tri
	 * @param points
	 */

	public void drawTriangulationEdges2D(final WB_Triangulation2D tri, final List<? extends WB_Coord> points) {
		final int[] edges = tri.getEdges();
		for (int i = 0; i < edges.length; i += 2) {
			drawSegment2D(points.get(edges[i]), points.get(edges[i + 1]));
		}
	}

	/**
	 *
	 *
	 * @param tri
	 * @param points
	 */

	public void drawTriangulation2D(final WB_Triangulation2D tri, final WB_Coord[] points) {
		final int[] triangles = tri.getTriangles();
		home.beginShape(PConstants.TRIANGLES);
		for (int i = 0; i < triangles.length; i += 3) {
			vertex2D(points[triangles[i]]);
			vertex2D(points[triangles[i + 1]]);
			vertex2D(points[triangles[i + 2]]);
		}
		home.endShape();
	}

	/**
	 *
	 *
	 * @param tri
	 * @param points
	 */

	public void drawTriangulationEdges2D(final WB_Triangulation2D tri, final WB_Coord[] points) {
		final int[] edges = tri.getEdges();
		for (int i = 0; i < edges.length; i += 2) {
			drawSegment2D(points[edges[i]], points[edges[i + 1]]);
		}
	}

	/**
	 *
	 *
	 * @param triangles
	 */

	public void drawTriangle2D(final Collection<? extends WB_Triangle> triangles) {
		final Iterator<? extends WB_Triangle> triItr = triangles.iterator();
		while (triItr.hasNext()) {
			drawTriangle2D(triItr.next());
		}
	}

	/**
	 *
	 *
	 * @param triangles
	 */

	public void drawTriangle2D(final WB_Triangle[] triangles) {
		for (final WB_Triangle triangle : triangles) {
			drawTriangle2D(triangle);
		}
	}

	/**
	 *
	 *
	 * @param triangles
	 */

	public void drawTriangle2DEdges(final Collection<? extends WB_Triangle> triangles) {
		final Iterator<? extends WB_Triangle> triItr = triangles.iterator();
		while (triItr.hasNext()) {
			drawTriangle2DEdges(triItr.next());
		}
	}

	/**
	 *
	 *
	 * @param triangle
	 */

	public void drawTriangle2DEdges(final WB_Triangle triangle) {
		line2D(triangle.p1(), triangle.p2());
		line2D(triangle.p3(), triangle.p2());
		line2D(triangle.p1(), triangle.p3());
	}

	/**
	 *
	 *
	 * @param triangles
	 */

	public void drawTriangle2DEdges(final WB_Triangle[] triangles) {
		for (final WB_Triangle triangle : triangles) {
			drawTriangle2DEdges(triangle);
		}
	}

	/**
	 *
	 *
	 * @param segments
	 */

	public void drawSegment2D(final Collection<? extends WB_Segment> segments) {
		final Iterator<? extends WB_Segment> segItr = segments.iterator();
		while (segItr.hasNext()) {
			drawSegment2D(segItr.next());
		}
	}

	/**
	 *
	 *
	 * @param segment
	 */

	public void drawSegment2D(final WB_Segment segment) {
		line2D(segment.getOrigin(), segment.getEndpoint());
	}

	/**
	 *
	 *
	 * @param segments
	 */

	public void drawSegment2D(final WB_Segment[] segments) {
		for (final WB_Segment segment : segments) {
			drawSegment2D(segment);
		}
	}

	public void drawSegment2D(final int[] segs, final List<? extends WB_Coord> coords) {
		for (int i = 0; i < segs.length; i += 2) {
			line2D(coords.get(segs[i]), coords.get(segs[i + 1]));
		}
	}

	/**
	 *
	 *
	 * @param p
	 */

	public void translate2D(final WB_Coord p) {
		home.translate(p.xf(), p.yf());
	}

	/**
	 *
	 *
	 * @param p
	 * @param q
	 */
	private void line2D(final WB_Coord p, final WB_Coord q) {
		home.beginShape(PConstants.LINES);
		vertex2D(p);
		vertex2D(q);
		home.endShape();
	}
}
