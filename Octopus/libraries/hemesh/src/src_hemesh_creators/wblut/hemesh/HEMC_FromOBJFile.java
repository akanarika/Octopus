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

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.text.MessageFormat;
import java.util.List;
import java.util.zip.GZIPInputStream;

import com.momchil_atanasov.data.front.parser.IOBJParser;
import com.momchil_atanasov.data.front.parser.OBJDataReference;
import com.momchil_atanasov.data.front.parser.OBJFace;
import com.momchil_atanasov.data.front.parser.OBJMesh;
import com.momchil_atanasov.data.front.parser.OBJModel;
import com.momchil_atanasov.data.front.parser.OBJObject;
import com.momchil_atanasov.data.front.parser.OBJParser;
import com.momchil_atanasov.data.front.parser.OBJVertex;

import wblut.geom.WB_Coord;
import wblut.geom.WB_Point;

/**
 *
 */
public class HEMC_FromOBJFile extends HEMC_MultiCreator {

	/**
	 *
	 */
	private String path;

	/**
	 *
	 */
	private double scale;

	/**
	 *
	 */
	public HEMC_FromOBJFile() {
		super();
		scale = 1;
		path = null;

	}

	/**
	 *
	 *
	 * @param path
	 */
	public HEMC_FromOBJFile(final String path) {
		super();
		this.path = path;
		scale = 1;

	}

	/**
	 *
	 *
	 * @param path
	 * @return
	 */
	public HEMC_FromOBJFile setPath(final String path) {
		this.path = path;
		return this;
	}

	/**
	 * 
	 *
	 * @param f 
	 * @return 
	 */
	public HEMC_FromOBJFile setScale(final double f) {
		scale = f;
		return this;
	}

	/* (non-Javadoc)
	 * @see wblut.hemesh.HEMC_MultiCreator#create()
	 */
	@Override
	public HE_MeshCollection create() {
		final HE_MeshCollection result = new HE_MeshCollection();

		if (path == null) {
			return result;
		}
		final File file = new File(path);
		final InputStream is = createInputStream(file);
		if (is == null) {
			return result;
		}
		try {
			final IOBJParser parser = new OBJParser();
			final OBJModel model = parser.parse(is);
			is.close();
			System.out.println(MessageFormat.format(
					"OBJ model has {0} vertices, {1} normals, {2} texture coordinates, and {3} objects.",
					model.getVertices().size(), model.getNormals().size(), model.getTexCoords().size(),
					model.getObjects().size()));

			List<OBJVertex> vertices = model.getVertices();
			//List<OBJTexCoord> texcoords = model.getTexCoords();
			WB_Coord[] newVertices = new WB_Coord[vertices.size()];
			int i = 0;
			for (OBJVertex v : vertices) {
				newVertices[i++] = new WB_Point(scale * v.x, scale * v.y, scale * v.z);

			}

			i = 0;
			for (OBJObject object : model.getObjects()) {
				int nof = 0;
				for (OBJMesh mesh : object.getMeshes()) {
					nof += mesh.getFaces().size();

				}
				int[][] faces = new int[nof][];
				for (OBJMesh mesh : object.getMeshes()) {
					for (OBJFace face : mesh.getFaces()) {
						faces[i] = new int[face.getReferences().size()];
						int j = 0;
						for (OBJDataReference reference : face.getReferences()) {
							faces[i][j] = reference.vertexIndex;
							j++;
						}
						i++;
					}
				}

				final HEC_FromFacelist creator = new HEC_FromFacelist();
				creator.setVertices(newVertices);
				creator.setFaces(faces);
				creator.setDuplicate(true);

				HE_Mesh obj = new HE_Mesh(creator);
				obj.cleanUnusedElementsByFace();
				obj.capHalfedges();
				result.add(obj);

			}
		} catch (final IOException e) {
			e.printStackTrace();
		}
		return result;
	}

	// Code excerpts form processing.core
	/**
	 *
	 *
	 * @param file
	 * @return
	 */
	private InputStream createInputStream(final File file) {
		if (file == null) {
			throw new IllegalArgumentException("file can't be null");
		}
		try {
			InputStream stream = new FileInputStream(file);
			if (file.getName().toLowerCase().endsWith(".gz")) {
				stream = new GZIPInputStream(stream);
			}
			return stream;
		} catch (final IOException e) {
			e.printStackTrace();
		}
		return null;
	}

}
