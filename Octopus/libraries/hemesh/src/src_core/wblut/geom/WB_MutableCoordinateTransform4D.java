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

/**
 * Interface for implementing mutable transformation operations on 4D
 * coordinates.
 *
 * All of the operators defined in the interface change the calling object. All
 * operators use the label "Self", such as "rotateXYSelf" to indicate this.
 *
 * @author Frederik Vanhoutte
 *
 */
public interface WB_MutableCoordinateTransform4D {
	public WB_Vector4D rotateXWSelf(double angle);

	public WB_Vector4D rotateXYSelf(double angle);

	public WB_Vector4D rotateXZSelf(double angle);

	public WB_Vector4D rotateYWSelf(double angle);

	public WB_Vector4D rotateYZSelf(double angle);

	public WB_Vector4D rotateZWSelf(double angle);
}
