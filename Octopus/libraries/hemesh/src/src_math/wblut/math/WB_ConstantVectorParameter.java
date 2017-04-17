/*
 * This file is part of HE_Mesh, a library for creating and manipulating meshes.
 * It is dedicated to the public domain. To the extent possible under law,
 * I , Frederik Vanhoutte, have waived all copyright and related or neighboring
 * rights.
 * 
 * This work is published from Belgium. (http://creativecommons.org/publicdomain/zero/1.0/)
 * 
 */
package wblut.math;

import wblut.geom.WB_Coord;
import wblut.geom.WB_SimpleVector;

/**
 * The Class WB_ConstantParameter.
 *
 * @author Frederik Vanhoutte, W:Blut
 *
 *         A parameter which is constant, i.e. a single unchanging value.
 */
public class WB_ConstantVectorParameter implements WB_VectorParameter {
    /** The value. */
    WB_Coord value;

    /**
     * 
     *
     * @param value 
     */
    public WB_ConstantVectorParameter(final WB_Coord value) {
	this.value = new WB_SimpleVector(value);
    }

    /* (non-Javadoc)
     * @see wblut.math.WB_VectorParameter#evaluate(double[])
     */
    @Override
    public WB_Coord evaluate(final double... x) {
	return value;
    }
}
