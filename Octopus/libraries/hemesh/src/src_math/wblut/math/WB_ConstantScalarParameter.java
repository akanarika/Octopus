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

/**
 * The Class WB_ConstantParameter.
 *
 * @author Frederik Vanhoutte, W:Blut
 *
 *         A parameter which is constant, i.e. a single unchanging value.
 */
public class WB_ConstantScalarParameter implements WB_ScalarParameter {
    /** The value. */
    double value;

    /**
     * 
     *
     * @param value 
     */
    public WB_ConstantScalarParameter(final double value) {
	this.value = value;
    }

    /* (non-Javadoc)
     * @see wblut.math.WB_ScalarParameter#evaluate(double[])
     */
    @Override
    public double evaluate(final double... x) {
	return value;
    }
}
