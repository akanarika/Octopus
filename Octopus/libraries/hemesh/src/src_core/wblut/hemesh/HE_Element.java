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

import java.util.concurrent.atomic.AtomicLong;

/**
 *
 */
public abstract class HE_Element {

	protected static AtomicLong currentKey = new AtomicLong(0);
	protected final long key;
	protected long labels;

	/**
	 *
	 */
	public HE_Element() {
		key = currentKey.getAndAdd(1);
		labels = mergeLabels(-1, -1);
	}

	private static long mergeLabels(final int internal, final int external) {
		return (long) internal << 32 | external & 0xffffffffL;

	}

	protected final void setInternalLabel(final int label) {
		labels = mergeLabels(label, getLabel());
	}

	/**
	 *
	 *
	 * @param label
	 */
	public final void setLabel(final int label) {
		labels = mergeLabels(getInternalLabel(), label);
	}

	/**
	 *
	 *
	 * @return
	 */
	public final long getKey() {
		return key;
	}

	/**
	 *
	 *
	 * @return
	 */
	public final int getInternalLabel() {
		return (int) (labels >> 32);

	}

	/**
	 *
	 *
	 * @return
	 */
	public final int getLabel() {
		return (int) labels;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode() {
		return (int) (key ^ key >>> 32);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(final Object other) {
		if (other == null) {
			return false;
		}
		if (other == this) {
			return true;
		}
		if (!(other instanceof HE_Element)) {
			return false;
		}
		return ((HE_Element) other).getKey() == key;
	}

	/**
	 *
	 *
	 * @param el
	 */
	public void copyProperties(final HE_Element el) {
		labels = mergeLabels(el.getInternalLabel(), el.getLabel());
	}

	/**
	 *
	 */
	protected abstract void clear();
}
