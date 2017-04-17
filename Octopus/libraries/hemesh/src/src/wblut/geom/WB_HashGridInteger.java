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

import gnu.trove.map.TLongIntMap;
import gnu.trove.map.hash.TLongIntHashMap;

public class WB_HashGridInteger {

	private final TLongIntMap values;

	private final int defaultValue;

	private final int K, L, M, KL;

	/**
	 *
	 *
	 * @param K
	 * @param L
	 * @param M
	 * @param defaultValue
	 */
	public WB_HashGridInteger(final int K, final int L, final int M, final int defaultValue) {
		this.K = K;
		this.L = L;
		this.M = M;
		KL = K * L;
		this.defaultValue = defaultValue;
		values = new TLongIntHashMap(10, 0.5f, -1L, defaultValue);
	}

	/**
	 *
	 *
	 * @param K
	 * @param L
	 * @param M
	 */
	public WB_HashGridInteger(final int K, final int L, final int M) {
		this.K = K;
		this.L = L;
		this.M = M;
		KL = K * L;
		defaultValue = Integer.MIN_VALUE;
		values = new TLongIntHashMap(10, 0.5f, -1L, defaultValue);
	}

	/**
	 *
	 *
	 * @param value
	 * @param i
	 * @param j
	 * @param k
	 * @return
	 */
	public boolean setValue(final int value, final int i, final int j, final int k) {
		final long id = safeIndex(i, j, k);
		if (id > 0) {
			values.put(id, value);
			return true;
		}
		return false;
	}

	/**
	 *
	 *
	 * @param value
	 * @param i
	 * @param j
	 * @param k
	 * @return
	 */
	public boolean addValue(final int value, final int i, final int j, final int k) {
		final long id = safeIndex(i, j, k);
		if (id > 0) {
			final int v = values.get(id);
			if (v == defaultValue) {
				values.put(id, value);
			} else {
				values.put(id, v + value);
			}
			return true;
		}
		return false;
	}

	/**
	 *
	 *
	 * @param i
	 * @param j
	 * @param k
	 * @return
	 */
	public boolean clearValue(final int i, final int j, final int k) {
		final long id = safeIndex(i, j, k);
		if (id > 0) {
			values.remove(id);
			return true;
		}
		return false;
	}

	/**
	 *
	 *
	 * @param i
	 * @param j
	 * @param k
	 * @return
	 */
	public double getValue(final int i, final int j, final int k) {
		final long id = safeIndex(i, j, k);
		if (id == -1) {
			return defaultValue;
		}
		if (id > 0) {
			return values.get(id);
		}
		return defaultValue;
	}

	/**
	 *
	 *
	 * @param i
	 * @param j
	 * @param k
	 * @return
	 */
	private long safeIndex(final int i, final int j, final int k) {
		if (i < 0) {
			return -1;
		}
		if (i > K - 1) {
			return -1;
		}
		if (j < 0) {
			return -1;
		}
		if (j > L - 1) {
			return -1;
		}
		if (k < 0) {
			return -1;
		}
		if (k > M - 1) {
			return -1;
		}
		return i + j * K + k * KL;
	}

	/**
	 *
	 *
	 * @return
	 */
	public int getW() {
		return K;
	}

	/**
	 *
	 *
	 * @return
	 */
	public int getH() {
		return L;
	}

	/**
	 *
	 *
	 * @return
	 */
	public int getD() {
		return M;
	}

	/**
	 *
	 *
	 * @return
	 */
	public int getDefaultValue() {
		return defaultValue;
	}

	/**
	 *
	 *
	 * @return
	 */
	public long[] getKeys() {
		return values.keys();
	}

	/**
	 *
	 *
	 * @return
	 */
	public int size() {
		return values.size();
	}
}
