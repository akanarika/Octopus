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

import gnu.trove.map.TLongObjectMap;
import gnu.trove.map.hash.TLongObjectHashMap;


public class WB_HashGridObject {


	@SuppressWarnings("rawtypes")
	private final TLongObjectMap values;


	private final Object defaultValue;


	private final int K, L, M, KL;


	/**
	 *
	 *
	 * @param K
	 * @param L
	 * @param M
	 * @param defaultValue
	 */
	@SuppressWarnings("rawtypes")
	public WB_HashGridObject(final int K, final int L, final int M, final Object defaultValue) {
		this.K = K;
		this.L = L;
		this.M = M;
		KL = K * L;
		this.defaultValue = defaultValue;
		values = new TLongObjectHashMap(10, 0.5f, -1L);
	}


	/**
	 *
	 *
	 * @param K
	 * @param L
	 * @param M
	 */
	@SuppressWarnings("rawtypes")
	public WB_HashGridObject(final int K, final int L, final int M) {
		this.K = K;
		this.L = L;
		this.M = M;
		KL = K * L;
		defaultValue = null;
		values = new TLongObjectHashMap(10, 0.5f, -1L);
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
	@SuppressWarnings("unchecked")
	public boolean setValue(final Object value, final int i, final int j, final int k) {
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
	public Object getValue(final int i, final int j, final int k) {
		final long id = safeIndex(i, j, k);
		if (id == -1) {
			return defaultValue;
		}
		if (id > 0) {
			final Object val = values.get(id);
			return val;
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
		if (i > (K - 1)) {
			return -1;
		}
		if (j < 0) {
			return -1;
		}
		if (j > (L - 1)) {
			return -1;
		}
		if (k < 0) {
			return -1;
		}
		if (k > (M - 1)) {
			return -1;
		}
		return i + (j * K) + (k * KL);
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
	public Object getDefaultValue() {
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
