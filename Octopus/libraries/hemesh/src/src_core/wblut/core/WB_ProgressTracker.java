/*
 * This file is part of HE_Mesh, a library for creating and manipulating meshes.
 * It is dedicated to the public domain. To the extent possible under law,
 * I , Frederik Vanhoutte, have waived all copyright and related or neighboring
 * rights.
 *
 * This work is published from Belgium. (http://creativecommons.org/publicdomain/zero/1.0/)
 *
 */

package wblut.core;

import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

public class WB_ProgressTracker {

	protected Queue<Status> statuses;
	protected volatile int depth;
	private static int indent = 3;
	protected volatile int maxdepth;

	/**
	 *
	 */
	protected WB_ProgressTracker() {
		statuses = new ConcurrentLinkedQueue<Status>();
		depth = 0;
		maxdepth = 2;
	}

	private static final WB_ProgressTracker tracker = new WB_ProgressTracker();

	/**
	 *
	 *
	 * @return
	 */
	public static WB_ProgressTracker instance() {
		return tracker;
	}

	/**
	 *
	 *
	 * @param indent
	 */
	public void setIndent(final int indent) {
		WB_ProgressTracker.indent = Math.max(0, indent);

	}

	/**
	 *
	 *
	 * @return
	 */
	protected String getStatus() {
		if (statuses.size() > 0) {
			return statuses.poll().getStatus();
		}
		return "";

	}

	/**
	 *
	 *
	 * @param caller
	 * @param status
	 * @param inc
	 */
	public void setStatus(final Object caller, final String status, final int inc) {
		if (inc < 0) {
			depth = Math.max(0, depth + inc);
		}
		if (depth <= maxdepth) {
			statuses.add(new Status(caller.getClass().getSimpleName(), status, depth));
		}
		if (inc > 0) {
			depth = Math.max(0, depth + inc);
		}
	}

	/**
	 *
	 *
	 * @param caller
	 * @param status
	 * @param inc
	 */
	public void setStatusByString(final String caller, final String status, final int inc) {
		if (inc < 0) {
			depth = Math.max(0, depth + inc);
		}
		if (depth <= maxdepth) {
			statuses.add(new Status(caller, status, depth));
		}
		if (inc > 0) {
			depth = Math.max(0, depth + inc);
		}
	}

	/**
	 *
	 *
	 * @param caller
	 * @param status
	 * @param counter
	 */
	public void setStatus(final Object caller, final String status, final WB_ProgressCounter counter) {

		if (counter.getLimit() > 0) {
			counter.caller = caller.getClass().getSimpleName();
			counter.text = status;
			if (depth <= maxdepth) {
				statuses.add(new Status(caller.getClass().getSimpleName(), status, counter, depth));
			}
		}
	}

	/**
	 *
	 *
	 * @param caller
	 * @param status
	 * @param counter
	 */
	public void setStatusByString(final String caller, final String status, final WB_ProgressCounter counter) {
		if (counter.getLimit() > 0) {
			if (depth <= maxdepth) {
				statuses.add(new Status(caller, status, counter, depth));
			}
		}
	}

	/**
	 *
	 *
	 * @return
	 */
	public boolean isUpdated() {

		return statuses.size() > 0;

	}

	class Status {
		String caller;
		String text;
		String counter;
		String depth;
		int level;

		/**
		 *
		 *
		 * @param caller
		 * @param text
		 * @param counter
		 * @param depth
		 */
		Status(final String caller, final String text, final WB_ProgressCounter counter, final int depth) {
			this.caller = caller;
			this.text = text;
			StringBuffer outputBuffer = new StringBuffer(depth);
			for (int i = 0; i < depth * indent; i++) {
				outputBuffer.append(" ");
			}
			this.depth = outputBuffer.toString();

			this.counter = counter.getLimit() > 0 ? " (" + counter.getCount() + " of " + counter.getLimit() + ")" : "";
			level = depth;
		}

		/**
		 *
		 *
		 * @param caller
		 * @param text
		 * @param depth
		 */
		Status(final String caller, final String text, final int depth) {
			this.caller = caller;
			this.text = text;
			StringBuffer outputBuffer = new StringBuffer(depth);
			for (int i = 0; i < depth * indent; i++) {
				outputBuffer.append(" ");
			}
			this.depth = outputBuffer.toString();
			this.counter = null;
		}

		/**
		 *
		 *
		 * @return
		 */
		String getStatus() {
			if (caller == null) {
				return null;
			}
			if (text == " ") {
				return "";
			}
			if (counter == null) {
				return depth + caller + ": " + text;
			}
			return depth + caller + ": " + text + counter;

		}

	}

	public void setDepth(final int depth) {
		maxdepth = depth;

	}
}
