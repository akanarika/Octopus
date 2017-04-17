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

public class WB_ProgressReporter extends Thread {

	WB_ProgressTracker tracker;
	String status;

	/**
	 *
	 */
	public WB_ProgressReporter(final int depth) {
		super();
		tracker = WB_ProgressTracker.instance();
		tracker.setDepth(depth);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Thread#start()
	 */
	@Override
	public void start() {
		super.start();
		System.out.println("Starting WB_ProgressTracker");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Thread#run()
	 */
	@Override
	public void run() {
		while (!Thread.interrupted()) {
			try {

				while (tracker.isUpdated()) {
					status = tracker.getStatus();
					if (status != null) {
						System.out.println(status);
					}
				}
			} catch (final Exception e) {

				e.printStackTrace();
			}
		}
	}
}
