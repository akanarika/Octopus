/**
 *
 */
package wblut.math;

/**
 * @author FVH
 *
 */
public class WB_SphericalHarmonics {

	public static WB_Complex Y(final int l, final int m, final double theta, final double phi) {
		int mPos = Math.abs(m);
		double factor = getNormFactor(l, mPos);
		double thetaPart = legendre(l, mPos, Math.cos(theta));
		WB_Complex result = new WB_Complex(Math.cos(mPos * phi), Math.sin(mPos * phi));

		result.mulSelf(factor * thetaPart);
		if (m < 0) {
			result.conjugateSelf();
			if (mPos % 2 == 1) {
				result.negateSelf();
			}
		}
		return result;
	}

	public static double Ylm(final int l, final int m, final double theta, final double phi) {
		if (m == 0) {
			return Y(l, m, theta, phi).re();
		} else if (m < 0) {
			return Math.sqrt(2.0) * Math.pow(-1, m) * Y(l, Math.abs(m), theta, phi).im();
		} else {
			return Math.sqrt(2.0) * Math.pow(-1, m) * Y(l, m, theta, phi).re();

		}
	}

	private static double getNormFactor(final int l, final int m) {
		double factor = 1.0;
		for (int i = l - m + 1; i <= l + m; i++) {
			factor *= i;
		}
		factor = (2 * l + 1) / factor;
		factor /= 4.0 * Math.PI;
		factor = Math.sqrt(factor);
		return factor;
	}

	public static double legendre(final int l, final int m, final double x) {
		double fact = 0.0;
		double pll = 0.0;
		double pmm = 0.0;
		double pmmp1 = 0.0;
		double somx2 = 0.0;
		int i = 0;
		int ll = 0;

		if (m < 0 || m > l || Math.abs(x) > 1.0) {
			return Double.NaN;
		}
		pmm = 1.0;
		if (m > 0) {
			somx2 = Math.sqrt((1.0 - x) * (1.0 + x));
			fact = 1.0;
			for (i = 1; i <= m; i++) {
				pmm *= -fact * somx2;
				fact += 2.0;
			}
		}
		if (l == m) {
			return pmm;
		} else {
			pmmp1 = x * (2 * m + 1) * pmm;
			if (l == m + 1) {
				return pmmp1;
			} else {
				for (ll = m + 2; ll <= l; ll++) {
					pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
					pmm = pmmp1;
					pmmp1 = pll;
				}
				return pll;
			}
		}
	}

}
