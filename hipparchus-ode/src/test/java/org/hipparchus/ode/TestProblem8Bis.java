package org.hipparchus.ode;

import static org.junit.Assert.*;

import org.junit.Test;
import org.hipparchus.geometry.euclidean.threed.Rotation;
import org.hipparchus.geometry.euclidean.threed.RotationConvention;
import org.hipparchus.geometry.euclidean.threed.RotationOrder;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.special.elliptic.jacobi.CopolarN;
import org.hipparchus.special.elliptic.jacobi.JacobiElliptic;
import org.hipparchus.special.elliptic.jacobi.JacobiEllipticBuilder;
import org.hipparchus.special.elliptic.legendre.LegendreEllipticIntegral;
import org.hipparchus.util.FastMath;

public class TestProblem8Bis extends TestProblemAbstract {

	/** Moments of inertia. */
	final double i1;
	final double i2;
	final double i3;

	/** Moments of inertia converted. */
	final double i1C;
	final double i2C;
	final double i3C;

	/** Substraction of inertias. */
	final double i32;
	final double i31;
	final double i21;

	/** Initial state. */
	final double[] y0;

	/** Converted state. */
	final double[] y0C;

	/** Twice the angular kinetic energy. */
	final double twoE;

	/** Square of kinetic momentum. */
	final double m2;

	/** State scaling factor. */
	final double o1Scale;

	/** State scaling factor. */
	final double o2Scale;

	/** State scaling factor. */
	final double o3Scale;

	/** Jacobi elliptic function. */
	final JacobiElliptic jacobi;

	/** Elliptic modulus k2 (k2 = m). */
	final double k2;

	/** Time scaling factor. */
	final double tScale;

	/** Time reference for rotation rate. */
	final double tRef;

	/**Offset rotation  between initial inertial frame and the frame with moment vector and Z axis aligned. */
	Rotation mAlignedToInert;

	/** Initial converted rotation. */
	final Rotation r0;

	/** Rotation to switch to the converted axes frame. */
	final Rotation convertAxes;

	public TestProblem8Bis() {
		//Arguments in the super constructor :
		//Intital time, Primary state (o1, o2, o3, q0, q1, q2, q3), Final time, Error scale
		super(0.0, new double[] {5.0, 0.0, 4.0, 0.9, 0.437, 0.0, 0.0}, 20.0, new double[] { 1.0, 1.0, 1.0 });
		i2 = 3.0 / 8.0;
		i1 = 1.0 / 2.0;
		i3 = 5.0 / 8.0;


		y0 = getInitialState().getPrimaryState();
		final double t0 = getInitialState().getTime();

		final double o12 = y0[0] * y0[0];
		final double o22 = y0[1] * y0[1];
		final double o32 = y0[2] * y0[2];

		twoE    =  i1 * o12 + i2 * o22 + i3 * o32;
		m2      =  i1 * i1 * o12 + i2 * i2 * o22 + i3 * i3 * o32;

		r0 = new Rotation(y0[3], y0[4], y0[5], y0[6], true);//InertToBody0, without any conversion, only with initial quaternion

		final Vector3D m0 = new Vector3D(i1 * y0[0], i2 * y0[1], i3 * y0[2]);
		final double theta0 = FastMath.acos(m0.getZ()/m0.getNorm());
		final double psi0 = FastMath.atan2(m0.getX(), m0.getY());
		final double phi0 = 0;

		final Rotation mAlignedToBody = new Rotation(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM, phi0, theta0, psi0);

		mAlignedToInert = r0.applyInverseTo(mAlignedToBody);

		
		//compute axes conversion
		final double[][] converted = sortInertiaAxis();
		y0C = converted[0];
		i1C = converted[1][0];
		i2C = converted[1][1];
		i3C = converted[1][2];
		
		//conversion axis offset
		final Vector3D X = new Vector3D(converted[2][0], converted[2][1], converted[2][2]);
		final Vector3D Y = new Vector3D(converted[3][0], converted[3][1], converted[3][2]);

		convertAxes = new Rotation(Vector3D.PLUS_I, Vector3D.PLUS_J, X, Y);
		
		//compute global calculus elements

		i32  = i3C - i2C;
		i31  = i3C - i1C;
		i21  = i2C - i1C;

		tScale  = FastMath.sqrt(i32 * (m2 - twoE * i1C) / (i1C * i2C * i3C));
		o1Scale = FastMath.sqrt((twoE * i3C - m2) / (i1C * i31));
		o2Scale = FastMath.sqrt((twoE * i3C - m2) / (i2C * i32));
		o3Scale = FastMath.sqrt((m2 - twoE * i1C) / (i3C * i31));
		
		k2     = i21 * (twoE * i3C - m2) / (i32 * (m2 - twoE * i1C));

		jacobi = JacobiEllipticBuilder.build(k2);

		tRef   = t0 - jacobi.arcsn(y0[1] / o2Scale) / tScale;
		System.out.println("Tref : "+tRef);
		
		
	}


	public double[][] sortInertiaAxis() {

		Vector3D[] axesP = {Vector3D.PLUS_I, Vector3D.PLUS_J, Vector3D.PLUS_K};

		double[] i = {i1, i2, i3};
		double[] y0P = y0;

		System.out.println("omega avant : "+y0P[0]+" "+ y0P[1]+" "+ y0P[2]);
		System.out.println("Initial : iA1 = "+i1+" "+axesP[0]);
		System.out.println("iA2 = "+i2 + " "+axesP[1]);
		System.out.println("iA3 = "+i3 + " "+axesP[2]);
		if (i[0] > i[1]) {
			Vector3D z = axesP[0];
			axesP[0] = axesP[1];
			axesP[1] = z;
			axesP[2] = axesP[2].negate();

			double y = i[0];
			i[0] = i[1];
			i[1] = y;

			double v = y0P[0];
			y0P[0] = y0P[1];
			y0P[1] = v;
			y0P[2] = -y0P[2];
			System.out.println("1ere boucle : iA1 = "+i[0]+" "+axesP[0]);
			System.out.println("iA2 = "+i[1] + " "+axesP[1]);
			System.out.println("iA3 = "+i[2] + " "+axesP[2]);
		}

		if (i[1] > i[2]) {
			Vector3D z = axesP[1];
			axesP[1] = axesP[2];
			axesP[2] = z;
			axesP[0] = axesP[0].negate();

			double y = i[1];
			i[1] = i[2];
			i[2] = y;

			double v = y0P[1];
			y0P[1] = y0P[2];
			y0P[2] = v;
			y0P[0] = -y0P[0];
			System.out.println("2ere boucle : iA1 = "+i[0]+" "+axesP[0]);
			System.out.println("iA2 = "+i[1] + " "+axesP[1]);
			System.out.println("iA3 = "+i[2] + " "+axesP[2]);
		}

		if (i[0] > i[1]) {
			Vector3D z = axesP[0];
			axesP[0] = axesP[1];
			axesP[1] = z;
			axesP[2] = axesP[2].negate();

			double y = i[0];
			i[0] = i[1];
			i[1] = y;

			double v = y0P[0];
			y0P[0] = y0P[1];
			y0P[1] = v;
			y0P[2] = -y0P[2];
			System.out.println("3ere boucle : iA1 = "+i[0]+" "+axesP[0]);
			System.out.println("iA2 = "+i[1] + " "+axesP[1]);
			System.out.println("iA3 = "+i[2] + " "+axesP[2]);
		}

		final double condition;
		if( y0P[0] == 0 && y0P[1] == 0 && y0P[2] == 0) {
			condition = 0.0;
		}else {
			condition = m2/twoE;
		}

		System.out.println("Condition : "+condition);
		if(condition < i[1]) {
			Vector3D z = axesP[0];
			axesP[0] = axesP[2];
			axesP[2] = z;
			axesP[1] = axesP[1].negate();

			double y = i[0];
			i[0] = i[2];
			i[2] = y;

			double v = y0P[0];
			y0P[0] = y0P[2];
			y0P[2] = v;
			y0P[1] = - y0P[1];
			System.out.println("repère final  : iA1 = "+i[0]+" "+axesP[0]);
			System.out.println("iA2 = "+i[1] + " "+axesP[1]);
			System.out.println("iA3 = "+i[2] + " "+axesP[2]);
		}
		System.out.println("omega après : "+y0P[0]+" "+ y0P[1]+" "+ y0P[2]);
		System.out.println("axes après : "+axesP[0]+" "+ axesP[1]+" "+ axesP[2]);
		System.out.println("inerties après : "+i[0]+" "+ i[1]+" "+ i[2]);

		return new double[][] { y0P, i, {axesP[0].getX(), axesP[0].getY(), axesP[0].getZ()}, {axesP[1].getX(), axesP[1].getY(), axesP[1].getZ()} };
	}

	
	public double[][] computeTorqueFreeMotion(double i1, double i2, double i3, double t0, double[] y0, double t) {

		/** Coefficients for φ model. */
		final double phiSlope;
		final double b;
		final double c;
		final double d;

		/**Period with respect to time. */
		final double period;

		/**DenseOutputModel of phi. */
		final DenseOutputModel phiQuadratureModel;
		final double integOnePeriod;

		// coefficients for φ model
		period = 4 * LegendreEllipticIntegral.bigK(k2) / tScale;
		phiSlope = FastMath.sqrt(m2) / i3;// =a, pente de la partie linéaire
		b = phiSlope * i32 * i31;
		c = i1 * i32;
		d = i3 * i21;

		//Integrate the quadrature phi term on one period
		final DormandPrince853Integrator integ = new DormandPrince853Integrator(1.0e-6 * period, 1.0e-2 * period,
				phiSlope * period * 1.0e-13, 1.0e-13);
		phiQuadratureModel = new DenseOutputModel();
		integ.addStepHandler(phiQuadratureModel);

		integ.integrate(new OrdinaryDifferentialEquation() {

			/** {@inheritDoc} */
			@Override
			public int getDimension() {
				return 1;
			}
			@Override
			public double[] computeDerivatives(final double t, final double[] y) {
				final double sn = jacobi.valuesN((t - tRef) * tScale).sn();
				return new double[] {
						b / (c + d * sn * sn)
				};
			}
		},new ODEState(t0, new double[1]), t0 + period);

		integOnePeriod = phiQuadratureModel.getInterpolatedState(phiQuadratureModel.getFinalTime()).getPrimaryState()[0];
		
		//Computation of omega
		final CopolarN valuesN = jacobi.valuesN((t - tRef) * tScale);

		final Vector3D omegaP = new Vector3D(
				o1Scale * valuesN.cn(),
                o2Scale * valuesN.sn(),
				o3Scale * valuesN.dn());

		//We need to convert again omega because the formula used isn't converted with the transformations before
		final Vector3D omega = convertAxes.applyTo(omegaP);

		System.out.println("omega fonction : "+o1Scale * valuesN.cn()+" "+o2Scale * valuesN.sn()+" "+o3Scale * valuesN.dn());
		System.out.println("omega fonction conversion : "+omega.getX()+" "+omega.getY()+" "+omega.getZ());
		
		//Computation of the euler angles
				//Compute rotation rate
				final double   o1            = omega.getX();
				final double   o2            = omega.getY();
				final double   o3            = omega.getZ();

				//Compute angles
				final double   psi           = FastMath.atan2(i1 * o1, i2 * o2);
				final double   theta         = FastMath.acos(o3 / phiSlope);
				final double   phiLinear     = phiSlope * t;

				System.out.println("PSI = "+ psi+"\n"+"THETA = "+theta+"\n"+"PHI LIN = "+phiLinear);
				//Integration for the computation of phi
				final int nbPeriods = (int) FastMath.floor((t - t0) / period);//floor = entier inférieur = nb période entière
				final double tStartInteg = t0 + nbPeriods * period;//partie de période à la fin entre tau Integ et tau end
				final double integPartial = phiQuadratureModel.getInterpolatedState(t - tStartInteg).getPrimaryState()[0];// a vérifier, partie de l'intégrale apres le nb entier de période
				final double phiQuadrature = nbPeriods * integOnePeriod + integPartial;

				final double phi = phiLinear + phiQuadrature;
				
				//Computation of the quaternion

				// Rotation between computation frame (aligned with momentum) and body
				//(It is simply the angles equations provided by L&L)
				final Rotation alignedToBody = new Rotation(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM,
						phi, theta, psi);
				
				final Rotation bodyToInert = mAlignedToInert.applyTo(convertAxes.applyInverseTo(alignedToBody));
				
				final Rotation inertToBody = (mAlignedToInert.applyInverseTo(convertAxes)).applyTo(alignedToBody);
				final double[] angles = inertToBody.getAngles(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM);
				
				return new double[][] {{angles[0], angles[1], angles[2]}};
	}

	public double[] doComputeDerivatives(double t, double[] y) {

		final  double[] yDot = new double[getDimension()];

		// compute the derivatives using Euler equations
		yDot[0] = y[1] * y[2] * (i2 - i3) / i1;
		yDot[1] = y[2] * y[0] * (i3 - i1) / i2;
		yDot[2] = y[0] * y[1] * (i1 - i2) / i3;

		// compute the derivatives using Qpoint = 0.5 * Omega_inertialframe * Q
		yDot[3] = 0.5 * (-y[0] * y[4] -y[1] * y[5] -y[2] * y[6]);
		yDot[4] = 0.5 * (y[0] * y[3] +y[2] * y[5] -y[1] * y[6]);
		yDot[5] = 0.5 * (y[1] * y[3] -y[2] * y[4] +y[0] * y[6]);
		yDot[6] = 0.5 * (y[2] * y[3] +y[1] * y[4] -y[0] * y[5]);
		return yDot;

	}

	public double[] computeTheoreticalState(double t) {

		double t0 = getInitialState().getTime();

				final double[][] tfm = computeTorqueFreeMotion(i1C, i2C, i3C, t0, y0C, t);
				
				final double[] angles = tfm[0];
				return new double[] {
						angles[0],
						angles[1],
						angles[2]
				};
//				final double[] omega = tfm[0];
//				final double[] quaternion = tfm[2];
//				final double[] angles = tfm[1];
//				final double[] X = tfm[3];
//				final double[] Y = tfm[4];
//		
//				return new double[] {
//						omega[0],
//						omega[1],
//						omega[2],
//						quaternion[0],
//						quaternion[1],
//						quaternion[2],
//						quaternion[3],
//						angles[0],
//						angles[1],
//						angles[2],
//						X[0],
//						X[1],
//						X[2],
//						Y[0],
//						Y[1],
//						Y[2]
//				};
	}



}
