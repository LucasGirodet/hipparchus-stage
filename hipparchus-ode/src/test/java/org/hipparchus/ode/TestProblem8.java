package org.hipparchus.ode;

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
import org.junit.Assert;

public class TestProblem8 extends TestProblemAbstract {

    /** Moments of inertia. */
    final double i1;
    final double i2;
    final double i3;

    /** Moments of inertia converted. */
    final double i1C;
    final double i2C;
    final double i3C;

    /**Substraction of inertias. */
    final double i32;
    final double i31;
    final double i21;

    /** Initial state. */
    final double[] y0;

    /** Initial time. */
    final double t0;

    /** Converted initial state. */
    final double[] y0C;

    /** Converted axes. */
    final Vector3D[] axes;

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
    final Rotation r0Conv;

    /** Rotation to switch to the converted axes frame. */
    final Rotation convertAxes;

    /** Period of rotation rate. */
    final double period;

    /** Slope of the linear part of the phi model. */
    final double phiSlope;

    /** DenseOutputModel of phi. */
    final DenseOutputModel phiQuadratureModel;

    /** Integral part of quadrature model over one period. */
    final double integOnePeriod;

    /**
     * Simple constructor.
     * @param t0 initial time
     * @param t1 final time
     * @param omega0 initial rotation rate
     * @param r0 initial rotation
     * @param i1 inertia along first axis
     * @param i2 inertia along second axis
     * @param i3 inertia along third axis
     */
    public TestProblem8(final double t0, final double t1, final Vector3D omega0, final Rotation r0,
            final double i1, final double i2, final double i3) {
        //Arguments in the super constructor :
        //Initial time, Primary state (o1, o2, o3, q0, q1, q2, q3), Final time, Error scale
        super(t0, new double[] {
                omega0.getX(), omega0.getY(), omega0.getZ(),
                r0.getQ0(), r0.getQ1(), r0.getQ2(), r0.getQ3()
        },
                t1,
                new double[] { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 });
        this.i1 = i1;
        this.i2 = i2;
        this.i3 = i3;

        this.t0 = t0;

        y0 = getInitialState().getPrimaryState();

        final double o12 = y0[0] * y0[0];
        final double o22 = y0[1] * y0[1];
        final double o32 = y0[2] * y0[2];

        twoE    =  i1 * o12 + i2 * o22 + i3 * o32;
        m2      =  i1 * i1 * o12 + i2 * i2 * o22 + i3 * i3 * o32;

        Vector3D[] axesP = {Vector3D.PLUS_I, Vector3D.PLUS_J, Vector3D.PLUS_K};

        double[] iP = {i1, i2, i3};
        double[] y0P = y0.clone();

        if (iP[0] > iP[1]) {
            Vector3D z = axesP[0];
            axesP[0] = axesP[1];
            axesP[1] = z;
            axesP[2] = axesP[2].negate();

            double y = iP[0];
            iP[0] = iP[1];
            iP[1] = y;

            double v = y0P[0];
            y0P[0] = y0P[1];
            y0P[1] = v;
            y0P[2] = -y0P[2];
        }

        if (iP[1] > iP[2]) {
            Vector3D z = axesP[1];
            axesP[1] = axesP[2];
            axesP[2] = z;
            axesP[0] = axesP[0].negate();

            double y = iP[1];
            iP[1] = iP[2];
            iP[2] = y;

            double v = y0P[1];
            y0P[1] = y0P[2];
            y0P[2] = v;
            y0P[0] = -y0P[0];
        }

        if (iP[0] > iP[1]) {
            Vector3D z = axesP[0];
            axesP[0] = axesP[1];
            axesP[1] = z;
            axesP[2] = axesP[2].negate();

            double y = iP[0];
            iP[0] = iP[1];
            iP[1] = y;

            double v = y0P[0];
            y0P[0] = y0P[1];
            y0P[1] = v;
            y0P[2] = -y0P[2];
        }


        final double condition = m2/twoE;

        if (condition < iP[1]) {
            Vector3D z = axesP[0];
            axesP[0] = axesP[2];
            axesP[2] = z;
            axesP[1] = axesP[1].negate();

            double y = iP[0];
            iP[0] = iP[2];
            iP[2] = y;

            double v = y0P[0];
            y0P[0] = y0P[2];
            y0P[2] = v;
            y0P[1] = - y0P[1];
        }

        i1C = iP[0];
        i2C = iP[1];
        i3C = iP[2];

        axes = axesP;
        convertAxes = new Rotation(Vector3D.PLUS_I, Vector3D.PLUS_J, axes[0], axes[1]);

        y0C = y0P.clone();

        i32  = i3C - i2C;
        i31  = i3C - i1C;
        i21  = i2C - i1C;

        tScale  = FastMath.sqrt(i32 * (m2 - twoE * i1C) / (i1C * i2C * i3C));
        o1Scale = FastMath.sqrt((twoE * i3C - m2) / (i1C * i31));
        o2Scale = FastMath.sqrt((twoE * i3C - m2) / (i2C * i32));
        o3Scale = FastMath.sqrt((m2 - twoE * i1C) / (i3C * i31));

        k2 = i21 * (twoE * i3C - m2) / (i32 * (m2 - twoE * i1C));

        jacobi = JacobiEllipticBuilder.build(k2);

        // convert initial conditions to Euler angles such the M is aligned with Z in computation frame
        final Vector3D omega0Body = new Vector3D(y0C[0], y0C[1], y0C[2]);

        r0Conv         = new Rotation(y0C[3], y0C[4], y0C[5], y0C[6], true); //Initial quaternion
        final Vector3D m0Body = new Vector3D(i1C * omega0Body.getX(), i2C * omega0Body.getY(), i3C * omega0Body.getZ());

        final double   theta0 =  FastMath.acos(m0Body.getZ() / m0Body.getNorm());
        final double   psi0       = FastMath.atan2(m0Body.getX(), m0Body.getY()); // it is really atan2(x, y), not atan2(y, x) as usual!
        final double   phi0       = 0.0; // this angle can be set arbitrarily, so 0 is a fair value (Eq. 37.13 - 37.14)

        //Compute offset rotation between inertial frame aligned with momentum and regular inertial frame
        final Rotation mAlignedToBody = new Rotation(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM,
                phi0, theta0, psi0);

        Rotation r0ConvertedAxis = convertAxes.applyInverseTo(r0Conv);

        mAlignedToInert = r0ConvertedAxis.applyInverseTo(mAlignedToBody);

        period             = 4 * LegendreEllipticIntegral.bigK(k2) / tScale;
        phiSlope           = FastMath.sqrt(m2) / i3C;
        phiQuadratureModel = computePhiQuadratureModel(t0);
        integOnePeriod     = phiQuadratureModel.getInterpolatedState(phiQuadratureModel.getFinalTime()).getPrimaryState()[0];

        tRef = getTRef(t0);
    }

    /** Provide a DenseOutputModel of phi quadrature to allow interpolation.
     *
     * @param t0
     * @return phiQuadratureModel
     */
    private DenseOutputModel computePhiQuadratureModel(final double t0) {

        final double i32 = i3C - i2C;
        final double i31 = i3C - i1C;
        final double i21 = i2C - i1C;

        // coefficients for φ model
        final double b = phiSlope * i32 * i31;
        final double c = i1C * i32;
        final double d = i3C * i21;

        //Integrate the quadrature phi term on one period
        final DormandPrince853Integrator integ = new DormandPrince853Integrator(1.0e-6 * period, 1.0e-2 * period,
                phiSlope * period * 1.0e-13, 1.0e-13);
        final DenseOutputModel model = new DenseOutputModel();
        integ.addStepHandler(model);

        integ.integrate(new OrdinaryDifferentialEquation() {

            /** {@inheritDoc} */
            @Override
            public int getDimension() {
                return 1;
            }

            /** {@inheritDoc} */
            @Override
            public double[] computeDerivatives(final double t, final double[] y) {
                final double sn = jacobi.valuesN((t - tRef) * tScale).sn();
                return new double[] {
                        b / (c + d * sn * sn)
                };
            }

        }, new ODEState(t0, new double[1]), t0 + period);

        return model;

    }
    /** Compute the tfmState elements.
     *
     * @param t
     * @return tfmState
     */
    public TfmState computeTorqueFreeMotion(double t) {

        final double[] omegaValues = omega(t);
        final Vector3D omegaModified = new Vector3D(omegaValues[0], omegaValues[1], omegaValues[2]);//Omega after sortInertiaAxis and modified equations
        final Vector3D omega    = convertAxes.applyInverseTo(omegaModified);//Omega with modified equations

        // compute angles
        final double psi         = FastMath.atan2(i1C * omegaModified.getX(), i2C * omegaModified.getY());
        final double   theta         = FastMath.acos(omegaModified.getZ() / phiSlope);
        final double   phiLinear     = phiSlope * t;

        // Integration for the computation of phi
        final double t0 = getInitialTime();
        final int nbPeriods = (int) FastMath.floor((t - t0) / period);//floor = previous integer = complete period number
        final double tStartInteg = t0 + nbPeriods * period;//part of the period between tau Integ and tau end
        final double integPartial = phiQuadratureModel.getInterpolatedState(t - tStartInteg).getPrimaryState()[0];// a vérifier, partie de l'intégrale apres le nb entier de période
        final double phiQuadrature = nbPeriods * integOnePeriod + integPartial;

        final double phi = phiLinear + phiQuadrature;

        // Computation of the quaternion

        // Rotation between computation frame (aligned with momentum) and body
        //(It is simply the angles equations provided by L&L)
        final Rotation alignedToBody = new Rotation(RotationOrder.ZXZ, RotationConvention.FRAME_TRANSFORM,
                phi, theta, psi);

        // combine with offset rotation to get back to regular inertial frame
        // Inert -> aligned + aligned -> body = inert -> body (What the user wants)
        Rotation inertToBody = alignedToBody.applyTo(mAlignedToInert.revert());

        Rotation originalFrameToBody = convertAxes.applyInverseTo(inertToBody);

        return new TfmState(t, omega, originalFrameToBody, phi, theta, psi, convertAxes, mAlignedToInert);

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
        final TfmState tfm = computeTorqueFreeMotion(t);
        return new double[] {
                tfm.getOmega().getX(),
                tfm.getOmega().getY(),
                tfm.getOmega().getZ(),
                tfm.getRotation().getQ0(),
                tfm.getRotation().getQ1(),
                tfm.getRotation().getQ2(),
                tfm.getRotation().getQ3()
        };
    }

    /** Simple class that contain elements useful for the torque free motion problem.
     *
     */
    public static class TfmState {
        private final double   t;
        private final Vector3D omega;
        private final Rotation rotation;
        private final double phi;
        private final double theta;
        private final double psi;
        private final Rotation convertAxes;
        private final Rotation mAlignedToInert;
        private TfmState(final double t, final Vector3D omega, final Rotation rotation,
                final double phi, final double theta, final double psi,
                final Rotation convertAxes, final Rotation mAlignedToInert) {
            this.t               = t;
            this.omega           = omega;
            this.rotation        = rotation;
            this.phi             = phi;
            this.theta           = theta;
            this.psi             = psi;
            this.convertAxes     = convertAxes;
            this.mAlignedToInert = mAlignedToInert;
        }
        public double getT() {
            return t;
        }
        public Vector3D getOmega() {
            return omega;
        }
        public Rotation getRotation() {
            return rotation;
        }
        public double getPhi() {
            return phi;
        }
        public double getTheta() {
            return theta;
        }
        public double getPsi() {
            return psi;
        }
        public Rotation getConvertAxes() {
            return convertAxes;
        }
        public Rotation getMAlignedToInert() {
            return mAlignedToInert;
        }
    }

    /** Compute a frame that apply correction to the instantaneous rotation vector
     *
     * @param t
     * @return omega the instantaneous rotation vector
     */
    public double[] omega(double t) {

        final double omegaSign1;
        final double omegaSign2;
        final double omegaSign3;

        final CopolarN valuesN = jacobi.valuesN((t - tRef) * tScale);

        final double omega1 = o1Scale * valuesN.cn();
        final double omega2 = o2Scale * valuesN.sn();
        final double omega3 = o3Scale * valuesN.dn();

        final double condition = m2/twoE;

        final double cas;

        if (condition > i2C) {
            if (i1 < i2 && i2 < i3) {//case 1 : i1 < i2 < condition < i3
                omegaSign1 = 1.0;
                omegaSign2 = 1.0;
                omegaSign3 = 1.0;

                cas = 1;
                return new double[] {
                        omegaSign1 * omega1,
                        omegaSign2 * omega2,
                        omegaSign3 * omega3,
                        cas
                };
            }
            if (i1 < i3 && i3 < i2) {//case 3 : i1 < i3 < condition < i2

                omegaSign1 = 1.0;
                omegaSign2 = 1.0;
                omegaSign3 = 1.0;

                cas = 3;
                return new double[] {
                        omegaSign1 * omega1,
                        omegaSign2 * omega2,
                        omegaSign3 * omega3,
                        cas
                };
            }
            if (i2 < i1 && i1 < i3) {//case 5 : i2 < i1 < condition < i3
                omegaSign1 = -1.0;
                omegaSign2 = 1.0;
                omegaSign3 = -1.0;

                cas = 5;
                return new double[] {
                        omegaSign1 * omega1,
                        omegaSign2 * omega2,
                        omegaSign3 * omega3,
                        cas
                };
            }
            if (i2 < i3 && i3 < i1) {//case 7 : i2 < i3 < condition < i1
                omegaSign1 = 1.0;
                omegaSign2 = -1.0;
                omegaSign3 = -1.0;

                cas = 7;
                return new double[] {
                        omegaSign1 * omega2,
                        omegaSign2 * omega3,
                        omegaSign3 * omega1,
                        cas
                };
            }
            if (i3 < i1 && i1 < i2) {//case 9 : i3 < i1 < condition < i2
                omegaSign1 = -1.0;
                omegaSign2 = -1.0;
                omegaSign3 = 1.0;

                cas = 9;
                return new double[] {
                        omegaSign1 * omega3,
                        omegaSign2 * omega1,
                        omegaSign3 * omega2,
                        cas
                };
            }
            if (i3 < i2 && i2 < i1) {//case 11 : i3 < i2 < condition < i1
                omegaSign1 = -1.0;
                omegaSign2 = 1.0;
                omegaSign3 = -1.0;

                cas = 11;
                return new double[] {
                        omegaSign1 * omega1,
                        omegaSign2 * omega2,
                        omegaSign3 * omega3,
                        cas
                };
            }
        }

        if(condition < i2C) {
            if (i1 < i2 && i2 < i3) {//case 2 : i1 < condition < i2 < i3
                omegaSign1 = 1.0;
                omegaSign2 = -1.0;
                omegaSign3 = 1.0;

                cas = 2;
                return new double[] {
                        omegaSign1 * omega1,
                        omegaSign2 * omega2,
                        omegaSign3 * omega3,
                        cas
                };
            }
            if (i1 < i3 && i3 < i2) {//case 4 : i1 < condition < i3 < i2
                omegaSign1 = 1.0;
                omegaSign2 = 1.0;
                omegaSign3 = -1.0;

                cas = 4;
                return new double[] {
                        omegaSign1 * omega2,
                        omegaSign2 * omega3,
                        omegaSign3 * omega1,
                        cas
                };
            }
            if (i2 < i1 && i1 < i3) {//case 6 : i2 < condition < i1<  i3
                omegaSign1 = -1.0;
                omegaSign2 = 1.0;
                omegaSign3 = 1.0;

                cas = 6;
                return new double[] {
                        omegaSign1 * omega3,
                        omegaSign2 * omega1,
                        omegaSign3 * omega2,
                        cas
                };
            }
            if (i2 < i3 && i3 < i1) {//case 8 : i2 < condition < i3 < i1
                omegaSign1 = 1.0;
                omegaSign2 = -1.0;
                omegaSign3 = 1.0;

                cas = 8;
                return new double[] {
                        omegaSign1 * omega1,
                        omegaSign2 * omega2,
                        omegaSign3 * omega3,
                        cas
                };
            }
            if (i3 < i1 && i1 < i2) {//case 10 : i3 < condition < i1 < i2

                omegaSign1 = 1.0;
                omegaSign2 = -1.0;
                omegaSign3 = 1.0;

                cas = 10;
                return new double[] {
                        omegaSign1 * omega1,
                        omegaSign2 * omega2,
                        omegaSign3 * omega3,
                        cas
                };
            }
            if (i3 < i2 && i2 < i1) {//case 12 : i3 < condition < i2 < i1
                omegaSign1 = 1.0;
                omegaSign2 = 1.0;
                omegaSign3 = -1.0;

                cas = 12;
                return new double[] {
                        omegaSign1 * omega1,
                        omegaSign2 * omega2,
                        omegaSign3 * omega3,
                        cas
                };
            }

        }
        return new double[] {Double.NaN};
    }

    /** Compute a time offset to modify the temporal origin of the problem.
     *
     * @param t
     * @return tRef
     */
    private double getTRef(double t) {

        final double cas = omega(t)[3];
        if (cas == 1  || cas == 4 || cas == 6
                ||cas == 8 || cas == 10 ||cas == 11 ) {;
                return t0 - jacobi.arcsn(y0C[1] / o2Scale) / tScale;
        }

        if (cas == 3 || cas == 12 || cas == 2|| cas == 9){
            return t0 - jacobi.arccn(y0C[0] / o1Scale) / tScale;
        }

        if (cas == 5||cas == 7) {
            return t0 - jacobi.arcdn(y0C[2] / o3Scale) / tScale;
        }

        return 0;

    }

}//Fin du programme
