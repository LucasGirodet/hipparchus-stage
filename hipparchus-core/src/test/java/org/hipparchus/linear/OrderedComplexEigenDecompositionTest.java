/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.linear;

import org.junit.Assert;

import org.hipparchus.complex.Complex;
import org.junit.Test;

public class OrderedComplexEigenDecompositionTest {

    private final RealMatrix A = MatrixUtils.createRealMatrix(new double[][] {
        { 3, -2 }, { 4, -1 } });

    @Test
    public void testDefinition() {
        OrderedComplexEigenDecomposition eigenDecomp = new OrderedComplexEigenDecomposition(A);

        FieldMatrix<Complex> A = MatrixUtils.createFieldMatrix(new Complex[][] {
            { new Complex(3, 0), new Complex(-2, 0) },
            { new Complex(4, 0), new Complex(-1, 0) } });

        // testing AV = lamba V - [0]
        Assert.assertEquals(A.operate(eigenDecomp.getEigenvector(0)),
                            eigenDecomp.getEigenvector(0).mapMultiply(eigenDecomp.getEigenvalues()[0]));

        // testing AV = lamba V - [1]
        Assert.assertEquals(A.operate(eigenDecomp.getEigenvector(1)),
                            eigenDecomp.getEigenvector(1).mapMultiply(eigenDecomp.getEigenvalues()[1]));

        // checking definition of the decomposition A*V = V*D
        Assert.assertEquals(A.multiply(eigenDecomp.getV()),
                            eigenDecomp.getV().multiply(eigenDecomp.getD()));
    }

    @Test
    public void testEqualEigenValues() {
        double[][] d = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
        Array2DRowRealMatrix matrix = new Array2DRowRealMatrix(d);
        ComplexEigenDecomposition ed = new OrderedComplexEigenDecomposition(matrix);
        for (Complex z : ed.getEigenvalues()) {
            Assert.assertEquals(1.0, z.getReal(),      1.0e-15);
            Assert.assertEquals(0.0, z.getImaginary(), 1.0e-15);
        }
    }

}
