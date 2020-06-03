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
package org.hipparchus.analysis.differentiation;

import org.hipparchus.RealFieldElement;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;

/** Abstract class representing both the value and the differentials of a function.
 * @param <S> the type of the field elements
 * @param <T> the type of the function derivative
 * @since 1.7
 */
public abstract class FieldUnivariateDerivative<S extends RealFieldElement<S>, T extends FieldUnivariateDerivative<S, T>>
    implements FieldDerivative<S, T> {

    /** {@inheritDoc} */
    @Override
    public int getFreeParameters() {
        return 1;
    }

    /** {@inheritDoc} */
    @Override
    public S getPartialDerivative(final int ... orders) throws MathIllegalArgumentException {
        if (orders.length != 1) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   orders.length, 1);
        }
        return getDerivative(orders[0]);
    }

    /** Get a derivative from the univariate derivative.
     * @param n derivation order (must be between 0 and {@link #getOrder()}, both inclusive)
     * @return n<sup>th</sup> derivative
     * @exception MathIllegalArgumentException if n is
     * either negative or strictly larger than {@link #getOrder()}
     */
    public abstract S getDerivative(int n) throws MathIllegalArgumentException;

    /** Convert the instance to a {@link DerivativeStructure}.
     * @return derivative structure with same value and derivative as the instance
     */
    public abstract FieldDerivativeStructure<S> toDerivativeStructure();

}
