/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  FOAM (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.2.0
|    o     o     |  ESI Ltd. <http://esi.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of FOAMcore.
    FOAMcore is based on OpenFOAM (R) <http://www.openfoam.org/>.

    FOAMcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FOAMcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with FOAMcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2013-2020 Esi Ltd.
    (c) 2009 Jovani Favero
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2009 Hrv Jasak

\*---------------------------------------------------------------------------*/

#include "primitives/Tensor/tensor/tensor.H"
#include "primitives/polynomialEqns/cubicEqn/cubicEqn.H"
#include "global/constants/mathematical/mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::tensor::vsType::typeName = "tensor";

template<>
const char* const Foam::tensor::vsType::componentNames[] =
{
    "xx", "xy", "xz",
    "yx", "yy", "yz",
    "zx", "zy", "zz"
};

template<>
const Foam::tensor Foam::tensor::vsType::zero(tensor::uniform(0));

template<>
const Foam::tensor Foam::tensor::vsType::one(tensor::uniform(1));

template<>
const Foam::tensor Foam::tensor::vsType::max(tensor::uniform(VGREAT));

template<>
const Foam::tensor Foam::tensor::vsType::min(tensor::uniform(-VGREAT));

template<>
const Foam::tensor Foam::tensor::vsType::rootMax(tensor::uniform(ROOTVGREAT));

template<>
const Foam::tensor Foam::tensor::vsType::rootMin(tensor::uniform(-ROOTVGREAT));

template<>
const Foam::tensor Foam::tensor::I
(
    1, 0, 0,
    0, 1, 0,
    0, 0, 1
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::vector Foam::eigenValues(const tensor& T)
{
    if (T == tensor::zero)
    {
        return vector::zero;
    }

    scalar i = 0;
    scalar ii = 0;
    scalar iii = 0;

    //Solving for the Ax^3 + 3Bx^2 + 3Cx + D = 0 with A = 1;

    //Creating coefficients
    scalar B = -T.xx() - T.yy() - T.zz();

    B = B/3;

    scalar C = T.xx()*T.yy() + T.xx()*T.zz() + T.yy()*T.zz()
        - T.xy()*T.yx() - T.xz()*T.zx() - T.yz()*T.zy();

    C = C/3;

    scalar D = - T.xx()*T.yy()*T.zz() - T.xy()*T.yz()*T.zx()
        - T.xz()*T.yx()*T.zy() + T.xz()*T.yy()*T.zx()
        + T.xy()*T.yx()*T.zz() + T.xx()*T.yz()*T.zy();


    //normalize the tensor if it contains small values
    scalar maxN = max(mag(B),mag(C));
    maxN = max(mag(D), maxN);
    maxN = max(SMALL, maxN);

    if (maxN < 1)
    {
        B /= maxN;
        C /= sqr(maxN);
        D /= pow3(maxN);
    }

    //Calculating covariant coefficients
    scalar d1 = C - B*B;
    scalar d2 = D - B*C;
    scalar d3 = B*D - C*C;

    //Calculating Discriminant
    scalar Di = 4*d1*d3-d2*d2;

    if ((fabs(Di))<SMALL)
    {
        if (fabs(d1)<SMALL)
        {
            i = -B;
            ii = i;
            iii = i;
        }
        else
        {
            i = -0.5*d2/d1;
            ii = d2/d1 - 3*B;
            iii = i;
        }
    }
    else if (Di<0)
    {
        // Case with complex roots.
        //
        // Depress the polynomial
        scalar p;
        scalar q;
        scalar T0;
        scalar T1;
        scalar C_b;
        scalar D_b;
        if (pow3(B)*D >= pow3(C))
        {
            D_b = -2*B*d1 + d2;
            C_b = d1;
            T0 = -sign(D_b)*sqrt(-Di);
            T1 = -D_b + T0;
        }
        else
        {
            D_b = -D*d2 + 2*C*d3;
            C_b = d3;
            T0 = -sign(D_b)*fabs(D)*sqrt(-Di);
            T1 = -D_b + T0;
        }

        p = cbrt(T1/2);

        if (fabs(T1-T0)<ROOTVSMALL)
        {
            q = -p;
        }
        else
        {
            q = -C_b/p;
        }

        if (C_b <= 0)
        {
            i = p + q;
        }
        else
        {
            i = -D_b/(sqr(p)+sqr(q)+C_b);
        }

        ii = -i/2;
        iii = ii;

        scalar ic = i + C;
        scalar iic = ii + C;
        scalar iiic = iii + C;

        if
        (
            pow3(B)*D >= pow3(C)
            || mag(ic) < ROOTVSMALL
            || mag(iic) < ROOTVSMALL
            || mag(iiic) < ROOTVSMALL
        )
        {
            i = i - B;
            ii = ii - B;
            iii = iii - B;
        }
        else
        {
            i = -D/ic;
            ii = -D/iic;
            iii = -D/iiic;
        }
    }
    else
    {
        d1 = min(d1, scalar(0));
        d3 = min(d3, scalar(0));
        scalar D_b = -2*B*d1 + d2;
        scalar th = 1/3.*fabs(atan2(sqrt(Di), -D_b));
        scalar x1 = 2*sqrt(-d1)*cos(th);
        scalar x3 = -sqrt(-d1)*(cos(th)+sqrt(3.)*sin(th));
        if (x1 + x3 > 2*B)
        {
            iii = x1-B;
        }
        else
        {
            iii = x3-B;
        }
        D_b = -D*d2 + 2*C*d3;
        th = 1/3.*fabs(atan2(D*sqrt(Di), -D_b));
        x1 = 2*sqrt(-d3)*cos(th);
        x3 = -sqrt(-d3)*(cos(th)+sqrt(3.)*sin(th));

        if (x1 + x3 < 2*C)
        {
            i = -D/(x1+C);
        }
        else
        {
            i = -D/(x3+C);
        }
        ii = -((i+iii)*C + B*i*iii)/(B*(i+iii)+C);
    }

    if (maxN < 1)
    {
        i *= maxN;
        ii *= maxN;
        iii *= maxN;
    }
    // Sort the eigenvalues into ascending order
    if (i > ii)
    {
        Swap(i, ii);
    }

    if (ii > iii)
    {
        Swap(ii, iii);
    }

    if (i > ii)
    {
        Swap(i, ii);
    }

    return vector(i, ii, iii);

//    static label nWarnings = 0;
//    const label maxWarnings = 100;
//
//    // Coefficients of the characteristic cubic polynomial (a = 1)
//    const scalar b =
//      - t.xx() - t.yy() - t.zz();
//    const scalar c =
//        t.xx()*t.yy() + t.xx()*t.zz() + t.yy()*t.zz()
//      - t.xy()*t.yx() - t.yz()*t.zy() - t.zx()*t.xz();
//    const scalar d =
//      - t.xx()*t.yy()*t.zz()
//      - t.xy()*t.yz()*t.zx() - t.xz()*t.zy()*t.yx()
//      + t.xx()*t.yz()*t.zy() + t.yy()*t.zx()*t.xz() + t.zz()*t.xy()*t.yx();
//
//    // Solve
//    Roots<3> roots = cubicEqn(1, b, c, d).roots();
//
//    // Check the root types
//    vector lambda = vector::zero;
//    forAll(roots, i)
//    {
//        switch (roots.type(i))
//        {
//            case roots::real:
//                lambda[i] = roots[i];
//                break;
//            case roots::complex:
//                if (nWarnings < maxWarnings)
//                {
//                    WarningInFunction
//                        << "Complex eigenvalues detected for tensor: " << t
//                        << endl;
//                    nWarnings++;
//                }
//                lambda[i] = 0;
//                break;
//            case roots::posInf:
//                lambda[i] = VGREAT;
//                break;
//            case roots::negInf:
//                lambda[i] = - VGREAT;
//                break;
//            case roots::nan:
//                FatalErrorInFunction
//                    << "Eigenvalue calculation failed for tensor: " << t
//                    << exit(FatalError);
//        }
//    }
//
//    // Sort the eigenvalues into ascending order
//    if (lambda.x() > lambda.y())
//    {
//        Swap(lambda.x(), lambda.y());
//    }
//    if (lambda.y() > lambda.z())
//    {
//        Swap(lambda.y(), lambda.z());
//    }
//    if (lambda.x() > lambda.y())
//    {
//        Swap(lambda.x(), lambda.y());
//    }
//
//    return lambda;
}


Foam::vector Foam::eigenVector
(
    const tensor& T,
    const scalar lambda,
    const vector& direction1,
    const vector& direction2
)
{
    // Construct the linear system for this eigenvalue
    tensor A(T - lambda*I);

    // Determinants of the 2x2 sub-matrices used to find the eigenvectors
    scalar sd0, sd1, sd2;
    scalar magSd0, magSd1, magSd2;

    // Sub-determinants for a unique eivenvalue
    sd0 = A.yy()*A.zz() - A.yz()*A.zy();
    sd1 = A.zz()*A.xx() - A.zx()*A.xz();
    sd2 = A.xx()*A.yy() - A.xy()*A.yx();
    magSd0 = mag(sd0);
    magSd1 = mag(sd1);
    magSd2 = mag(sd2);

    // Evaluate the eigenvector using the largest sub-determinant
    if (magSd0 > magSd1 && magSd0 > magSd2 && magSd0 > ROOTVSMALL)
    {
        vector ev
        (
            1,
            (A.yz()*A.zx() - A.zz()*A.yx())/sd0,
            (A.zy()*A.yx() - A.yy()*A.zx())/sd0
        );

        return ev/mag(ev);
    }
    else if (magSd1 > magSd2 && magSd1 > ROOTVSMALL)
    {
        vector ev
        (
            (A.xz()*A.zy() - A.zz()*A.xy())/sd1,
            1,
            (A.zx()*A.xy() - A.xx()*A.zy())/sd1
        );

        return ev/mag(ev);
    }
    else if (magSd2 > ROOTVSMALL)
    {
        vector ev
        (
            (A.xy()*A.yz() - A.yy()*A.xz())/sd2,
            (A.yx()*A.xz() - A.xx()*A.yz())/sd2,
            1
        );

        return ev/mag(ev);
    }

    // Sub-determinants for a repeated eigenvalue
    sd0 = A.yy()*direction1.z() - A.yz()*direction1.y();
    sd1 = A.zz()*direction1.x() - A.zx()*direction1.z();
    sd2 = A.xx()*direction1.y() - A.xy()*direction1.x();
    magSd0 = mag(sd0);
    magSd1 = mag(sd1);
    magSd2 = mag(sd2);

    // Evaluate the eigenvector using the largest sub-determinant
    if (magSd0 > magSd1 && magSd0 > magSd2 && magSd0 > ROOTVSMALL)
    {
        vector ev
        (
            1,
            (A.yz()*direction1.x() - direction1.z()*A.yx())/sd0,
            (direction1.y()*A.yx() - A.yy()*direction1.x())/sd0
        );

        return ev/mag(ev);
    }
    else if (magSd1 > magSd2 && magSd1 > ROOTVSMALL)
    {
        vector ev
        (
            (direction1.z()*A.zy() - A.zz()*direction1.y())/sd1,
            1,
            (A.zx()*direction1.y() - direction1.x()*A.zy())/sd1
        );

        return ev/mag(ev);
    }
    else if (magSd2 > ROOTVSMALL)
    {
        vector ev
        (
            (A.xy()*direction1.z() - direction1.y()*A.xz())/sd2,
            (direction1.x()*A.xz() - A.xx()*direction1.z())/sd2,
            1
        );

        return ev/mag(ev);
    }

    // Triple eigenvalue
    return direction1^direction2;
}


Foam::tensor Foam::eigenVectors(const tensor& T, const vector& lambdas)
{
    vector Ux(1, 0, 0), Uy(0, 1, 0), Uz(0, 0, 1);

    Ux = eigenVector(T, lambdas.x(), Uy, Uz);
    Uy = eigenVector(T, lambdas.y(), Uz, Ux);
    Uz = eigenVector(T, lambdas.z(), Ux, Uy);

    return tensor(Ux, Uy, Uz);
}


Foam::tensor Foam::eigenVectors(const tensor& T)
{
    const vector lambdas(eigenValues(T));

    return eigenVectors(T, lambdas);
}


Foam::vector Foam::eigenValues(const symmTensor& T)
{
    return eigenValues(tensor(T));
}


Foam::vector Foam::eigenVector
(
    const symmTensor& T,
    const scalar lambda,
    const vector& direction1,
    const vector& direction2
)
{
    return eigenVector(tensor(T), lambda, direction1, direction2);
}


Foam::tensor Foam::eigenVectors(const symmTensor& T, const vector& lambdas)
{
    return eigenVectors(tensor(T), lambdas);
}


Foam::tensor Foam::eigenVectors(const symmTensor& T)
{
    return eigenVectors(tensor(T));
}


// Matrix inversion with singular value decomposition
Foam::tensor Foam::hinv(const tensor& t)
{
    static const scalar large = 1e14;
    static const scalar small = 1e-14;

    scalar dett(det(t));
    if (dett > small)
    {
        return inv(t, dett);
    }
    else
    {
        //check to filter out small tensors
        scalar tmagPow3 =
        sign(t.xx())*t.xx()*t.xx()*t.xx()+
        sign(t.xy())*t.xy()*t.xy()*t.xy()+
        sign(t.xz())*t.xz()*t.xz()*t.xz()+
        sign(t.yx())*t.yx()*t.yx()*t.yx()+
        sign(t.yy())*t.yy()*t.yy()*t.yy()+
        sign(t.yz())*t.yz()*t.yz()*t.yz()+
        sign(t.zx())*t.zx()*t.zx()*t.zx()+
        sign(t.zy())*t.zy()*t.zy()*t.zy()+
        sign(t.zz())*t.zz()*t.zz()*t.zz();

        if (dett > tmagPow3*small)
        {
            return inv(t, dett);
        }
        else //degenerate
        {
            vector eig = eigenValues(t);
            tensor eigVecs = eigenVectors(t);

            tensor zeroInv(tensor::zero);

            // Test if all eigen values are zero.
            // If this happens then eig.z() = SMALL, and hinv(t)
            // returns a zero tensor.
            // Jovani Favero, 18/Nov/2009
            // Further bug fix: replace > with == and add SMALL to zeroInv
            // Dominik Christ, 7/Aug/2012
            if (mag(eig.z()) == large*mag(eig.z()))
            {
                // Three zero eigen values (z is largest in magnitude).
                // Return zero inverse
                return zeroInv;
            }

            // Compare smaller eigen values and if out of range of large
            // consider them singular

            if (mag(eig.z()) > large*mag(eig.x()))
            {
                // Make a tensor out of symmTensor sqr.  HJ, 24/Oct/2009
                zeroInv += tensor(sqr(eigVecs.x()));
            }

            if (mag(eig.z()) > large*mag(eig.y()))
            {
                // Make a tensor out of symmTensor sqr.  HJ, 24/Oct/2009
                zeroInv += tensor(sqr(eigVecs.y()));
            }

            return inv(t + zeroInv) - zeroInv;
        }
    }
}


Foam::symmTensor Foam::hinv(const symmTensor& t)
{
    static const scalar large = 1e14;
    static const scalar small = 1e-14;

    scalar dett(det(t));
    if (dett > small)
    {
        return inv(t, dett);
    }
    else
    {
        //check to filter out small tensors
        //by scaling with tensor magnitude
        //check to filter out small tensors
        scalar tmagPow3 =
        sign(t.xx())*t.xx()*t.xx()*t.xx()+
        sign(t.xy())*t.xy()*t.xy()*t.xy()+
        sign(t.xz())*t.xz()*t.xz()*t.xz()+
        sign(t.yy())*t.yy()*t.yy()*t.yy()+
        sign(t.yz())*t.yz()*t.yz()*t.yz()+
        sign(t.zz())*t.zz()*t.zz()*t.zz();

        if (dett > tmagPow3*small)
        {
            return inv(t, dett);
        }
        else //degenerate
        {
            vector eig = eigenValues(t);
            tensor eigVecs = eigenVectors(t);

            symmTensor zeroInv(symmTensor::zero);

            // Test if all eigen values are zero,
            // If this happens then eig.z() = SMALL
            // and hinv(t) return a zero tensor.
            // Jovani Favero, 18/Nov/2009
            // Further bug fix: replace > with == and add SMALL to zeroInv
            // Dominik Christ, 7/Aug/2012
            if (mag(eig.z()) == large*mag(eig.z()))
            {
                // Three zero eigen values (z is largest in magnitude).
                // Return zero inverse
                return zeroInv;
            }

            // Compare smaller eigen values and if out of range of large
            // consider them singular

            if (mag(eig.z()) > large*mag(eig.x()))
            {
                zeroInv += sqr(eigVecs.x());
            }

            if (mag(eig.z()) > large*mag(eig.y()))
            {
                zeroInv += sqr(eigVecs.y());
            }

            return inv(t + zeroInv) - zeroInv;
        }
    }
}


// ************************************************************************* //
