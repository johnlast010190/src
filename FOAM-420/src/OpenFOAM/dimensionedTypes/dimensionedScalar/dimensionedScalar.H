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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2017-2022 Esi Ltd.

Typedef
    Foam::dimensionedScalar

Description
    Dimensioned scalar obtained from generic dimensioned type.

SourceFiles
    dimensionedScalar.C

\*---------------------------------------------------------------------------*/

#ifndef dimensionedScalar_H
#define dimensionedScalar_H

#include "dimensionedTypes/dimensionedType/dimensionedType.H"
#include "primitives/Scalar/scalar/scalar.H"
#include "dimensionedTypes/dimensionedScalar/dimensionedScalarFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensionedScalar operator+(const dimensionedScalar&, const scalar);
dimensionedScalar operator+(const scalar, const dimensionedScalar&);

dimensionedScalar operator-(const dimensionedScalar&, const scalar);
dimensionedScalar operator-(const scalar, const dimensionedScalar&);

dimensionedScalar operator*(const dimensionedScalar&, const scalar);
dimensionedScalar operator/(const scalar, const dimensionedScalar&);

dimensionedScalar pow(const dimensionedScalar&, const dimensionedScalar&);
dimensionedScalar pow3(const dimensionedScalar&);
dimensionedScalar pow4(const dimensionedScalar&);
dimensionedScalar pow5(const dimensionedScalar&);
dimensionedScalar pow6(const dimensionedScalar&);
dimensionedScalar pow025(const dimensionedScalar&);

dimensionedScalar sqrt(const dimensionedScalar&);
dimensionedScalar cbrt(const dimensionedScalar&);
dimensionedScalar hypot(const dimensionedScalar&, const dimensionedScalar&);

dimensionedScalar sign(const dimensionedScalar&);
dimensionedScalar sign101(const dimensionedScalar&); //- dcombest 8/15/2012
//dimensionedScalar sign101(const dimensionedScalar&, const dimensionedScalar&);
dimensionedScalar pos(const dimensionedScalar&);
dimensionedScalar pos0(const dimensionedScalar&);
dimensionedScalar neg(const dimensionedScalar&);
dimensionedScalar neg0(const dimensionedScalar&);
dimensionedScalar ptanh(const dimensionedScalar&); //- edevilliers 3/4/2013
dimensionedScalar posPart(const dimensionedScalar&);
dimensionedScalar negPart(const dimensionedScalar&);

dimensionedScalar exp(const dimensionedScalar&);
dimensionedScalar log(const dimensionedScalar&);
dimensionedScalar log10(const dimensionedScalar&);
dimensionedScalar sin(const dimensionedScalar&);
dimensionedScalar cos(const dimensionedScalar&);
dimensionedScalar tan(const dimensionedScalar&);
dimensionedScalar asin(const dimensionedScalar&);
dimensionedScalar acos(const dimensionedScalar&);
dimensionedScalar atan(const dimensionedScalar&);
dimensionedScalar atan2(const dimensionedScalar&, const dimensionedScalar&);
dimensionedScalar sinh(const dimensionedScalar&);
dimensionedScalar cosh(const dimensionedScalar&);
dimensionedScalar tanh(const dimensionedScalar&);
dimensionedScalar asinh(const dimensionedScalar&);
dimensionedScalar acosh(const dimensionedScalar&);
dimensionedScalar atanh(const dimensionedScalar&);
dimensionedScalar erf(const dimensionedScalar&);
dimensionedScalar erfc(const dimensionedScalar&);
dimensionedScalar lgamma(const dimensionedScalar&);
dimensionedScalar j0(const dimensionedScalar&);
dimensionedScalar j1(const dimensionedScalar&);
dimensionedScalar jn(const int, const dimensionedScalar&);
dimensionedScalar y0(const dimensionedScalar&);
dimensionedScalar y1(const dimensionedScalar&);
dimensionedScalar yn(const int, const dimensionedScalar&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
