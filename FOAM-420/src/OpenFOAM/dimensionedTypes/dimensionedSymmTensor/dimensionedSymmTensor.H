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

Typedef
    Foam::dimensionedSymmTensor

Description
    Dimensioned tensor obtained from generic dimensioned type.

SourceFiles
    dimensionedSymmTensor.C

\*---------------------------------------------------------------------------*/

#ifndef dimensionedSymmTensor_H
#define dimensionedSymmTensor_H

#include "dimensionedTypes/dimensionedVector/dimensionedVector.H"
#include "primitives/SymmTensor/symmTensor/symmTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef dimensioned<symmTensor> dimensionedSymmTensor;


// global functions

dimensionedSymmTensor sqr(const dimensionedVector&);
dimensionedSymmTensor innerSqr(const dimensionedSymmTensor&);

dimensionedScalar tr(const dimensionedSymmTensor&);
dimensionedSymmTensor symm(const dimensionedSymmTensor&);
dimensionedSymmTensor twoSymm(const dimensionedSymmTensor&);
dimensionedSymmTensor dev(const dimensionedSymmTensor&);
dimensionedSymmTensor dev2(const dimensionedSymmTensor&);
dimensionedScalar det(const dimensionedSymmTensor&);
dimensionedSymmTensor cof(const dimensionedSymmTensor&);
dimensionedSymmTensor inv(const dimensionedSymmTensor&);


// global operators

//- Hodge Dual operator (tensor -> vector)
dimensionedVector operator*(const dimensionedSymmTensor&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
