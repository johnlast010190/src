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
    (c) 2021-2022 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/constraint/nonConformalError/nonConformalErrorFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::nonConformalErrorFvPatchField<Type>::nonConformalErrorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    zeroGradientFvPatchField<Type>(p, iF),
    nonConformalErrorFvPatch_(refCast<const nonConformalErrorFvPatch>(p))
{}


template<class Type>
Foam::nonConformalErrorFvPatchField<Type>::nonConformalErrorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchField<Type>(p, iF, dict),
    nonConformalErrorFvPatch_(refCast<const nonConformalErrorFvPatch>(p))
{}


template<class Type>
Foam::nonConformalErrorFvPatchField<Type>::nonConformalErrorFvPatchField
(
    const nonConformalErrorFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchField<Type>(ptf, p, iF, mapper),
    nonConformalErrorFvPatch_(refCast<const nonConformalErrorFvPatch>(p))
{}


template<class Type>
Foam::nonConformalErrorFvPatchField<Type>::nonConformalErrorFvPatchField
(
    const nonConformalErrorFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    zeroGradientFvPatchField<Type>(ptf, iF),
    nonConformalErrorFvPatch_(ptf.nonConformalErrorFvPatch_)
{}


// ************************************************************************* //
