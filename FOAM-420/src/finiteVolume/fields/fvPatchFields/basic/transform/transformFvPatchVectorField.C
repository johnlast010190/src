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
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/basic/transform/transformFvPatchField.H"
#include "fields/Fields/symmTransformField/symmTransformField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
Foam::tmp<Foam::CoeffField<Foam::vector>>
Foam::transformFvPatchField<Foam::vector>::gradientInternalBCoeffs() const
{
    //To revisit
    typedef CoeffField<vector> TypeCoeffField;
    typedef typename TypeCoeffField::squareTypeField
        squareTypeField;

    tmp<TypeCoeffField> bct
    (
        new TypeCoeffField(this->size())
    );
    TypeCoeffField& bc = bct.ref();

    squareTypeField& bcSq = bc.asSquare();

    const fvPatch& p = this->patch();

    tmp<Field<vector>> snGradTrt(snGradTransformDiag());

    bcSq = -p.deltaCoeffs()*snGradTrt()*snGradTrt();

    return bct;
}


template<>
Foam::tmp<Foam::Field<Foam::vector>>
Foam::transformFvPatchField<Foam::vector>::gradientBoundaryBCoeffs() const
{
    typedef CoeffField<vector> TypeCoeffField;
    typedef typename TypeCoeffField::squareTypeField
        squareTypeField;

    tmp<TypeCoeffField> bct
    (
        gradientInternalBCoeffs()
    );
    TypeCoeffField& bc = bct.ref();

    squareTypeField& bcSq = bc.asSquare();

    return snGrad() - (bcSq&this->patchInternalField());
}


// ************************************************************************* //
