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
    (c) 2021-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "weightedAverageMixture.H"


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::weightedAverageMixture::weightedAverageMixtureTypeInternal
(
    const UPtrList<baseModels<Type>>& mod
) const
{
    const label fieldSize
    (
        frac_.fractions()
        [
            frac_.passiveIndex() == 0 ? 1 : 0
        ].primitiveField().size()
    );
    tmp<Field<Type>> tMixture(new Field<Type>(fieldSize, Zero));
    Field<Type>& mixture = tMixture.ref();
    scalarField fractionsSum(fieldSize, 0.0);

    forAll(mod, i)
    {
        if (frac_.active()[i] && i != frac_.passiveIndex())
        {
            mixture +=
                mod[i].primitiveField()*frac_.fractions()[i].primitiveField();
            fractionsSum += frac_.fractions()[i].primitiveField();
        }
    }
    if (frac_.passiveIndex() > -1)
    {
        mixture += mod[frac_.passiveIndex()].primitiveField()*(1-fractionsSum);
        return tMixture;
    }
    else
    {
        return tMixture/fractionsSum;
    }
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::weightedAverageMixture::weightedAverageMixtureTypeGeometric
(
    const UPtrList<baseModels<Type>>& mod
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> volField;
    tmp<volField> tMixture(mod[0]()*0.0);
    volField& mixture = tMixture.ref();
    volScalarField YSum(frac_.fractions()[0]*0.0);
    forAll(mod, i)
    {
        if (frac_.active()[i] && i != frac_.passiveIndex())
        {
            mixture += mod[i]()*frac_.fractions()[i];
            YSum += frac_.fractions()[i];
        }
    }
    if (frac_.passiveIndex() > -1)
    {
        mixture += mod[frac_.passiveIndex()]()*(1-YSum);
        return tMixture;
    }
    else
    {
        return tMixture/YSum;
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::weightedAverageMixture::weightedAverageMixtureTypePatch
(
    const UPtrList<baseModels<Type>>& mod,
    const label patchi
) const
{
    const label fieldSize
    (
        frac_.fractions()
        [
            frac_.passiveIndex() == 0 ? 1 : 0
        ].boundaryField()[patchi].size()
    );
    tmp<Field<Type>> tMixture(new Field<Type>(fieldSize, Zero));
    Field<Type>& mixture = tMixture.ref();
    scalarField fractionsSum(fieldSize, 0.0);
    forAll(mod, i)
    {
        if (frac_.active()[i] && i != frac_.passiveIndex())
        {
            mixture +=
                mod[i].boundaryField()[patchi]
               *frac_.fractions()[i].boundaryField()[patchi];
            fractionsSum += frac_.fractions()[i].boundaryField()[patchi];
        }
    }
    if (frac_.passiveIndex() > -1)
    {
        mixture +=
            mod[frac_.passiveIndex()].boundaryField()[patchi]*(1-fractionsSum);
        return tMixture;
    }
    else
    {
        return tMixture/fractionsSum;
    }
}


template<class Type>
Type Foam::weightedAverageMixture::weightedAverageMixtureTypeCell
(
    const UPtrList<baseModels<Type>>& mod,
    const label celli
) const
{
    Type mixture = Zero;
    scalar fractionsSum = 0;
    forAll(mod, i)
    {
        if (frac_.active()[i] && i != frac_.passiveIndex())
        {
            mixture += mod[i][celli]*frac_.fractions()[i][celli];
            fractionsSum += frac_.fractions()[i][celli];
        }
    }
    if (frac_.passiveIndex() > -1)
    {
        mixture += mod[frac_.passiveIndex()][celli]*(1-fractionsSum);
        return mixture;
    }
    else
    {
        return mixture/fractionsSum;
    }
}


// ************************************************************************* //
