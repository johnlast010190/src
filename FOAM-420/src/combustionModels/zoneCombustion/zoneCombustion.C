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
    (c) 2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "zoneCombustion/zoneCombustion.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::zoneCombustion<Type>::filter
(
    const tmp<fvScalarMatrix>& tR
) const
{
    fvScalarMatrix& R = tR.ref();
    scalarField& Su = R.source();
    scalarField filteredField(Su.size(), 0);

    forAll(zoneNames_, zonei)
    {
        const labelList& cells = this->mesh().cellZones()[zoneNames_[zonei]];

        forAll(cells, i)
        {
            filteredField[cells[i]] = Su[cells[i]];
        }
    }

    Su = filteredField;

    if (R.hasDiag())
    {
        scalarField& Sp = R.diag();

        forAll(zoneNames_, zonei)
        {
            const labelList& cells =
                this->mesh().cellZones()[zoneNames_[zonei]];

            forAll(cells, i)
            {
                filteredField[cells[i]] = Sp[cells[i]];
            }
        }

        Sp = filteredField;
    }

    return tR;
}


template<class Type>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::zoneCombustion<Type>::filter
(
    const tmp<volScalarField>& tS
) const
{
    scalarField& S = tS.ref();
    scalarField filteredField(S.size(), 0);

    forAll(zoneNames_, zonei)
    {
        const labelList& cells = this->mesh().cellZones()[zoneNames_[zonei]];

        forAll(cells, i)
        {
            filteredField[cells[i]] = S[cells[i]];
        }
    }

    S = filteredField;

    return tS;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::combustionModels::zoneCombustion<Type>::zoneCombustion
(
    const word& modelType,
    const objectRegistry& obr,
    const word& combustionProperties,
    const word& phaseName
)
:
    Type(modelType, obr, combustionProperties, phaseName),
    combustionModelPtr_
    (
        Type::New
        (
            obr,
            "zoneCombustionProperties",
            phaseName
        )
    ),
    zoneNames_(this->coeffs().lookup("zones"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::combustionModels::zoneCombustion<Type>::~zoneCombustion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
typename Type::ReactionThermo&
Foam::combustionModels::zoneCombustion<Type>::thermo()
{
    return combustionModelPtr_->thermo();
}


template<class Type>
const typename Type::ReactionThermo&
Foam::combustionModels::zoneCombustion<Type>::thermo() const
{
    return combustionModelPtr_->thermo();
}


template<class Type>
void Foam::combustionModels::zoneCombustion<Type>::correct()
{
    combustionModelPtr_->correct();
}


template<class Type>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::zoneCombustion<Type>::R(volScalarField& Y) const
{
    return filter(combustionModelPtr_->R(Y));
}


template<class Type>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::zoneCombustion<Type>::Qdot() const
{
    return filter(combustionModelPtr_->Qdot());
}


template<class Type>
bool Foam::combustionModels::zoneCombustion<Type>::read()
{
    if (Type::read())
    {
        combustionModelPtr_->read();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
