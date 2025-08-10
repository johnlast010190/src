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

#include "functionObjects/moleFractions/moleFractions.H"
#include "basicThermo/basicThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
void Foam::moleFractions<ThermoType>::calculateMoleFractions()
{
    const ThermoType& thermo =
        obr().template lookupObject<ThermoType>
        (
            IOobject::groupName(basicThermo::dictName, phaseName_)
        );

    const PtrList<volScalarField>& Y = thermo.composition().Y();

    const volScalarField W(thermo.composition().W());

    forAll(Y, i)
    {
        X_[i] = W*Y[i]/thermo.composition().W(i);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::moleFractions<ThermoType>::moleFractions
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    phaseName_(dict.lookupOrDefault("phaseName", word::null))
{
    const word dictName =
        IOobject::groupName(basicThermo::dictName, phaseName_);
    if (obr().template foundObject<ThermoType>(dictName))
    {
        const ThermoType& thermo =
            obr().template lookupObject<ThermoType>(dictName);

        const PtrList<volScalarField>& Y = thermo.composition().Y();

        X_.setSize(Y.size());

        forAll(Y, i)
        {
            X_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "X_" + Y[i].name(),
                        mesh_.time().timeName(),
                        obr(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar("X", dimless, 0)
                )
            );
        }

        calculateMoleFractions();
    }
    else
    {
        FatalErrorInFunction
            << "Cannot find thermodynamics model of type "
            << ThermoType::typeName
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::moleFractions<ThermoType>::~moleFractions()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
bool Foam::moleFractions<ThermoType>::read
(
    const dictionary& dict
)
{
    return true;
}


template<class ThermoType>
bool Foam::moleFractions<ThermoType>::execute()
{
    calculateMoleFractions();
    return true;
}


template<class ThermoType>
bool Foam::moleFractions<ThermoType>::write()
{
    return true;
}


// ************************************************************************* //
