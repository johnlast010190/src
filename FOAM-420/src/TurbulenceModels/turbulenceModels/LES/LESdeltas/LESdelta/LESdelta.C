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

\*---------------------------------------------------------------------------*/

#include "LESdelta.H"
#include "fields/fvPatchFields/basic/calculated/calculatedFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LESdelta, 0);
    defineRunTimeSelectionTable(LESdelta, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LESdelta::LESdelta
(
    const word& name,
    const turbulenceModel& turbulence
)
:
    turbulenceModel_(turbulence),
    delta_
    (
        IOobject
        (
            name,
            turbulence.mesh().time().timeName(),
            turbulence.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        turbulence.mesh(),
        dimensionedScalar(name, dimLength, SMALL),
        calculatedFvPatchScalarField::typeName
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::LESdelta> Foam::LESdelta::New
(
    const word& name,
    const turbulenceModel& turbulence,
    const dictionary& dict
)
{
    const word deltaType(dict.lookup("delta"));

    Info<< "Selecting LES delta type " << deltaType << endl;

    const auto ctor = ctorTableLookup("LESdelta type", dictionaryConstructorTable_(), deltaType);
    return autoPtr<LESdelta>(ctor(name, turbulence, dict));
}


Foam::autoPtr<Foam::LESdelta> Foam::LESdelta::New
(
    const word& name,
    const turbulenceModel& turbulence,
    const dictionary& dict,
    const dictionaryConstructorTable& additionalConstructors
)
{
    const word deltaType(dict.lookup("delta"));

    Info<< "Selecting LES delta type " << deltaType << endl;

    // First on additional ones
    dictionaryConstructorTable::const_iterator cstrIter =
        additionalConstructors.find(deltaType);

    if (cstrIter != additionalConstructors.end())
    {
        return autoPtr<LESdelta>(cstrIter->second(name, turbulence, dict));
    }
    else
    {
        const auto ctor = ctorTableLookup("LESdelta type", dictionaryConstructorTable_(), deltaType);
        return autoPtr<LESdelta>(ctor(name, turbulence, dict));
    }
}


// ************************************************************************* //
