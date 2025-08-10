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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2021-2023 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "speciesMassFractions.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(speciesMassFractions, 0);
}

const Foam::fvMesh& Foam::speciesMassFractions::getMesh(const objectRegistry& obr) const
{
    if (isA<fvSolutionRegistry>(obr))
    {
        return dynamic_cast<const fvSolutionRegistry&>(obr).mesh();
    }
    return refCast<const fvMesh>(obr);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void Foam::speciesMassFractions::createFields
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const objectRegistry& obr,
    const word& phaseName
)
{
    word YdefaultName(IOobject::groupName("Ydefault", phaseName));
    tmp<volScalarField> tYdefault;
    const fvMesh& mesh = getMesh(obr);
    forAll(species_, i)
    {
        IOobject header
        (
            IOobject::groupName(species_[i], phaseName),
            obr.time().timeName(),
            obr,
            IOobject::NO_READ
        );

        // Check if field exists and can be looked up or read
        word YName(IOobject::groupName(species_[i], phaseName));
        if
        (
            obr.foundObject<volScalarField>(YName)
         && obr.lookupObject<volScalarField>(YName).ownedByRegistry()
        )
        {
            // Take ownership
            volScalarField& Y = obr.lookupObjectRef<volScalarField>(YName);
            Y.release();
            Y_.set(i, &Y);
        }
        else if (header.typeHeaderOk<volScalarField>(true))
        {
            DebugInformation
                << "speciesMassFractions: reading " << species_[i]
                << endl;

            Y_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        YName,
                        obr.time().timeName(),
                        obr,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        }
        else
        {
            // Read Ydefault if not already read
            if (!tYdefault.valid())
            {
                if (!obr.foundObject<volScalarField>(YdefaultName))
                {
                    IOobject timeIO
                    (
                        YdefaultName,
                        obr.time().timeName(),
                        obr,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    );

                    IOobject constantIO
                    (
                        YdefaultName,
                        obr.time().constant(),
                        obr,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    );

                    IOobject time0IO
                    (
                        YdefaultName,
                        Time::timeName(0),
                        obr,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    );

                    if (timeIO.typeHeaderOk<volScalarField>(true))
                    {
                        tYdefault = new volScalarField(timeIO, mesh);
                    }
                    else if (constantIO.typeHeaderOk<volScalarField>(true))
                    {
                        tYdefault = new volScalarField(constantIO, mesh);
                    }
                    else
                    {
                        tYdefault = new volScalarField(time0IO, mesh);
                    }
                }
            }

            Y_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName(species_[i], phaseName),
                        obr.time().timeName(),
                        obr,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    obr.lookupObject<volScalarField>(YdefaultName)
                )
            );
        }
    }

    // Do not enforce constraint of sum of mass fractions to equal 1 here
    // - not applicable to all models
}


Foam::speciesMassFractions::speciesMassFractions
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const objectRegistry& obr,
    const word& phaseName
)
:
    volumeMassFractions(obr, phaseName),
    species_(specieNames),
    active_(species_.size(), true),
    Y_(species_.size())
{
    createFields(thermoDict, specieNames, obr, phaseName);
}



Foam::speciesMassFractions::speciesMassFractions
(
    const dictionary& thermoDict,
    const objectRegistry& obr,
    const word& phaseName
)
:
    volumeMassFractions(obr, phaseName),
    species_
    (
        (phaseName != word::null)
      ? thermoDict.optionalSubDict
        (
            phaseName
        ).lookupOrDefault<wordList>("species", {})
      : thermoDict.lookupOrDefault<wordList>("species", {})
    ),
    active_(species_.size(), true),
    Y_(species_.size())
{
    if (species_.size())
    {
        createFields(thermoDict, species_, obr, phaseName);
    }
}

// ************************************************************************* //
