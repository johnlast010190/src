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
    (c) 2023 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "phaseVolumeFractions.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseVolumeFractions, 0);
}

const Foam::fvMesh& Foam::phaseVolumeFractions::getMesh
(
    const objectRegistry& obr
) const
{
    if (isA<fvSolutionRegistry>(obr))
    {
        return dynamic_cast<const fvSolutionRegistry&>(obr).mesh();
    }
    return refCast<const fvMesh>(obr);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void Foam::phaseVolumeFractions::createFields
(
    const dictionary& materialPropertiesDict,
    const objectRegistry& obr
)
{
    const fvMesh& mesh = getMesh(obr);

    word passivePhase =
        materialPropertiesDict.lookupOrDefault
        (
            "passivePhase",
            phases_.size() == 2 ? phases_[phases_.size()-1] : word::null
        );
    if (passivePhase == "none")
    {
        passivePhase = word::null;
    }
    if (passivePhase != word::null)
    {
        passivePhaseIndex_ = phases_.find(passivePhase);
        if (passivePhaseIndex_ == -1)
        {
            FatalErrorInFunction
                << "Passive phase " << passivePhase
                << " not found in available phases "
                << phases_
                << exit(FatalError);
        }
    }

    autoPtr<volScalarField> sumAlphas;
    forAll(phases_, phasei)
    {
        // Check if field exists and can be looked up or read
        const word alphaName(IOobject::groupName("alpha", phases_[phasei]));
        IOobject header
        (
            alphaName,
            obr.time().timeName(),
            obr,
            IOobject::NO_READ
        );
        dgdts_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("dgdt", phases_[phasei]),
                    obr.time().timeName(),
                    obr,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("0", dimless/dimTime, 0)
            )
        );
        if
        (
            obr.foundObject<volScalarField>(alphaName)
         && obr.lookupObject<volScalarField>(alphaName).ownedByRegistry()
        )
        {
            // Take ownership
            volScalarField& alpha =
                obr.lookupObjectRef<volScalarField>(alphaName);
            alpha.release();
            alphas_.set(phasei, &alpha);
        }
        else if (header.typeHeaderOk<volScalarField>(true))
        {
            DebugInformation
                << "phaseVolumeFractions: reading " << phases_[phasei]
                << endl;

            alphas_.set
            (
                phasei,
                new volScalarField
                (
                    IOobject
                    (
                        alphaName,
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
            alphas_.set
            (
                phasei,
                new volScalarField
                (
                    IOobject
                    (
                        alphaName,
                        obr.time().timeName(),
                        obr,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    dimensionedScalar("alpha", dimless, phasei == 0 ? 1 : 0)
                )
            );
        }

        if (phasei != passivePhaseIndex_)
        {
            if (!sumAlphas.valid())
            {
                sumAlphas.set
                (
                    new volScalarField(alphas_[phasei])
                );
            }
            else
            {
                sumAlphas() += alphas_[phasei];
            }
        }
    }

    if (passivePhaseIndex_ != -1)
    {
        alphas_[passivePhaseIndex_] = 1-sumAlphas();
    }

    if (phases_.size() > 2)
    {
        combinedAlphas_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "alphas",
                    obr.time().timeName(),
                    obr,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                getMesh(obr),
                dimensionedScalar("alphas", dimless, 0.0)
            )
        );
    }
}


Foam::phaseVolumeFractions::phaseVolumeFractions
(
    const dictionary& materialPropertiesDict,
    const wordList& phaseNames,
    const objectRegistry& obr,
    const word& dummy
)
:
    volumeMassFractions(obr, word::null),
    phases_(phaseNames),
    active_(phases_.size(), true),
    passivePhaseIndex_(-1),
    alphas_(phases_.size()),
    dgdts_(phases_.size())
{
    createFields(materialPropertiesDict, obr);
}


Foam::phaseVolumeFractions::phaseVolumeFractions
(
    const dictionary& materialPropertiesDict,
    const objectRegistry& obr,
    const word& dummy
)
:
    volumeMassFractions(obr, word::null),
    phases_(materialPropertiesDict.lookupOrDefault<wordList>("phases", {})),
    active_(phases_.size(), true),
    passivePhaseIndex_(-1),
    alphas_(phases_.size()),
    dgdts_(phases_.size())
{
    if (phases_.size())
    {
        createFields(materialPropertiesDict, obr);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::phaseVolumeFractions::recomputeCombinedAlphas()
{
    if (phases_.size() <= 2)
    {
        return;
    }

    scalar level = 1.0;
    combinedAlphas_().forceAssign(0.0);

    label excludedIndex = 0;
    if (passivePhaseIndex_ > -1)
    {
        excludedIndex = passivePhaseIndex_;
    }
    forAll(phases_, phasei)
    {
        if (phasei != excludedIndex)
        {
            combinedAlphas_() += level*alphas_[phasei];
            level += 1.0;
        }
    }
}

// ************************************************************************* //
