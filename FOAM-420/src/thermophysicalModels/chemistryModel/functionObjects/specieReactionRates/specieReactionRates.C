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
    (c) 2019-2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "functionObjects/specieReactionRates/specieReactionRates.H"
#include "fields/volFields/volFields.H"
#include "finiteVolume/fvc/fvcVolumeIntegrate.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchField.H"
#include "chemistryModel/chemistryModel/chemistryModel.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ChemistryModelType>
void Foam::functionObjects::specieReactionRates<ChemistryModelType>::
writeFileHeader
(
    Ostream& os
) const
{
    writeHeader(os, "Specie reaction rates");
    volRegion::writeFileHeader(*this, os);
    writeHeaderValue(os, "nSpecie", chemistryModel_.nSpecie());
    writeHeaderValue(os, "nReaction", chemistryModel_.nReaction());

    writeCommented(os, "Time");
    writeDelimited(os, "Reaction");

    const wordList& speciesNames =
        chemistryModel_.thermo().composition().species();

    forAll(speciesNames, si)
    {
        writeDelimited(os, speciesNames[si] + "_min");
        writeDelimited(os, speciesNames[si] + "_max");
        writeDelimited(os, speciesNames[si] + "_mean");
    }
    writeDelimited(os, "Qdot_min");
    writeDelimited(os, "Qdot_max");
    writeDelimited(os, "Qdot_mean");

    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ChemistryModelType>
Foam::functionObjects::specieReactionRates<ChemistryModelType>::
specieReactionRates
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    volRegion(fvMeshFunctionObject::mesh_, dict),
    writeFile(obr_, name, typeName, dict),
    chemistryModel_
    (
        fvMeshFunctionObject::obr().template lookupObject<ChemistryModelType>
        (
            IOobject::groupName
            (
                "chemistryProperties",
                dict.lookupOrDefault("phaseName", word::null)
            )
        )
    ),
    writeRRFields_
    (
        dict.lookupOrDefault("writeRRFields", List<Tuple2<word, label>>())
    ),
    writeQdotFields_
    (
        dict.lookupOrDefault("writeQdotFields", labelList())
    )
{
    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ChemistryModelType>
Foam::functionObjects::specieReactionRates<ChemistryModelType>::
~specieReactionRates()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ChemistryModelType>
bool Foam::functionObjects::specieReactionRates<ChemistryModelType>::read
(
    const dictionary& dict
)
{
    regionFunctionObject::read(dict);

    return true;
}


template<class ChemistryModelType>
bool Foam::functionObjects::specieReactionRates<ChemistryModelType>::execute()
{
    return true;
}


template<class ChemistryModelType>
bool Foam::functionObjects::specieReactionRates<ChemistryModelType>::write()
{
    const fvMesh& mesh = this->fvMeshFunctionObject::mesh_;
    const label nSpecie = chemistryModel_.nSpecie();
    const label nReaction = chemistryModel_.nReaction();

    // Region volume
    const scalar V = this->V();

    for (label ri=0; ri<nReaction; ri++)
    {
        writeTime(file());
        file() << token::TAB << ri;

        // Sum heat release for this reaction
        volScalarField::Internal Qdot
        (
            IOobject
            (
                "Qdot",
                mesh.time().timeName(),
                obr()
            ),
            mesh,
            dimensionedScalar("0", dimPower/dimVolume, scalar(0))
        );

        for (label si=0; si<nSpecie; si++)
        {
            volScalarField::Internal RR
            (
                chemistryModel_.calculateRR(ri, si)
            );

            word speciesName =
                chemistryModel_.thermo().composition().species()[si];

            scalar sumVRRi = 0;
            scalar minVRRi = GREAT;
            scalar maxVRRi = -GREAT;

            if (isNull(cellIDs()))
            {
                sumVRRi = fvc::domainIntegrate(RR).value();
                minVRRi = gMin(RR);
                maxVRRi = gMax(RR);
            }
            else
            {
                sumVRRi = gSum
                (
                    scalarField(mesh.V()*RR, cellIDs())
                );
                minVRRi =
                    gMin
                    (
                        scalarField
                        (
                            RR,
                            cellIDs()
                        )
                    );
                maxVRRi =
                    gMax
                    (
                        scalarField
                        (
                            RR,
                            cellIDs()
                        )
                    );
            }

            file() << token::TAB << minVRRi;
            file() << token::TAB << maxVRRi;
            file() << token::TAB << sumVRRi/V;

            Qdot += chemistryModel_.calculateQdot(ri, si);

            forAll(writeRRFields_, i)
            {
                if
                (
                    writeRRFields_[i]
                 == Tuple2<word, label>
                    (
                        speciesName,
                        ri
                    )
                )
                {
                    OStringStream name;
                    name << "RR_" << speciesName << "_" << ri;
                    if (!obr().template foundObject<volScalarField>(name.str()))
                    {
                        regIOobject::store
                        (
                            new volScalarField
                            (
                                IOobject
                                (
                                    name.str(),
                                    mesh.time().timeName(),
                                    obr(),
                                    IOobject::NO_READ,
                                    IOobject::AUTO_WRITE
                                ),
                                mesh,
                                RR.dimensions(),
                                zeroGradientFvPatchField<scalar>::typeName
                            )
                        );
                    }
                    volScalarField& RRField =
                        obr().template lookupObjectRef<volScalarField>(name.str());
                    RRField.ref() = RR;
                    RRField.correctBoundaryConditions();
                }
            }
        }

        // Write heat release
        scalar sumVQdoti = 0;
        scalar minVQdoti = GREAT;
        scalar maxVQdoti = -GREAT;

        if (isNull(cellIDs()))
        {
            minVQdoti = gMin(Qdot);
            maxVQdoti = gMax(Qdot);
            sumVQdoti = fvc::domainIntegrate(Qdot).value();
        }
        else
        {
            sumVQdoti = gSum
            (
                scalarField(fvMeshFunctionObject::mesh_.V()*Qdot, cellIDs())
            );
            minVQdoti =
                gMin
                (
                    scalarField
                    (
                        Qdot,
                        cellIDs()
                    )
                );
            maxVQdoti =
                gMax
                (
                    scalarField
                    (
                        Qdot,
                        cellIDs()
                    )
                );
        }

        file() << token::TAB << minVQdoti;
        file() << token::TAB << maxVQdoti;
        file() << token::TAB << sumVQdoti/V;

        forAll(writeQdotFields_, i)
        {
            if (writeQdotFields_[i] == ri)
            {
                OStringStream name;
                name << "Qdot_" << ri;
                if (!obr().template foundObject<volScalarField>(name.str()))
                {
                    regIOobject::store
                    (
                        new volScalarField
                        (
                            IOobject
                            (
                                name.str(),
                                mesh.time().timeName(),
                                obr(),
                                IOobject::NO_READ,
                                IOobject::AUTO_WRITE
                            ),
                            mesh,
                            Qdot.dimensions(),
                            zeroGradientFvPatchField<scalar>::typeName
                        )
                    );
                }
                volScalarField& QdotField =
                    obr().template lookupObjectRef<volScalarField>(name.str());
                QdotField.ref() = Qdot;
                QdotField.correctBoundaryConditions();
            }
        }

        file() << nl;
    }

    file() << endl;

    return true;
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "chemistryModel/rhoChemistryModel/rhoChemistryModel.H"
#include "chemistryModel/psiChemistryModel/psiChemistryModel.H"

namespace Foam
{
    typedef functionObjects::specieReactionRates<psiChemistryModel>
        psiSpecieReactionRates;

    defineTemplateTypeNameAndDebugWithName
    (
        psiSpecieReactionRates,
        "psiSpecieReactionRates",
        0
    );

    addToRunTimeSelectionTable
    (
        functionObject,
        psiSpecieReactionRates,
        dictionary
    );


    typedef functionObjects::specieReactionRates<rhoChemistryModel>
        rhoSpecieReactionRates;

    defineTemplateTypeNameAndDebugWithName
    (
        rhoSpecieReactionRates,
        "rhoSpecieReactionRates",
        0
    );

    addToRunTimeSelectionTable
    (
        functionObject,
        rhoSpecieReactionRates,
        dictionary
    );
}


// ************************************************************************* //
