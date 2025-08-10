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
    (c) 2019-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "powerLoss/powerLossObjectiveFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/derived/resistivePressure/resistivePressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(powerLossObjectiveFunctionObject, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        powerLossObjectiveFunctionObject,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::powerLossObjectiveFunctionObject::
powerLossObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict,
    bool useAdjointFileFormat
)
:
    objectiveFunctionObject(runTime, const_cast<dictionary&>(objectiveDict)),
    areaObjPatchInlet_(objectivePatchInletArea()),
    areaObjPatchOutlet_(objectivePatchOutletArea()),
    adiabatic_(objectiveDict.lookupOrDefault<Switch>("adiabatic", false)),
    stressWork_(objectiveDict.lookupOrDefault<Switch>("stressWork", false)),
    compressibleForm_
    (
         objectiveDict.lookupOrDefault<word>
         (
             "compressibleForm",
            "pressureKinetic"
         )
    )
{
    createFiles(useAdjointFileFormat);
}


Foam::functionObjects::powerLossObjectiveFunctionObject::
powerLossObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict
)
:
    powerLossObjectiveFunctionObject(name, runTime, objectiveDict, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::powerLossObjectiveFunctionObject::
~powerLossObjectiveFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::functionObjects::powerLossObjectiveFunctionObject::read
(
    const dictionary& dict
)
{
    objectiveFunctionObject::read(dict);

    areaObjPatchInlet_ = objectivePatchInletArea();
    areaObjPatchOutlet_ = objectivePatchOutletArea();
    adiabatic_ = dict.lookupOrDefault<Switch>("adiabatic", false);
    stressWork_ = dict.lookupOrDefault<Switch>("stressWork", false);

    return true;
}


bool
Foam::functionObjects::powerLossObjectiveFunctionObject::execute()
{
    scalar powerLoss = 0;
    scalar mass = 0;

    scalar powerIn = 0;
    scalar powerOut = 0;

    tmp<volSymmTensorField> tdevRhoReff = devRhoReff();
    const volSymmTensorField::Boundary& devRhoReffb =
        tdevRhoReff().boundaryField();

    const surfaceVectorField::Boundary& Sfb =
        mesh_.Sf().boundaryField();

    const volScalarField* h(nullptr);

    if (compressible())
    {
        if (adiabatic_)
        {
            tmp<volScalarField> hh = p()/rho();
            h = hh.ptr();
        }
        else if (foundObject<volScalarField>("h"))
        {
            h = &lookupObject<volScalarField>("h");
        }
        else if (foundObject<volScalarField>("e"))
        {
            tmp<volScalarField> hh =
                lookupObject<volScalarField>("e") + p()/rho();
            h = hh.ptr();
        }
        else
        {
            FatalErrorInFunction
                << "Could not find sensible enthalpy or internal "
                << "energy field in database." << exit(FatalError);
        }
    }

    forAll(mesh_.boundary(), patchI)
    {
        if (objectivePatch_[patchI])
        {
            const fvPatch& patch = mesh_.boundary()[patchI];

            const scalarField rhop = rho()().boundaryField()[patchI];
            const vectorField& Vp = U().boundaryField()[patchI];

            scalarField sw( (Sfb[patchI] & devRhoReffb[patchI]) & Vp );

            scalarField pp = p().boundaryField()[patchI];
            scalarField phip = phi().boundaryField()[patchI];

            //grab outside pressure from resistive outlet
            if
            (
                isA<resistivePressureFvPatchScalarField>
                (
                    p().boundaryField()[patchI]
                )
            )
            {
                pp =
                    dynamic_cast<const resistivePressureFvPatchScalarField&>
                    (
                        p().boundaryField()[patchI]
                    ).p0();
            }

            //different formulation based on flow solver
            if (compressible())
            {
                if (h == nullptr)
                {
                    FatalErrorInFunction << "h not initialized" << endl;
                }
                //TODO: the formula based on enthalpy is wrong, don't use it
                const scalarField& hp = h->boundaryField()[patchI];

                //different treatment based on type of thermo used
                mass += sum(phip);

                if
                (
                    patch.type() == "inlet"
                 || patch.patch().physicalType() == "inlet"
                 || patch.type() == "outlet"
                 || patch.patch().physicalType() == "outlet"
                )
                {

                    scalar tmpPower = 0;

                    if (compressibleForm_ == "enthalpyPressure")
                    {
                        tmpPower =
                            sum
                            (
                                phip*(hp - pp/rhop)
                            );
                    }
                    else
                    {
                        tmpPower =
                            sum
                            (
                                phip*(pp/rhop + 0.5*magSqr(Vp))
                            );
                    }

                    if (stressWork_)
                    {
                        tmpPower += sum(sw);
                    }

                    powerLoss -= tmpPower;

                    if
                    (
                        patch.type() == "inlet"
                     || patch.patch().physicalType() == "inlet"
                    )
                    {
                        powerIn -= tmpPower;
                    }
                    else if
                    (
                        patch.type() == "outlet"
                     || patch.patch().physicalType() == "outlet"
                    )
                    {
                        powerOut -= tmpPower;
                    }
                }
            }
            else
            {
                if
                (
                    patch.type() == "inlet"
                 || patch.patch().physicalType() == "inlet"
                )
                {
                    scalar tmpPower = sum
                    (
                        rhop*phip*(pp + 0.5*magSqr(Vp))
                    );

                    if (stressWork_)
                    {
                        tmpPower += sum(sw);
                    }

                    powerIn -= tmpPower;
                    powerLoss -= tmpPower;
                }
                else if
                (
                    patch.type() == "outlet"
                 || patch.patch().physicalType() == "outlet"
                )
                {
                    scalar tmpPower = sum
                    (
                        rhop*phip*(pp + 0.5*magSqr(Vp))
                    );

                    if (stressWork_)
                    {
                        tmpPower += sum(sw);
                    }

                    powerOut -= tmpPower;
                    powerLoss -= tmpPower;
                }

                mass += sum(phip);
            }
        }
    }

    reduce(
        std::tie(powerLoss, mass, powerIn, powerOut),
        UniformParallelOp<sumOp<scalar>, 4>{}
    );

    objectiveValue_ = powerLoss;

    Info<< type() << " " << name() << " execute:" << nl
        << "Power losses = " << objectiveValue_ << " [W]" << nl << endl;

    return true;
}


bool
Foam::functionObjects::powerLossObjectiveFunctionObject::write()
{
    writeOutputToFile();

    return true;
}

// ************************************************************************* //
