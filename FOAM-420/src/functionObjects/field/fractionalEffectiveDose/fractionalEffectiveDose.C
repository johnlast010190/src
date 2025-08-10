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
    (c) 2022-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fractionalEffectiveDose/fractionalEffectiveDose.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fractionalEffectiveDose, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        fractionalEffectiveDose,
        dictionary
    );
}
}

Foam::wordList Foam::functionObjects::fractionalEffectiveDose::
speciesFields = {"CO", "CO2", "O2"};

const Foam::Enum
<
    Foam::functionObjects::fractionalEffectiveDose::activityLevel
>
Foam::functionObjects::fractionalEffectiveDose::activityLevelNames
{
    { activityLevel::activityRest, "atRest"},
    { activityLevel::activityLight, "lightWork"},
    { activityLevel::activityHeavy, "heavyWork"},
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::fractionalEffectiveDose::calc()
{
    if (!initialised_)
    {
        // Ensure first calculation works unconditionally
        prevTimeIndex_ = -1;

        initialised_ = true;
    }

    const label currentTimeIndex = obr().time().timeIndex();
    const scalar currentTime = obr().time().value();

    if (prevTimeIndex_ == currentTimeIndex)
    {
        return true;
    }
    else
    {
        prevTimeIndex_ = currentTimeIndex;
    }

    bool doReset = false;
    if (periodicReset_ && currentTime >= resetPeriod_*periodIndex_)
    {
        doReset = true;
        periodIndex_++;
    }

    if (currentTime >= resetTime_ && resetTime_ > SMALL)
    {
        doReset = true;        // Reset is overdue
        resetTime_ = GREAT;    // Avoid triggering again
    }

    if (doReset)
    {
        reset();
    }

    // Main driver function
    calculateFED();

    // Reset previous time
    prevTime_ = currentTime;

    return true;
}


void Foam::functionObjects::fractionalEffectiveDose::calculateFED()
{
    const volScalarField& alphaCO
    (
        lookupObject<volScalarField>(speciesFields[0])
    );

    const volScalarField& alphaCO2
    (
        lookupObject<volScalarField>(speciesFields[1])
    );

    const volScalarField& alphaO2
    (
        lookupObject<volScalarField>(speciesFields[2])
    );

    // Species concentrations
    const volScalarField C_CO(alphaCO*1e6);      // [ppm]
    const volScalarField C_CO2(alphaCO2*100);    // [%]
    const volScalarField C_O2(alphaO2*100);      // [%]

    const scalar dt = (obr().time().value() - prevTime_)/60;    // [min]

    // Calculate FED field
    volScalarField& fed =
        const_cast<volScalarField&>
        (
            lookupObject<volScalarField>(resultName_)
        );

    forAll(fed, celli)
    {
        fedCO_[celli] += FEDcoeff_*Foam::pow(C_CO[celli], 1.036)*dt;
        fedO2_[celli] += dt/(Foam::exp(8.13 - 0.51*(20.9 - C_O2[celli])));
        scalar HVCO2 = Foam::exp(0.1903*C_CO2[celli] + 2.0004)/7.1;

        fed[celli] = fedCO_[celli]*HVCO2 + fedO2_[celli];
    }

    fed.correctBoundaryConditions();

    volScalarField& fedCO =
        const_cast<volScalarField&>
        (
            lookupObject<volScalarField>(FEDFieldNames_[0])
        );

    volScalarField& fedO2 =
        const_cast<volScalarField&>
        (
            lookupObject<volScalarField>(FEDFieldNames_[1])
        );

    forAll(fed, celli)
    {
        fedCO[celli] = fedCO_[celli];
        fedO2[celli] = fedO2_[celli];
    }

    fedCO.correctBoundaryConditions();
    fedO2.correctBoundaryConditions();

    // Get averaged quantities
    updateAverages();
}


void Foam::functionObjects::fractionalEffectiveDose::updateAverages()
{
    // Retrieve species phase-fractions (optionally filtered for cellZones)
    scalarField alphaCO
    (
        filterField(lookupObject<volScalarField>(speciesFields[0]))
    );

    scalarField alphaCO2
    (
        filterField(lookupObject<volScalarField>(speciesFields[1]))
    );

    scalarField alphaO2
    (
        filterField(lookupObject<volScalarField>(speciesFields[2]))
    );

    // Cell volumes
    scalarField VField(filterField<scalar>(fvMeshFunctionObject::mesh_.V()));
    scalar V(gSum(VField));

    // Volume-averaged species concentrations
    const scalar C_CO(gSum((alphaCO*1e6)*VField)/V);      // [ppm]
    const scalar C_CO2(gSum((alphaCO2*100)*VField)/V);    // [%]
    const scalar C_O2(gSum((alphaO2*100)*VField)/V);      // [%]

    const scalar dt = (obr().time().value() - prevTime_)/60;    // [min]

    // Expressions for relevant averaged quantities
    avgFedCO_.value() += FEDcoeff_*Foam::pow(C_CO, 1.036)*dt;
    avgFedO2_.value() += dt/(Foam::exp(8.13 - 0.51*(20.9 - C_O2)));
    avgHVCO2_.value() = Foam::exp(0.1903*C_CO2 + 2.0004)/7.1;

    // Update average FED of incapacitation
    avgFedIN_ = avgFedCO_*avgHVCO2_ + avgFedO2_;
}


void Foam::functionObjects::fractionalEffectiveDose::reset()
{
    Info<< "    Resetting FED calculation at time "
        << obr().time().timeOutputValue() << nl << endl;

    volScalarField& fed =
        const_cast<volScalarField&>
        (
            lookupObject<volScalarField>(resultName_)
        );

    fed.primitiveFieldRef() = Zero;
    fedCO_.primitiveFieldRef() = Zero;
    fedO2_.primitiveFieldRef() = Zero;

    avgFedIN_ = avgFedCO_ = avgFedO2_ = avgHVCO2_ = 0;

    volScalarField& fedCO =
        const_cast<volScalarField&>
        (
            lookupObject<volScalarField>(FEDFieldNames_[0])
        );

    fedCO.primitiveFieldRef() = Zero;

    volScalarField& fedO2 =
        const_cast<volScalarField&>
        (
            lookupObject<volScalarField>(FEDFieldNames_[1])
        );

    fedO2.primitiveFieldRef() = Zero;

    // Ensure first calculation works unconditionally
    prevTimeIndex_ = -1;

    initialised_ = true;
}


template<typename Type>
Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::fractionalEffectiveDose::filterField
(
    const Field<Type>& field
) const
{
    if (isNull(cellIDs()))
    {
        return field;
    }
    else
    {
        return tmp<Field<Type>>(new Field<Type>(field, cellIDs()));
    }
}


void Foam::functionObjects::fractionalEffectiveDose::writeFileHeader
(
    Ostream& os
) const
{
    writeHeader(os, "Fractional effective dose (average)");
    writeCommented(os, "Time");
    writeDelimited(os, "FED (incapacitation)");
    writeDelimited(os, "FED (CO)");
    writeDelimited(os, "HV_CO2");
    writeDelimited(os, "FED (O2)");
    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fractionalEffectiveDose::fractionalEffectiveDose
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict),
    volRegion(fvMeshFunctionObject::mesh_, dict),
    writeFile(obr_, name, typeName, dict),
    avgFedIN_("avgFedIN", dimTime, Zero),
    avgFedCO_("avgFedCO", dimTime, Zero),
    avgFedO2_("avgFedO2", dimTime, Zero),
    avgHVCO2_("avgHVCO2", dimless, Zero),
    fedCO_
    (
        IOobject
        (
            "fedCO",
            time_.timeName(),
            fvMeshFunctionObject::mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        fvMeshFunctionObject::mesh_,
        dimensionedScalar("zero", dimTime, Zero)
    ),
    fedO2_
    (
        IOobject
        (
            "fedO2",
            time_.timeName(),
            fvMeshFunctionObject::mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        fvMeshFunctionObject::mesh_,
        dimensionedScalar("zero", dimTime, Zero)
    ),
    activityLevel_(activityLevel::activityLight),
    FEDcoeff_(0),
    prevTime_(0),
    prevTimeIndex_(-1),
    periodicReset_(false),
    periodIndex_(1),
    resetPeriod_(GREAT),
    resetTime_(GREAT),
    initialised_(false),
    FEDFieldNames_()
{
    read(dict);
    writeFileHeader(file());
    setResultName("fed", "fed");

    tmp<volScalarField> fedFieldPtr
    (
        new volScalarField
        (
            IOobject
            (
                resultName_,
                time_.timeName(),
                fvMeshFunctionObject::mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvMeshFunctionObject::mesh_,
            dimensionedScalar("zero", dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    store(resultName_, fedFieldPtr);

    FEDFieldNames_ = {"fedCO", "fedO2"};

    // Create and store fields for FED_CO and FED_O2
    tmp<volScalarField> fedCOFieldPtr(new volScalarField(fedCO_));
    tmp<volScalarField> fedO2FieldPtr(new volScalarField(fedO2_));

    store(FEDFieldNames_[0], fedCOFieldPtr);
    store(FEDFieldNames_[1], fedO2FieldPtr);

    if (obr().time().timeOutputValue() > 0 && avgFedCO_.value() < SMALL)
    {
        // Calculate averages for a restart

        // Cell volumes
        scalarField VField
        (
            filterField<scalar>(fvMeshFunctionObject::mesh_.V())
        );
        scalar V(gSum(VField));

        // Expressions for relevant averaged quantities
        avgFedCO_.value() = gSum(fedCO_*VField)/V;
        avgFedO2_.value() = gSum(fedO2_*VField)/V;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fractionalEffectiveDose::~fractionalEffectiveDose()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fractionalEffectiveDose::read
(
    const dictionary& dict
)
{
    if (fieldExpression::read(dict) && writeFile::read(dict))
    {
        // Determine FED coefficients according to the activity level
        if (dict.found("activityLevel"))
        {
            activityLevel_ = activityLevelNames.read
            (
                dict.lookup("activityLevel")
            );
        }

        switch (activityLevel_)
        {
            case activityLevel::activityRest:
            {
                FEDcoeff_ = 7.04863E-6;
                break;
            }
            case activityLevel::activityLight:
            {
                FEDcoeff_ = 2.76417E-5;
                break;
            }
            case activityLevel::activityHeavy:
            {
                FEDcoeff_ = 8.29250E-6;
                break;
            }
            default:
            {
                FatalErrorInFunction
                    << "    Unknown activity level. Valid types are:"
                    << activityLevelNames << nl << exit(FatalError);
            }
        }

        // Reset control parameters to their defaults
        initialised_ = false;
        periodicReset_ = false;
        resetPeriod_ = GREAT;
        resetTime_ = GREAT;

        Info<< type() << " " << name() << ":" << nl;

        dict.readIfPresent("periodicReset", periodicReset_);

        const scalar currentTime = obr().time().value();

        if (periodicReset_)
        {
            scalar userResetPeriod = readScalar(dict.lookup("resetPeriod"));
            resetPeriod_ = obr().time().userTimeToTime(userResetPeriod);

            if (resetPeriod_ > 0)
            {
                // Determine the appropriate interval for the next reset
                periodIndex_ = 1;
                while (currentTime > resetPeriod_*periodIndex_)
                {
                    ++periodIndex_;
                }

                Info<< "    Resetting period " << userResetPeriod
                    << ": next reset at " << (userResetPeriod*periodIndex_)
                    << nl << endl;
            }
            else
            {
                periodicReset_ = false;

                Info<< "    Resetting period " << userResetPeriod
                    << ": ignored" << nl << endl;
            }
        }

        scalar userResetTime = 0;
        if (dict.readIfPresent("resetTime", userResetTime))
        {
            resetTime_ = obr().time().userTimeToTime(userResetTime);

            if (resetTime_ > SMALL)
            {
                if (currentTime > resetTime_)
                {
                    // The reset time is already in the past: ignore
                    resetTime_ = GREAT;
                }
                else
                {
                    Info<< "    Resetting scheduled at time " << userResetTime
                        << nl << endl;
                }
            }
            else
            {
                Info<< "    Reset time " << userResetTime
                    << ": ignored" << nl << endl;
            }
        }

        Info<< "    Starting FED calculation at time "
            << obr().time().timeOutputValue() << endl;

        if (regionType_ == vrtCellZone)
        {
            Info<< "    Average quantities calculated for cellZone "
            << regionName_ << nl;
        }

        // Set initial time
        prevTime_ = currentTime;

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::functionObjects::fractionalEffectiveDose::write()
{
    if (fieldExpression::write())
    {
        // Output to file
        writeTime(file());

        file()
            << token::TAB << avgFedIN_.value()
            << token::TAB << avgFedCO_.value()
            << token::TAB << avgHVCO2_.value()
            << token::TAB << avgFedO2_.value() << endl;
    }

    writeObject(FEDFieldNames_[0]);
    writeObject(FEDFieldNames_[1]);

    return true;
}

// ************************************************************************* //
