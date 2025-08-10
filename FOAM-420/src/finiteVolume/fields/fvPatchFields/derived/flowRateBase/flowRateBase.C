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
    (c) 2017-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/flowRateBase/flowRateBase.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * Static Members * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(flowRateBase, 0);
}

const Foam::Enum<Foam::flowRateBase::massFlowUnits>
Foam::flowRateBase::massFlowUnitNames_
({
    {massFlowUnits::kg_s, "kg_s"},
    {massFlowUnits::kg_min, "kg_min"}
});

const Foam::Enum<Foam::flowRateBase::volumetricFlowUnits>
Foam::flowRateBase::volumetricFlowUnitNames_
({
    {volumetricFlowUnits::m3_s, "m3_s"},
    {volumetricFlowUnits::m3_hr, "m3_hr"},
    {volumetricFlowUnits::l_min, "l_min"},
    {volumetricFlowUnits::l_s, "l_s"}
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Null constructor
Foam::flowRateBase::flowRateBase
(
    const fvPatch& patch,
    const Foam::objectRegistry& db,
    bool targetType
)
:
    db_(db),
    patch_(patch),
    volumetric_(false),
    flowRate_(),
    rhoName_("rho"),
    rhoInlet_(0.0),
    rhoAverage_(true),
    massFlowUnits_(kg_s),
    volumetricFlowUnits_(m3_s),
    inflowPatchesSpecified_(false),
    curInletMassFlow_(0),
    phaseName_("none"),
    targetType_(targetType)
{}


//- Construct from components
Foam::flowRateBase::flowRateBase
(
    const fvPatch& patch,
    const Foam::objectRegistry& db,
    const dictionary& dict,
    const bool targetType
)
:
    db_(db),
    patch_(patch),
    rhoInlet_(dict.lookupOrDefault<scalar>("rhoInlet", -VGREAT)),
    massFlowUnits_(kg_s),
    volumetricFlowUnits_(m3_s),
    inflowPatchesSpecified_(false),
    curInletMassFlow_(dict.lookupOrDefault<scalar>("curInletMassFlow", 0.0)),
    phaseName_(dict.lookupOrDefault<word>("phaseName", "none")),
    targetType_(targetType)
{
    word massFlowSplitKeyword = flowRateKeyword("massFlowSplit");
    if
    (
        dict.found(flowRateKeyword("massFlowRate"))
     || dict.found(massFlowSplitKeyword)
    )
    {
        volumetric_ = false;
        if (dict.found(massFlowSplitKeyword))
        {
            flowRate_ = Function1<scalar>::New(massFlowSplitKeyword, dict);
            massFlowUnits_ = kg_s;

            inflowPatchSet_.clear();
            if (dict.found("patches"))
            {
                inflowPatchSet_ =
                    patch_.patch().boundaryMesh().patchSet
                    (
                        wordReList(dict.lookup("patches")), true, true
                    );
                inflowPatchesSpecified_ = true;
            }
        }
        else
        {
            flowRate_ =
                Function1<scalar>::New
                (
                    flowRateKeyword("massFlowRate"), dict
                );

            if (dict.found("units"))
            {
                massFlowUnits_ = massFlowUnitNames_.lookup("units", dict);
            }
            else
            {
                massFlowUnits_ = massFlowUnits::kg_s;
            }
        }
        rhoName_ = dict.lookupOrDefault<word>("rho", "rho");
        rhoAverage_
            = Switch(dict.lookupOrDefault<Switch>("rhoAverage", true));
    }
    else if
    (
        dict.found(flowRateKeyword("volumetricFlowRate"))
     || dict.found(flowRateKeyword("flowRate"))
    )
    {
        volumetric_ = true;
        if (dict.found(flowRateKeyword("flowRate")))
        {
            flowRate_ =
                Function1<scalar>::New
                (
                    flowRateKeyword("flowRate"), dict
                );
        }
        else
        {
            flowRate_ =
                Function1<scalar>::New
                (
                    flowRateKeyword("volumetricFlowRate"), dict
                );
        }
        rhoName_ = "rho";
        if (dict.found("lpmUnits"))
        {
            DeprecationWarningInFunction("lpmUnits", "keyword", 40200)
                << "Please specify 'units' instead." << nl << endl;
            volumetricFlowUnits_ = volumetricFlowUnits::l_min;
        }
        else
        {
            volumetricFlowUnits_ = volumetricFlowUnits::m3_s;
        }

        if (dict.found("units"))
        {
            volumetricFlowUnits_ =
                volumetricFlowUnitNames_.lookup("units", dict);
        }
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Please supply either 'volumetricFlowRate' or"
            << " 'massFlowRate'." << nl
            << exit(FatalIOError);
    }
}


//- Copy construct
Foam::flowRateBase::flowRateBase
(
    const fvPatch& patch,
    const objectRegistry& db,
    const flowRateBase& frb,
    const bool targetType
)
:
    db_(db),
    patch_(patch),
    volumetric_(frb.volumetric_),
    flowRate_(frb.flowRate_, false),
    rhoName_(frb.rhoName_),
    rhoInlet_(frb.rhoInlet_),
    rhoAverage_(frb.rhoAverage_),
    massFlowUnits_(frb.massFlowUnits_),
    volumetricFlowUnits_(frb.volumetricFlowUnits_),
    inflowPatchSet_(frb.inflowPatchSet_),
    inflowPatchesSpecified_(frb.inflowPatchesSpecified_),
    curInletMassFlow_(frb.curInletMassFlow_),
    phaseName_(frb.phaseName_),
    targetType_(targetType)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::flowRateBase::~flowRateBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::flowRateBase::currentFlowRate() const
{
    scalar fr = flowRate_->value(db_.time().timeOutputValue());

    if (isSplitBC())
    {
        fr *= curInletMassFlow_;
    }
    else
    {
        if (volumetric())
        {
            switch(volumetricFlowUnits_)
            {
                case m3_s:
                {
                    break;
                }
                case m3_hr:
                {
                    fr /= scalar(3600);
                    break;
                }
                case l_min:
                {
                    fr /= scalar(60000);
                    break;
                }
                case l_s:
                {
                    fr *= scalar(0.001);
                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Missing support for volumetric units: "
                        << volumetricFlowUnitNames_[volumetricFlowUnits_]
                        << exit(FatalError);
                }
            }
        }
        else
        {
            switch(massFlowUnits_)
            {
                case kg_s:
                {
                    break;
                }
                case kg_min:
                {
                    fr /= 60;
                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Missing support for mass flow units: "
                        << massFlowUnitNames_[massFlowUnits_]
                        << exit(FatalError);
                }
            }
        }
    }

    // update target flow rate based on wetted area
    if (phaseName_ != "none")
    {
        IOdictionary transProp
        (
            IOobject
            (
                "transportProperties",
                patch_.boundaryMesh().mesh().time().constant(),
                patch_.boundaryMesh().mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        word phaseName(phaseName_);
        // Accept 'alpha.<phaseName>' as well
        if (phaseName.substr(0, 6) == "alpha.")
        {
            phaseName = phaseName.substr(6, phaseName.size());
        }

        dimensionedScalar rhoPhase(transProp.subDict(phaseName).lookup("rho"));
        scalar rhoMean = averageRhoOrOne();

        if (debug)
        {
            Info<< "scaling multiphase massflow rate by : "
                << rhoMean/rhoPhase.value() << nl
                << "  target flow rate unscaled : " << fr << nl
                << "  target flow rate scaled : " << fr*rhoMean/rhoPhase.value() << endl;
        }

        fr *= rhoMean/rhoPhase.value();
    }

    return fr;
}


Foam::scalar Foam::flowRateBase::averageRhoOrOne() const
{
    scalar areaWeightedDensityAverage =
        gSum(multiplyByRhoOrOne(patch_.magSf(), false))
      / gSum(patch_.magSf());
    return areaWeightedDensityAverage;
}


Foam::tmp<Foam::scalarField> Foam::flowRateBase::multiplyByRhoOrOne
(
    const scalarField& sF, bool useRhoAverage
) const
{
    if (volumetric_ || rhoName() == "none")
    {
        return tmp<scalarField>(new scalarField(sF));
    }
    else
    {
        if (db_.foundObject<volScalarField>(rhoName()))
        {
            const fvPatchField<scalar>& rhoPatch =
                patch_.lookupPatchFieldInDb<volScalarField, scalar>
                (
                    this->db_,
                    rhoName()
                );

            if (useRhoAverage && rhoAverage_)
            {
                scalar rhoAve = gSum(rhoPatch*patch_.magSf())
                                    / gSum(patch_.magSf());

                return tmp<scalarField>(new scalarField(sF*rhoAve));
            }
            else
            {
                return tmp<scalarField>(new scalarField(sF*rhoPatch));
            }
        }
        else
        {
            // Use constant density
            if (rhoInlet() < 0)
            {
                FatalErrorInFunction
                    << "Did not find registered density field "
                    << rhoName()
                    << " and no constant density 'rhoInlet' specified"
                    << exit(FatalError);
            }

            return tmp<scalarField>(new scalarField(sF*rhoInlet_));
        }
    }
}


void Foam::flowRateBase::updateInletFlowRate(const List<scalar>& inletFluxes)
{
    // Do sum based on selected patches if specified
    curInletMassFlow_ = scalar(0);
    forAll(inletFluxes, patchi)
    {
        // If no patches specified, include them all
        if (!inflowPatchesSpecified_ || inflowPatchSet_.found(patchi))
        {
            curInletMassFlow_ += inletFluxes[patchi];
        }
    }
}


void Foam::flowRateBase::write(Ostream& os) const
{
    flowRate_->writeData(os);

    if (volumetric_)
    {
        if (!isSplitBC())
        {
            os.writeEntry
            (
                "units",
                word(volumetricFlowUnitNames_[volumetricFlowUnits_])
            );
        }
    }
    else
    {
        if (!isSplitBC())
        {
            os.writeEntry("units", word(massFlowUnitNames_[massFlowUnits_]));
        }

        fvPatchVectorField::writeEntryIfDifferent<word>
        (
            os, "rho", "rho", rhoName_
        );
        fvPatchVectorField::writeEntryIfDifferent<scalar>
        (
            os, "rhoInlet", -VGREAT, rhoInlet_
        );
        fvPatchVectorField::writeEntryIfDifferent<Switch>
        (
            os, "rhoAverage", true, rhoAverage_
        );
    }

    if (isSplitBC())
    {
        if (inflowPatchesSpecified_)
        {
            DynamicList<word> patchList;
            forAllConstIters(inflowPatchSet_, patchi)
            {
                patchList.append(patch_.boundaryMesh()[patchi()].name());
            }
            os.writeEntry("patches", patchList);
        }
        fvPatchVectorField::writeEntryIfDifferent<scalar>
        (
            os, "curInletMassFlow", 0.0, curInletMassFlow_
        );
    }

    fvPatchVectorField::writeEntryIfDifferent<word>
    (
        os, "phaseName", "none", phaseName_
    );
}


// ************************************************************************* //