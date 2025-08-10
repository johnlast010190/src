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
    (c) 2017-2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "disperseEulerian.H"
#include "solverOption/SolverOption.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "db/dictionary/dictionary.H"
#include "db/Time/Time.H"
#include "finiteVolume/fvc/fvc.H"
#include "fields/Fields/oneField/oneField.H"
#include "cfdTools/general/boundarySetup/boundarySetup.H"
#include "sets/topoSetSource/topoSetSource.H"
#include "cfdTools/general/bound/bound.H"
#include "db/IOobjectList/IOobjectList.H"
#include "interpolation/surfaceInterpolation/limitedSchemes/upwind/upwind.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(disperseEulerian, 0);
}
}

makeFvSolverOption(disperseEulerian);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::disperseEulerian::disperseEulerian
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict),

    solveTD_(false),
    alphaTD_(false),
    turb_(nullptr),
    stabDivAlpha_(0),
    UName_("U"),
    phiName_("phiwdr"),
    phases_(),
    forcesD_(),
    forcesL_(),
    forcesTD_(),
    sMass_(),
    gMass_(),
    stochasticDispersionModel_(false),
    phidPtr_(nullptr),
    timeIntOn_(false),
    startTime_(0),
    reportOn_(false),
    nNonOrthCorr_(0),
    fieldDependency_("none"),
    scalingFactor_(0)
{
    read(dict);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::disperseEulerian::~disperseEulerian()
{
    deleteDemandDrivenData(phidPtr_);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::disperseEulerian::initialise()
{
    // moved from constructor body
    if (debug && solveTD_)
    {
        Info<< "Solving with turbulent dispersion!" << endl;
    }

    turb_ = obr_.lookupObjectRefPtr<turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

    IOdictionary disperseEulerianProperties
    (
        IOobject
        (
            "disperseEulerianProperties",
            mesh_.time().constant(),
            obr_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    PtrList<entry> phaseEntries(disperseEulerianProperties.lookup("phases"));
    phases_.setSize(phaseEntries.size());
    forcesD_.setSize(phaseEntries.size());
    forcesL_.setSize(phaseEntries.size());
    forcesTD_.setSize(phaseEntries.size());
    sMass_.setSize(phaseEntries.size());
    gMass_.setSize(phaseEntries.size());

    forAll(phases_, phaseI)
    {
        phases_.set
        (
            phaseI,
            new decoupledEulerian::phase
            (
                phaseEntries[phaseI].keyword(),
                phaseEntries[phaseI].dict(),
                obr_,
                mesh_,
                true
            )
        );

        forcesD_.set
        (
            phaseI,
            decoupledEulerian::dragModel::New
            (
                phaseEntries[phaseI].dict(),
                phases_[phaseI],
                true
            )
        );

        forcesL_.set
        (
            phaseI,
            decoupledEulerian::liftModel::New
            (
                phaseEntries[phaseI].dict(),
                phases_[phaseI],
                true
            )
        );

        forcesTD_.set
        (
            phaseI,
            decoupledEulerian::turbulentDispersionModel::New
            (
                phaseEntries[phaseI].dict(),
                phases_[phaseI],
                true
            )
        );

        sMass_.set
        (
            phaseI,
            new volScalarField
            (
                IOobject
                (
                    "sMass."+Foam::word(phases_[phaseI].name()),
                    mesh_.time().timeName(),
                    obr_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE //AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar
                (
                    "sMass",
                    timeIntOn_ ? dimMass/dimArea : dimMass/(dimArea*dimTime),
                    0.0
                )
            )
        );

        gMass_.set
        (
            phaseI,
            new volScalarField
            (
                IOobject
                (
                    "gMass."+Foam::word(phases_[phaseI].name()),
                    mesh_.time().timeName(),
                    obr_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE //AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar
                (
                    "gMass",
                    timeIntOn_ ? dimMass : dimMass/dimTime,
                    0.0
                )
            )
        );

        //fluxRequired needed to compute report fluxed later
        mesh().schemes().setFluxRequired(phases_[phaseI].alphad().name());
    }

    phidPtr_ =
    (
        new surfaceScalarField
        (
            IOobject
            (
                "phid",
                mesh_.time().timeName(),
                obr_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::flux(Ud(0))
        )
    );

    return true;
}


void Foam::fv::disperseEulerian::read(const dictionary& dict)
{
    solveTD_ = dict.lookupOrDefault<Switch>("solveTD", false);
    alphaTD_ = dict.lookupOrDefault<Switch>("alphaTD", false);
    stabDivAlpha_ = dict.lookupOrDefault<scalar>("stabDivAlpha", 0.001);
    UName_ = dict.lookupOrDefault<word>("UName", "");
    phiName_ = dict.lookupOrDefault<word>("phiName", "phiwdr");
    stochasticDispersionModel_ =
        dict.lookupOrDefault<Switch>("stochasticDispersionModel", false);
    timeIntOn_ = dict.lookupOrDefault<Switch>("ReportSumOverTime", true);
    startTime_ = dict.lookupOrDefault<scalar>("startTime", 0.0);
    reportOn_ = dict.lookupOrDefault<Switch>("ReportOn", false);
    scalingFactor_ = dict.lookupOrDefault<scalar>("scalingFactor", scalar(1.0));

    nNonOrthCorr_ = dict.lookupOrDefault("nNonOrthogonalCorrectors", 0);
    fieldDependency_ = dict.lookupOrDefault<word>("fieldDependency", "none");
}


void Foam::fv::disperseEulerian::correct
(
    const word& solveName,
    const word& regionName
)
{
    if (debug)
    {
        Info<< "    " << "Solving for disperse Eulerian particles. " << endl;
    }

    const Time& runTime(mesh_.time());
    const fvMesh& mesh(mesh_);

    // using header to create g does not work as it tries to lookup in wrong place
    //#include "cfdTools/general/include/readGravitationalAcceleration.H"
    Info<< "\nReading g" << endl;
    uniformDimensionedVectorField g
    (
        IOobject
        (
            "g",
            mesh_.time().constant(),
            obr_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    // for all non-orth corrcetions
    for (int nonOrth=0; nonOrth <= nNonOrthCorr_; nonOrth++)
    {
        // for all spray fields
        for (label phase_no = 0; phase_no < phases_.size(); phase_no++)
        {
            volScalarField* nutd = nullptr;
            volScalarField* Ctd = nullptr;

            const dimensionedScalar& rhop = phases_[phase_no].rhod();
            //const dimensionedScalar& dp = phases_[phase_no].diam();

            volScalarField CdRe(forcesD_[phase_no].CdRe());

            // compute the flux
            Ud(phase_no).correctBoundaryConditions();

            if (phiName_ == "")
            {
                *phidPtr_ = turb_->phi() + fvc::flux(Ud(phase_no)-U());

                // flux upwind interpolation based on linear interpolated flux
                //surfaceScalarField phiLin = linearInterpolate(Ud(phase_no)) & mesh_.Sf();
                //*phidPtr_ = upwind<vector>(mesh_, phiLin).interpolate(Ud(phase_no)) & mesh_.Sf();

                // alpha-weighted flux interpolation
                //*phidPtr_ =
                //    (
                //        (linearInterpolate(alphad(phase_no) * Ud(phase_no))
                //      & mesh_.Sf())
                //      / linearInterpolate(alphad(phase_no) + VSMALL)
                //    );
            }
            else if (phiName_ == "phiwdr")
            {
                *phidPtr_ = linearInterpolate(Ud(phase_no)) & mesh.Sf();
            }
            else
            {
                const surfaceScalarField& phi = obr_.lookupObject<surfaceScalarField>(phiName_);
                *phidPtr_ = phi + fvc::flux(Ud(phase_no)-U());
            }

            #include "alphaCourantNo.H"
            #include "computePhysicalQuantities.H"
            #include "alphaEqns.H"
            #include "UEqns.H"

            delete nutd;
            delete Ctd;
        }
    }

    if (runTime.writeTime())
    {
        // write result - fields not written by default
        if (reportOn_)
        {
            for (int phase_no = 0; phase_no < phases_.size(); phase_no++)
            {
                sMass_[phase_no].write();
                gMass_[phase_no].write();
            }
        }
    }
}

// ************************************************************************* //
