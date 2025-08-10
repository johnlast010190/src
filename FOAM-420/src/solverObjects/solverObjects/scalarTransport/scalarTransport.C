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
    (c) 2019 Esi Ltd.
    (c) 2012-2017 OpenFOAM Foundation
    (c) 2015-2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "scalarTransport.H"
#include "solverOption/SolverOption.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "finiteVolume/fvm/fvmDdt.H"
#include "finiteVolume/fvm/fvmDiv.H"
#include "finiteVolume/fvm/fvmLaplacian.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "fvMatrices/solvers/MULES/CMULES.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "algorithms/subCycle/subCycle.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(scalarTransport, 0);
}
}

makeFvSolverOption(scalarTransport);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::volScalarField& Foam::fv::scalarTransport::transportedField()
{
    if (!mesh_.foundObject<volScalarField>(fieldName_))
    {
        tmp<volScalarField> tfldPtr
        (
            new volScalarField
            (
                IOobject
                (
                    fieldName_,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
        store(fieldName_, tfldPtr);

        if (phaseName_ != "none")
        {
            mesh_.schemes().setFluxRequired(fieldName_);
        }
    }

    return const_cast<volScalarField&>
    (
        mesh_.lookupObject<volScalarField>(fieldName_)
    );
}


Foam::tmp<Foam::volScalarField> Foam::fv::scalarTransport::D
(
    const volScalarField& s,
    const surfaceScalarField& phi
) const
{
    typedef incompressible::turbulenceModel icoModel;
    typedef compressible::turbulenceModel cmpModel;

    word Dname("D" + s.name());

    if (constantD_)
    {
        tmp<volScalarField> tD
        (
            new volScalarField
            (
                IOobject
                (
                    Dname,
                    mesh_.time().timeName(),
                    mesh_.time(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(Dname, phi.dimensions()/dimLength, D_)
            )
        );

        return tD;
    }
    else if (nutName_ != "none")
    {
        const volScalarField& nutMean =
            mesh_.lookupObject<volScalarField>(nutName_);

        return tmp<volScalarField>
        (
            new volScalarField(Dname, nutMean)
        );
    }
    else if (mesh_.foundObject<icoModel>(turbulenceModel::propertiesName))
    {
        const icoModel& model = mesh_.lookupObject<icoModel>
        (
            turbulenceModel::propertiesName
        );

        return tmp<volScalarField>
        (
             new volScalarField
             (
                 Dname,
                 alphaD_*model.nu() + alphaDt_*model.nut()
             )
        );
    }
    else if (mesh_.foundObject<cmpModel>(turbulenceModel::propertiesName))
    {
        const cmpModel& model = mesh_.lookupObject<cmpModel>
        (
            turbulenceModel::propertiesName
        );

        return tmp<volScalarField>
        (
             new volScalarField
             (
                 Dname,
                 alphaD_*model.mu() + alphaDt_*model.mut()
             )
        );
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    Dname,
                    mesh_.time().timeName(),
                    mesh_.time(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(Dname, phi.dimensions()/dimLength, 0.0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::scalarTransport::scalarTransport
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict),
    fieldName_(dict.lookupOrDefault<word>("field", "s")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    nutName_(dict.lookupOrDefault<word>("nut", "none")),
    phaseName_(dict.lookupOrDefault<word>("phase", "none")),
    phasePhiCompressedName_
    (
        dict.lookupOrDefault<word>("phasePhiCompressed", "alphaPhiUn")
    ),
    D_(0),
    constantD_(false),
    nCorr_(0),
    resetOnStartUp_(false),
    schemesField_("unknown-schemesField"),
    bounded01_(dict.lookupOrDefault<Switch>("bounded01", true)),
    solveSpeedup_(dict.lookupOrDefault<label>("solveSpeedup", 1)),
    nSubCycles_(dict.lookupOrDefault<label>("nSubCycles", solveSpeedup_)),
    residualPhaseField_
    (
        "residualPhaseField",
        dimless,
        dict.lookupOrDefault<scalar>("residualPhaseField", 1e-5)
    ),
    fieldDependency_(dict.lookupOrDefault<word>("fieldDependency", "none")),
    transient_(dict.lookupOrDefault<Switch>("transient", true)),
    bounded_(dict.lookupOrDefault<Switch>("bounded", false)),
    zeroWallDiffusion_(dict.lookupOrDefault<Switch>("zeroWallDiffusion", true)),
    stablePhasic_(dict.lookupOrDefault<Switch>("stablePhasic", false)),
    phasicSources_(dict.lookupOrDefault<Switch>("phasicSources", true))
{
    read(dict);

    // Force creation of transported field so any BCs using it can
    // look it up
    volScalarField& s = transportedField();

    if (resetOnStartUp_)
    {
        s.forceAssign(dimensionedScalar("zero", dimless, 0.0));
    }
    if (stablePhasic_)
    {
        residualPhaseField_ = dimensionedScalar("zero", dimless, 0.0);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::scalarTransport::~scalarTransport()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::scalarTransport::read(const dictionary& dict)
{
    dict.readIfPresent("phi", phiName_);
    dict.readIfPresent("rho", rhoName_);
    dict.readIfPresent("nut", nutName_);
    dict.readIfPresent("phase", phaseName_);
    dict.readIfPresent("bounded01", bounded01_);
    dict.readIfPresent("solveSpeedup", solveSpeedup_);

    schemesField_ = dict.lookupOrDefault("schemesField", fieldName_);
    constantD_ = dict.readIfPresent("D", D_);
    alphaD_ = dict.lookupOrDefault("alphaD", 1.0);
    alphaDt_ = dict.lookupOrDefault("alphaDt", 1.0);

    dict.readIfPresent("nCorr", nCorr_);
    dict.readIfPresent("resetOnStartUp", resetOnStartUp_);

    dict.readIfPresent("fieldDependency", fieldDependency_);
    dict.readIfPresent("transient", transient_);
    dict.readIfPresent("bounded", bounded_);
    dict.readIfPresent("zeroWallDiffusion", zeroWallDiffusion_);
    dict.readIfPresent("stablePhasic", stablePhasic_);
    dict.readIfPresent("phasicSources", phasicSources_);
}


void Foam::fv::scalarTransport::correct
(
    const word& solveName,
    const word& regionName
)
{
    fv::options& fvOptions = this->fvOptions();

    volScalarField& s = transportedField();

    //Log << type() << " execute: " << s.name() << endl;

    const surfaceScalarField& phi =
        mesh_.lookupObject<surfaceScalarField>(phiName_);

    // Calculate the diffusivity
    volScalarField D(this->D(s, phi));

    if (zeroWallDiffusion_)
    {
        // set wall diffusion to zero for conservation
        forAll(D.boundaryField(), pI)
        {
            if (!D.boundaryField()[pI].coupled())
            {
                D.boundaryFieldRef()[pI] = 0.0;
            }
        }
    }

    // modify deltaT
    Time& time = const_cast<Time&>(mesh_.time());
    scalar deltaT = time.deltaTValue();

    if (nSubCycles_ > 1)
    {
        time.setDeltaT(deltaT*solveSpeedup_);

        for
        (
            subCycle<volScalarField> scalarSubCycle(s, nSubCycles_);
            !(++scalarSubCycle).end();
        )
        {
            solveEquation
            (
                s,
                phi,
                D,
                fvOptions
            );
        }

        // reset deltaT
        time.setDeltaT(deltaT);
    }
    else
    {
        solveEquation
        (
            s,
            phi,
            D,
            fvOptions
        );
    }

    //Log << endl;
}


void Foam::fv::scalarTransport::solveEquation
(
    volScalarField& s,
    const surfaceScalarField& phi,
    volScalarField& D,
    fv::options& fvOptions
)
{
    word divScheme("div(phi," + schemesField_ + ")");
    word laplacianScheme("laplacian(" + D.name() + "," + schemesField_ + ")");

    // Set under-relaxation coeff
    scalar relaxCoeff = 0.0;
    if (mesh_.solution().relaxEquation(schemesField_))
    {
        relaxCoeff = mesh_.solution().equationRelaxationFactor(schemesField_);
    }

    // Two phase scalar transport
    if (phaseName_ != "none")
    {
        const volScalarField& alpha =
            mesh_.lookupObject<volScalarField>(phaseName_);

        const surfaceScalarField& limitedPhiAlpha =
            mesh_.lookupObject<surfaceScalarField>(phasePhiCompressedName_);

        if (stablePhasic_)
        {
            D *= pos(alpha - 0.99);
        }
        else
        {
            D *= alpha;
        }

        // Reset D dimensions consistent with limitedPhiAlpha
        D.dimensions().reset(limitedPhiAlpha.dimensions()/dimLength);

        // multiplicator for fvOption sources
        volScalarField phaseField
        (
            IOobject
            (
                "phase",
                mesh_.time().timeName(),
                obr_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("phase", dimless, 1.0)
        );
        if (phasicSources_)
        {
            phaseField.forceAssign(alpha);
        }

        // Solve
        tmp<surfaceScalarField> tTPhiUD;
        for (label i = 0; i <= nCorr_; i++)
        {
            fvScalarMatrix sEqn
            (
                fvm::div(limitedPhiAlpha, s, divScheme)
              - fvm::laplacian(D, s, laplacianScheme)
             ==
                phaseField*fvOptions(s)
            );


            if (transient_)
            {
                if (stablePhasic_)
                {
                    sEqn += fvm::ddt(s);
                }
                else
                {
                    sEqn += fvm::ddt(alpha, s);

                    if (bounded_)
                    {
                        // use implicit correction for mass continuity error
                        volScalarField continuityError
                        (
                            fvc::ddt(alpha) + fvc::div(limitedPhiAlpha)
                        );

                        sEqn -= fvm::Sp(continuityError, s);
                    }
                }
            }
            else
            {
                if (bounded_)
                {
                    // use implicit "bounded" scheme
                    sEqn -= fvm::Sp(fvc::div(limitedPhiAlpha), s);
                }
            }

            if (residualPhaseField_.value() > 0)
            {
                sEqn += fvm::ddt(residualPhaseField_, s)
                      - fvc::ddt(residualPhaseField_, s);
            }

            sEqn.relax(relaxCoeff);
            fvOptions.constrain(sEqn);
            sEqn.solve(mesh_.solution().solverDict(schemesField_));

            tTPhiUD = sEqn.flux();
        }

        if (bounded01_)
        {
            MULES::explicitSolve(s, phi, tTPhiUD.ref(), 1, 0);
        }
    }
    else if (phi.dimensions() == dimMass/dimTime)
    {
        const volScalarField& rho = mesh_.lookupObject<volScalarField>(rhoName_);

        for (label i = 0; i <= nCorr_; i++)
        {

            fvScalarMatrix sEqn
            (
                fvm::div(phi, s, divScheme)
              - fvm::laplacian(D, s, laplacianScheme)
             ==
                fvOptions(rho, s)
            );

            if (transient_)
            {
                sEqn += fvm::ddt(rho, s);

                if (bounded_)
                {
                    // use implicit correction for mass continuity error
                    volScalarField continuityError(fvc::ddt(rho) + fvc::div(phi));

                    sEqn -= fvm::Sp(continuityError, s);
                }
            }
            else
            {
                if (bounded_)
                {
                    // use implicit "bounded" scheme
                    sEqn -= fvm::Sp(fvc::div(phi), s);
                }
            }

            sEqn.relax(relaxCoeff);

            fvOptions.constrain(sEqn);

            sEqn.solve(mesh_.solution().solverDict(schemesField_));
        }
    }
    else if (phi.dimensions() == dimVolume/dimTime)
    {
        for (label i = 0; i <= nCorr_; i++)
        {
            fvScalarMatrix sEqn
            (
                fvm::div(phi, s, divScheme)
              - fvm::laplacian(D, s, laplacianScheme)
             ==
                fvOptions(s)
            );

            if (transient_)
            {
                sEqn += fvm::ddt(s);

                if (bounded_)
                {
                    // use implicit correction for mass continuity error
                    volScalarField continuityError = fvc::div(phi)();

                    sEqn -= fvm::Sp(continuityError, s);
                }
            }
            else
            {
                if (bounded_)
                {
                    // use implicit "bounded" scheme
                    sEqn -= fvm::Sp(fvc::div(phi), s);
                }
            }

            sEqn.relax(relaxCoeff);

            fvOptions.constrain(sEqn);

            sEqn.solve(mesh_.solution().solverDict(schemesField_));
        }
    }
    else
    {
        FatalErrorInFunction
            << "Incompatible dimensions for phi: " << phi.dimensions() << nl
            << "Dimensions should be " << dimMass/dimTime << " or "
            << dimVolume/dimTime << exit(FatalError);
    }
}


void Foam::fv::scalarTransport::write()
{}



template<class ObjectType>
bool Foam::fv::scalarTransport::store
(
    word& fieldName,
    const tmp<ObjectType>& tfield,
    bool cacheable
)
{
    if (cacheable && fieldName == tfield().name())
    {
        WarningInFunction
            << "Cannot store cache-able field with the name used in the cache."
            << nl
            << "    Either choose a different name or cache the field"
            << "    and use the 'writeObjects' functionObject."
            << endl;

        return false;
    }

    if (fieldName.size() && mesh_.foundObject<ObjectType>(fieldName))
    {
        const ObjectType& field = mesh_.lookupObject<ObjectType>(fieldName);

        // If there is a result field already registered, assign to the new
        // result field. Otherwise transfer ownership of the new result field to
        // the object registry
        if (&field != &tfield())
        {
            const_cast<ObjectType&>(field) = tfield;
        }
        else
        {
            mesh_.thisDb().objectRegistry::store(tfield.ptr());
        }
    }
    else
    {
        if (fieldName.size() && fieldName != tfield().name())
        {
            tfield.ref().rename(fieldName);
        }
        else
        {
            fieldName = tfield().name();
        }

        mesh_.thisDb().objectRegistry::store(tfield.ptr());
    }

    return true;
}

// ************************************************************************* //
