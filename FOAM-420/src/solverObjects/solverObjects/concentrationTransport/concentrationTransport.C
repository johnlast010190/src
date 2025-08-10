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
    (c) 2011-2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "concentrationTransport.H"
#include "solverOption/SolverOption.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "db/dictionary/dictionary.H"
#include "db/Time/Time.H"
#include "finiteVolume/fvc/fvc.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "fields/Fields/oneField/oneField.H"
#include "cfdTools/general/boundarySetup/boundarySetup.H"
#include "sets/topoSetSource/topoSetSource.H"
#include "cfdTools/general/bound/bound.H"
#include "derivedFvPatchFields/phaseChangeHumidity/phaseChangeHumidityFvPatchScalarField.H"
#include "algorithms/subCycle/subCycle.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(concentrationTransport, 0);
}
}

makeFvSolverOption(concentrationTransport);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::concentrationTransport::setFvDicts()
{
    fvMesh& refMesh = const_cast<fvMesh&>(mesh_);

    // add solver options to relevant dictionaries as required

    dictionary fvSchemes = refMesh.schemes().localSchemeDict();
    dictionary fvSolution = refMesh.solution().localSolutionDict();

    // add time scheme
    word timeType("none");
    if (solverOptions_.found("ddtScheme"))
    {
        timeType = word(solverOptions_.lookup("ddtScheme"));
        string timeSchemeKeyword = "ddt(" + fieldName_ + ")";
        if (!fvSchemes.found(word("ddtSchemes")))
        {
            fvSchemes.add(word("ddtSchemes"), dictionary(), false);
        }
        fvSchemes.subDict("ddtSchemes").add
        (
            word(timeSchemeKeyword),
            ITstream(solverOptions_.lookup("ddtScheme")),
            true
        );
    }

    // add relaxationFactor
    if (!fvSolution.found(word("relaxationFactors")))
    {
        fvSolution.add(word("relaxationFactors"), dictionary(), false);
    }
    else if (!fvSolution.subDict("relaxationFactors").found(word("equations")))
    {
        fvSolution.subDict("relaxationFactors").add
        (
            word("equations"),
            dictionary(),
            false
        );
    }

    if (timeType == "steadyState")
    {

        dictionary& relaxFactor = fvSolution.subDict("relaxationFactors");
        if (relaxFactor.found("equations"))
        {
            relaxFactor.subDict("equations").add
            (
                fieldName_,
                solverOptions_.lookupOrDefault<scalar>("relaxationFactor", 0.7),
                false
            );
        }
        else
        {
            //backward compatibility, just add to base level
            relaxFactor.add
            (
                fieldName_,
                solverOptions_.lookupOrDefault<scalar>("relaxationFactor", 0.7),
                false
            );
        }
    }
    else //transient
    {
        dictionary& relaxFactor = fvSolution.subDict("relaxationFactors");
        if (relaxFactor.found("equations"))
        {
            relaxFactor.subDict("equations").add
            (
                fieldName_,
                solverOptions_.lookupOrDefault<scalar>("relaxationFactor", 1.0),
                false
            );
        }
        else
        {
            //backward compatibility, just add to base level
            relaxFactor.add
            (
                fieldName_,
                solverOptions_.lookupOrDefault<scalar>("relaxationFactor", 1.0),
                false
            );
        }
    }

    // add solver
    if (!fvSolution.found(word("solvers")))
    {
        fvSolution.add(word("solvers"), dictionary(), false);
    }

    if (solverOptions_.found("solver"))
    {
        fvSolution.subDict("solvers").add
        (
            fieldName_,
            solverOptions_.subDict("solver"),
            true
        );
    }
    else
    {
        dictionary solver;
        solver.add("solver", "smoothSolver");
        solver.add("smoother", "GaussSeidel");
        solver.add("tolerance", "1e-6");

        if (timeType == "steadyState")
        {
           solver.add("relTol", "0.1");
        }
        else
        {
           solver.add("relTol", "0");
        }

        fvSolution.subDict("solvers").add
        (
            fieldName_,
            solver,
            false
        );
    }


    //convection scheme
    if (!fvSchemes.found(word("divSchemes")))
    {
        fvSchemes.add(word("divSchemes"), dictionary(), false);
    }
    string convSchemeKeyword = "div(phi,"+ fieldName_ +")";

    if (solverOptions_.found("divScheme"))
    {
        fvSchemes.subDict("divSchemes").add
        (
            word(convSchemeKeyword),
            ITstream(solverOptions_.lookup("divScheme")),
            true
        );
    }
    else
    {
        char pSdiv[] = "Gauss SFCD";

        fvSchemes.subDict("divSchemes").add
        (
            word(convSchemeKeyword),
            pSdiv,
            false
        );
    }

    // Laplacian
    if (!fvSchemes.found(word("laplacianSchemes")))
    {
        fvSchemes.add(word("laplacianSchemes"), dictionary(), false);
    }

    string diffSchemeKeyword
        = "laplacian(Deff,"+ fieldName_ +")";

    if (solverOptions_.found("laplacianScheme"))
    {
        fvSchemes.subDict("laplacianSchemes").add
        (
            word(diffSchemeKeyword),
            ITstream(solverOptions_.lookup("laplacianScheme")),
            true
        );
    }
    else
    {
        char lapScheme[] = "Gauss linear limited 0.333";
        fvSchemes.subDict("laplacianSchemes").add
        (
            word(diffSchemeKeyword),
            lapScheme,
            true
        );
    }

    // gradSchemes
    if (!fvSchemes.found("gradSchemes"))
    {
        fvSchemes.add(word("gradSchemes"), dictionary(), false);
    }
    string gradSchemeKeyword = "grad("+fieldName_+")";
    if (solverOptions_.found("gradScheme"))
    {
        fvSchemes.subDict("gradSchemes").add
        (
            word(gradSchemeKeyword),
            ITstream(solverOptions_.lookup("gradScheme")),
            true
        );
    }
    else
    {
        char Grad[] = "cellLimited Gauss linear 1";
        fvSchemes.subDict("gradSchemes").add
        (
            word(gradSchemeKeyword), Grad, false
        );
    }

    // snGradSchemes
    if (!fvSchemes.found("snGradSchemes"))
    {
        fvSchemes.add(word("snGradSchemes"), dictionary(), false);
    }
    string snGradSchemeKeyword = "snGrad("+fieldName_+")";
    if (solverOptions_.found("snGradScheme"))
    {
        fvSchemes.subDict("snGradSchemes").add
        (
            word(snGradSchemeKeyword),
            ITstream(solverOptions_.lookup("snGradScheme")),
            true
        );
    }
    else
    {
        char snGrad[] = "limited 0.333";
        fvSchemes.subDict("snGradSchemes").add
        (
            word(snGradSchemeKeyword), snGrad, false
        );
    }

    // interpolationSchemes
    if (!fvSchemes.found("interpolationSchemes"))
    {
        fvSchemes.add(word("interpolationSchemes"), dictionary(), false);
    }
    string interpolationSchemeKeyword = "interpolate("+fieldName_+")";
    if (solverOptions_.found("interpolationScheme"))
    {
        fvSchemes.subDict("interpolationSchemes").add
        (
            word(interpolationSchemeKeyword),
            ITstream(solverOptions_.lookup("snGradScheme")),
            true
        );
    }
    else
    {
        char interpolation[] = "linear";
        fvSchemes.subDict("interpolationSchemes").add
        (
            word(interpolationSchemeKeyword), interpolation, false
        );
    }

    refMesh.schemes().setLocalSchemeDict(fvSchemes);
    refMesh.solution().setLocalSolutionDict(fvSolution);
}


tmp<volScalarField> Foam::fv::concentrationTransport::diffusivity()
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    volScalarField& Dt
        = const_cast<volScalarField&>
        (obr_.lookupObject<volScalarField>("Dt"+fieldName_));

    if (obr_.foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const cmpTurbModel& turb =
            obr_.lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

        //Dt.internalField() = (turb.mut()/Sct_())->internalField();
        Dt = turb.mut()/Sct_();

        Dt.correctBoundaryConditions();

        tmp<volScalarField> Deff
        (
            new volScalarField
            (
                IOobject
                (
                    "Deff",
                    mesh_.time().timeName(),
                    obr_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                Dt + turb.rho() * D_()
            )
        );
    // wetAreaRatio corr for humidity
    // do not move as diffusivity is called
    // by fvOptions thermalHumiditySource
    volScalarField& field
        = const_cast<volScalarField&>
        (obr_.lookupObject<volScalarField>(fieldName_));

    forAll(field.boundaryField(), patchI)
    {
        if
        (
            isA< phaseChangeHumidityFvPatchScalarField >
            ( field.boundaryField( )[patchI] )
        )
        {
        Deff.ref().boundaryFieldRef()[patchI] *=
            (
            dynamic_cast<const phaseChangeHumidityFvPatchScalarField& >
            (field.boundaryField()[patchI]).wetAreaRatio()
            );
        }
    }

        return Deff;
    }
    else if (obr_.foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        const incompressible::turbulenceModel& turb =
            obr_.lookupObject<icoTurbModel>(icoTurbModel::propertiesName);

        //Dt.internalField() = (turb.nut()/Sct_())->internalField();
        Dt = turb.nut()/Sct_();

        Dt.correctBoundaryConditions();

        tmp<volScalarField> Deff
        (
            new volScalarField
            (
                IOobject
                (
                    "Deff",
                    mesh_.time().timeName(),
                    obr_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                Dt + D_()
            )
        );

    // wetAreaRatio corr for humidity
    // do not move as diffusivity is called
    // by fvOptions thermalHumiditySource
    volScalarField& field
        = const_cast<volScalarField&>
        (obr_.lookupObject<volScalarField>(fieldName_));

    forAll(field.boundaryField(), patchI)
    {

        if
        (
            isA< phaseChangeHumidityFvPatchScalarField >
            ( field.boundaryField( )[patchI] )
        )
        {
        Deff.ref().boundaryFieldRef()[patchI] *=
            (
            dynamic_cast<const phaseChangeHumidityFvPatchScalarField& >
            (field.boundaryField()[patchI]).wetAreaRatio()
            );
        }
    }

        return Deff;
    }
    else
    {
        //no turb model available, return warning
        WarningInFunction
            << "A valid turbulence model could not be found in the database."
            << nl << "Continue with zero diffusivity" << endl;
        ::abort();
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::concentrationTransport::concentrationTransport
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict),
    fieldName_(dict.lookup("fieldName")),
    phiName_(dict.lookupOrDefault(word("phi"), word("phi"))),
    rhoName_(dict.lookupOrDefault(word("rho"), word("rho"))),
    D_(nullptr),
    Sct_(nullptr),
    solverOptions_(),
    sourceCellSets_(0),
    fixedValueSets_(0),
    fixedCellValues_(0),
    SuSpSets_(0),
    SuSpCoeffs_(0),
    SuSets_(0),
    SuValues_(0),
    gHat_(vector::zero),
    uTerminal_
    (
        dimensionedScalar::lookupOrDefault
        (
            "terminalVelocity",
            dict,
            dimVelocity,
            0.0
         )
    ),
    phaseNamePtr_(),
    solveSpeedup_(0),
    nSubCycles_(0),
    residualPhaseField_("res", dimless, -1),
    fieldDependency_(dict.lookupOrDefault<word>("fieldDependency", "none"))
{
    read(dict);

    Info<< endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::concentrationTransport::~concentrationTransport()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::concentrationTransport::initialise()
{
    // Read molecular diffusivity and turbulent Schmidt number
    const word subDictName(fieldName_ + type());
    if (obr_.foundObject<fluidThermo>(basicThermo::dictName))
    {
        const fluidThermo& thermo =
             obr_.lookupObject<fluidThermo>(basicThermo::dictName);

        if (basicThermo::dictName == basicThermo::matDictName)
        {
            D_.reset
            (
                new dimensionedScalar
                (
                    "D",
                    thermo.optionalSubDict(subDictName).lookup("D")
                )
            );
            Sct_.reset
            (
                new dimensionedScalar
                (
                    "Sct",
                    thermo.optionalSubDict(subDictName).lookup("Sct")
                )
            );
        }
        else
        {
            D_.reset
            (
                new dimensionedScalar
                (
                    "D" + fieldName_,
                    thermo.lookup("D" + fieldName_)
                )
            );
            Sct_.reset
            (
                new dimensionedScalar
                (
                    "Sct" + fieldName_,
                    thermo.lookup("Sct" + fieldName_)
                )
            );
        }
    }
    else if (obr_.foundObject<dictionary>("transportProperties"))
    {
        const dictionary& transport =
             obr_.lookupObject<dictionary>("transportProperties");

        D_.reset
        (
            new dimensionedScalar
            (
                "D" + fieldName_,
                transport.lookup("D" + fieldName_)
            )
        );
        Sct_.reset
        (
            new dimensionedScalar
            (
                "Sct" + fieldName_,
                transport.lookup("Sct" + fieldName_)
            )
        );
    }
    else
    {
        //no turb model available, exit with FatalError
        FatalErrorInFunction
            << "A valid turbulence model could not be found in the database."
            << abort(FatalError);
    }

    return true;
}


void Foam::fv::concentrationTransport::correct
(
    const word& solveName,
    const word& regionName
)
{
    if (debug)
    {
        Info<< "    " << "Solving passive scalar transport variable: "
             << fieldName_ <<  endl;
    }

    fv::options& fvOptions = this->fvOptions();

    volScalarField& field =
        const_cast<volScalarField&>
        (
            obr_.lookupObject<volScalarField>(fieldName_)
        );

    const surfaceScalarField& phi
        = obr_.lookupObject<surfaceScalarField>(phiName_);

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
    if (phaseNamePtr_.valid())
    {
        phaseField
            = obr_.lookupObject<volScalarField>(phaseNamePtr_());
    }

    // modify deltaT
    Time& time = const_cast<Time&>(mesh_.time());
    scalar deltaT = time.deltaTValue();

    if (nSubCycles_ > 1)
    {
        time.setDeltaT(deltaT*solveSpeedup_);

        for
        (
            subCycle<volScalarField> scalarSubCycle(field, nSubCycles_);
            !(++scalarSubCycle).end();
        )
        {
            tmp<volScalarField> D = diffusivity();

            if (phi.dimensions() == dimVelocity * dimArea)
            {
                tmp<fvScalarMatrix> Eq
                (
                    fvm::ddt(phaseField, field)
                  + fvm::div(phi,field,"div(phi," + Foam::word(field.name()) + ")")
                  - fvm::laplacian(phaseField*D,field,"laplacian(Deff," + Foam::word(field.name()) + ")")
                 ==
                    fvOptions(field)
                );

                if (residualPhaseField_.value() > 0)
                {
                    Eq.ref() += fvm::ddt(residualPhaseField_, field)
                              - fvc::ddt(residualPhaseField_, field);
                }

                if (mag(uTerminal_.value()) > SMALL)
                {
                    if (mag(gHat_) > SMALL)
                    {
                        surfaceScalarField phi_mod
                        (
                            "phi_mod",
                            (gHat_ & mesh_.Sf()) * uTerminal_
                        );

                        forAll(phi_mod.boundaryField(), pI)
                        {
                            if (!phi_mod.boundaryField()[pI].coupled())
                            {
                                phi_mod.boundaryFieldRef()[pI] = 0;
                            }
                        }

                        Eq.ref() += fvm::div(phi_mod,field,"div(phi,"+ fieldName_ +")");
                    }
                    else
                    {
                        WarningInFunction
                            << "Finite terminal velocity " << uTerminal_ << " but "
                            << " no gravitational direction specified " << gHat_
                            << endl;
                    }
                }

                Eq->relax();
                fvOptions.constrain(Eq.ref());

                Eq->solve();
                fvOptions.correct(field);
            }
            else //variable density
            {

                const volScalarField& rho
                    = obr_.lookupObject<volScalarField>(rhoName_);

                tmp<fvScalarMatrix> Eq
                (
                    fvm::ddt(rho,field)
                  + fvm::div(phi,field,"div(phi," + Foam::word(field.name()) + ")")
                  - fvm::laplacian(D,field,"laplacian(Deff," + Foam::word(field.name()) + ")")
                 ==
                    fvOptions(rho, field)
                );

                if (mag(uTerminal_.value()) > SMALL)
                {
                    if (mag(gHat_) > SMALL)
                    {
                        surfaceScalarField phi_mod
                        (
                            fvc::interpolate(rho)*((uTerminal_ * gHat_)& mesh_.Sf())
                        );
                        Eq.ref() += fvm::div(phi_mod,field,"div(phi,"+ fieldName_ +")");
                    }
                    else
                    {
                        WarningInFunction
                            << "Finite terminal velocity " << uTerminal_ << " but "
                            << " no gravitational direction specified " << gHat_
                            << endl;
                    }
                }

                Eq->relax();
                fvOptions.constrain(Eq.ref());

                Eq->solve();
                fvOptions.correct(field);
            }
        }
    }
    else
    {
        tmp<volScalarField> D = diffusivity();

        if (phi.dimensions() == dimVelocity * dimArea)
        {
            tmp<fvScalarMatrix> Eq
            (
                fvm::ddt(phaseField, field)
              + fvm::div(phi,field,"div(phi," + Foam::word(field.name()) + ")")
              - fvm::laplacian(phaseField*D,field,"laplacian(Deff," + Foam::word(field.name()) + ")")
             ==
                fvOptions(field)
            );

            if (residualPhaseField_.value() > 0)
            {
                Eq.ref() += fvm::ddt(residualPhaseField_, field)
                          - fvc::ddt(residualPhaseField_, field);
            }

            if (mag(uTerminal_.value()) > SMALL)
            {
                if (mag(gHat_) > SMALL)
                {
                    surfaceScalarField phi_mod
                    (
                        "phi_mod",
                        (gHat_ & mesh_.Sf()) * uTerminal_
                    );

                    forAll(phi_mod.boundaryField(), pI)
                    {
                        if (!phi_mod.boundaryField()[pI].coupled())
                        {
                            phi_mod.boundaryFieldRef()[pI] = 0;
                        }
                    }

                    Eq.ref() += fvm::div(phi_mod,field,"div(phi,"+ fieldName_ +")");
                }
                else
                {
                    WarningInFunction
                        << "Finite terminal velocity " << uTerminal_ << " but "
                        << " no gravitational direction specified " << gHat_
                        << endl;
                }
            }

            Eq->relax();
            fvOptions.constrain(Eq.ref());

            Eq->solve();
            fvOptions.correct(field);
        }
        else //variable density
        {

            const volScalarField& rho
                = obr_.lookupObject<volScalarField>(rhoName_);

            tmp<fvScalarMatrix> Eq
            (
                fvm::ddt(rho,field)
              + fvm::div(phi,field,"div(phi," + Foam::word(field.name()) + ")")
              - fvm::laplacian(D,field,"laplacian(Deff," + Foam::word(field.name()) + ")")
             ==
                fvOptions(rho, field)
            );

            if (mag(uTerminal_.value()) > SMALL)
            {
                if (mag(gHat_) > SMALL)
                {
                    surfaceScalarField phi_mod
                    (
                        fvc::interpolate(rho)*((uTerminal_ * gHat_)& mesh_.Sf())
                    );
                    Eq.ref() += fvm::div(phi_mod,field,"div(phi,"+ fieldName_ +")");
                }
                else
                {
                    WarningInFunction
                        << "Finite terminal velocity " << uTerminal_ << " but "
                        << " no gravitational direction specified " << gHat_
                        << endl;
                }
            }

            Eq->relax();
            fvOptions.constrain(Eq.ref());

            Eq->solve();
            fvOptions.correct(field);
        }
    }
/*
    bound
    (
        field,
        dimensionedScalar(fieldName_, field.dimensions(), 0.0)
    );
*/
}


void Foam::fv::concentrationTransport::write()
{
    //TODO: surely should not check for writeTime here?
    if (mesh_.time().writeTime())
    {
        obr_.lookupObject<volScalarField>(fieldName_).write();
    }

    Info<< endl;
}


void Foam::fv::concentrationTransport::read(const dictionary& dict)
{
    //Read gravity if present
    {
        IOobject fieldHeader
        (
            "g",
            mesh_.time().constant(),
            obr_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );
        if (fieldHeader.typeHeaderOk<uniformDimensionedVectorField>(false))
        {
            uniformDimensionedVectorField g
            (
               fieldHeader
            );
            gHat_ = g.value() / (mag(g.value()) + SMALL);
        }
    }

    //Read terminal velocity if present
    uTerminal_ = dimensionedScalar::lookupOrDefault
    (
        "terminalVelocity",
        dict,
        dimVelocity,
        0.0
    );

    // create field if it does not already exist
    if (!obr_.foundObject<volScalarField>(fieldName_))
    {

        IOobject fieldHeader
        (
            fieldName_,
            mesh_.time().timeName(),
            obr_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        // read if present
        if (fieldHeader.typeHeaderOk<volScalarField>(true))
        {
            autoPtr<volScalarField> field
            (
                new volScalarField(fieldHeader, mesh_)
            );
            field->store(field);
        }
        else //otherwise create from functionObject dict
        {
            fieldHeader.readOpt() = IOobject::NO_READ;

            const dictionary& fieldDict(dict.subDict("fieldDefinition"));

            autoPtr<volScalarField> field
            (
                new volScalarField
                (
                    fieldHeader,
                    mesh_,
                    dimensionSet(fieldDict.lookup("dimensions"))
                )
            );

            field->primitiveFieldRef()
                = scalarField("internalField", fieldDict, mesh_.nCells());

            field->boundaryFieldRef().reset
            (
                volScalarField::Boundary
                (
                    mesh_.boundary(),
                    field(),
                    boundarySetup<scalar>
                    (
                        mesh_,
                        field(),
                        fieldDict
                    )
                )
            );
            field->correctBoundaryConditions();

            field->store(field);
        }
    }


    // Create turbulent diffusivity if it doesn't already exist
    if (!obr_.foundObject<volScalarField>("Dt"+fieldName_))
    {
        autoPtr<volScalarField> field
        (
            new volScalarField
            (
                IOobject
                (
                    "Dt"+fieldName_,
                    mesh_.time().timeName(),
                    obr_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
        field->store(field);
    }


    if (dict.found("sets"))
    {
        PtrList<entry> regions(dict.lookup("sets"));

        sourceCellSets_.setSize(regions.size());

        fixedValueSets_.setSize(regions.size());
        fixedCellValues_.setSize(regions.size());
        SuSets_.setSize(regions.size());
        SuValues_.setSize(regions.size());
        SuSpSets_.setSize(regions.size());
        SuSpCoeffs_.setSize(regions.size());

        label nFixedValue = 0;
        label nFixedSource = 0;
        label nLinearSource = 0;

        forAll(regions, regioni)
        {
            const entry& region = regions[regioni];

            autoPtr<topoSetSource> cellSelector =
                topoSetSource::New(region.keyword(), mesh_, region.dict());

            sourceCellSets_.set
            (
                regioni,
                new cellSet
                (
                    mesh_,
                    "cellSet"+Foam::name(regioni),
                    mesh_.nCells()/10+1  // Reasonable size estimate.
                )
            );

            cellSelector->applyToSet
            (
                topoSetSource::NEW,
                sourceCellSets_[regioni]
            );

            if (region.dict().found("fixedValue"))
            {
                fixedValueSets_[nFixedValue] = regioni;
                fixedCellValues_[nFixedValue++]
                    = readScalar(region.dict().lookup("fixedValue"));
            }
            if (region.dict().found("fixedSource"))
            {
                SuSets_[nFixedSource] = regioni;
                SuValues_[nFixedSource++]
                    = readScalar(region.dict().lookup("fixedSource"));
            }
            if (region.dict().found("linearSource"))
            {
                SuSpSets_[nLinearSource] = regioni;
                SuSpCoeffs_[nLinearSource++]
                    = readScalar(region.dict().lookup("linearSource"));
            }
        }

        fixedValueSets_.setSize(nFixedValue);
        fixedCellValues_.setSize(nFixedValue);
        SuSets_.setSize(nFixedSource);
        SuValues_.setSize(nFixedSource);
        SuSpSets_.setSize(nLinearSource);
        SuSpCoeffs_.setSize(nLinearSource);
    }

    solverOptions_ = dict.subDict("solverOptions");

    // deactivated to read exclusively from
    // system/fvSchemes and system/fvSchemes/fvSolution
    // setFvDicts();

    //check if a phase name has been specified
    if (dict.found("phaseName"))
    {
        phaseNamePtr_.reset(new word(dict.lookup("phaseName")));

        residualPhaseField_ = dimensionedScalar
        (
            "residualPhaseField",
            dimless,
            dict.lookupOrDefault<scalar>("residualPhaseField", 1e-5)
        );
    }

    solveSpeedup_ = dict.lookupOrDefault<label>("solveSpeedup", 1);
    nSubCycles_ = dict.lookupOrDefault<label>("nSubCycles", solveSpeedup_);

    dict.readIfPresent("fieldDependency", fieldDependency_);
}



// ************************************************************************* //
