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
    (c) 2010-2019 Esi Ltd.
    (c) 1991-2008 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "passiveScalarTransport.H"
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


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(passiveScalarTransport, 0);
}
}

makeFvSolverOption(passiveScalarTransport);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::passiveScalarTransport::setFvDicts()
{
    fvMesh& refMesh = const_cast<fvMesh&>(mesh_);

    // add solver options to relevant dictionaries as required

    dictionary fvSchemes = refMesh.schemes().localSchemeDict();

    word timeType = "transient";
    string timeSchemeKeyword = "ddt(" + fieldName_ + ")";

    if (solveTime_)
    {
        if (solverOptions_.found("timeScheme"))
        {
            if
            (
                word(solverOptions_.lookup("timeScheme")) == word("steadyState")
            )
            {
                timeType = word("steadyState");
            }
        }
        else
        {
            const dictionary& fvTime = mesh_.schemes().dict().subDict("ddtSchemes");

            if (fvTime.found(timeSchemeKeyword))
            {
                if
                (
                    word(fvTime.lookup(timeSchemeKeyword))
                    == word("steadyState")
                )
                {
                    timeType = word("steadyState");
                }
            }
            else if (fvTime.found("default"))
            {
                if
                (
                    word(fvTime.lookup("default"))
                    == word("steadyState")
                )
                {
                    timeType = word("steadyState");
                }
            }
            else
            {
                FatalErrorInFunction
                    << "Time derivative term activated for "
                     << "solverObject:" << type() << ":" << name_
                     << ", but no valid time scheme is specified"
                     << exit(FatalError);
            }
        }
    }
    else
    {
        timeType = word("steadyState");
    }

    if //replace existing field specific time schemes if present
    (
        solveTime_
    )
    {
        if (!fvSchemes.found(word("ddtSchemes")))
        {
            fvSchemes.add(word("ddtSchemes"), dictionary(), false);
        }

        if (solverOptions_.found("timeScheme"))
        {
            fvSchemes.subDict("ddtSchemes").add
            (
                word(timeSchemeKeyword),
                ITstream(solverOptions_.lookup("timeScheme")),
                true
            );
        }
        else if
        (
            !fvSchemes.subDict("ddtSchemes").found("default")
            || !fvSchemes.subDict("ddtSchemes").found(timeSchemeKeyword)
        )
        {
            fvSchemes.subDict("ddtSchemes").add
            (
                word(timeSchemeKeyword),
                word("Euler"),
                false
            );
        }
    }

    if (timeType == word("steadyState"))
    {
        solveTime_ = false;
    }

    //add relaxation factor and solver
    dictionary fvSolution = refMesh.solution().localSolutionDict();

    if (!solveTime_)
    {
        if (!fvSolution.found(word("relaxationFactors")))
        {
            fvSolution.add(word("relaxationFactors"), dictionary(), false);
            fvSolution.subDict("relaxationFactors").add
            (
                word("equations"),
                dictionary(),
                false
             );
        }
        else if
        (
            !fvSolution.subDict("relaxationFactors").found(word("equations"))
        )
        {
            fvSolution.subDict("relaxationFactors").add
            (
                word("equations"),
                dictionary(),
                false
            );
        }

        if (solverOptions_.found("relaxationFactor"))
        {
            dictionary& relaxFactor = fvSolution.subDict("relaxationFactors");
            relaxFactor.subDict("equations").add
            (
                fieldName_,
                readScalar(solverOptions_.lookup("relaxationFactor")),
                true
            );
        }
        else
        {
            dictionary& relaxFactor = fvSolution.subDict("relaxationFactors");
            relaxFactor.subDict("equations").add
            (
                fieldName_,
                scalar(0.7),
                false
            );
        }
    }
    else
    {
        if (!fvSolution.found(word("relaxationFactors")))
        {
            fvSolution.add(word("relaxationFactors"), dictionary(), false);
            fvSolution.subDict("relaxationFactors").add
            (
                word("equations"),
                dictionary(),
                false
             );
        }
        else if
        (
            !fvSolution.subDict("relaxationFactors").found(word("equations"))
        )
        {
            fvSolution.subDict("relaxationFactors").add
            (
                word("equations"),
                dictionary(),
                false
             );
        }

        dictionary& relaxFactor = fvSolution.subDict("relaxationFactors");
        relaxFactor.subDict("equations").add
        (
            fieldName_,
            scalar(1.0),
            false
        );
    }

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

        if (!solveTime_)
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

    refMesh.solution().setLocalSolutionDict(fvSolution);


    if //replace existing field specific comvection schemes if necessary
    (
        solveConvection_
    )
    {
        if (!fvSchemes.found(word("divSchemes")))
        {
            fvSchemes.add(word("divSchemes"), dictionary(), false);
        }
        string convSchemeKeyword = "div(phi,"+ fieldName_ +")";

        if (solverOptions_.found("convectionScheme"))
        {
            fvSchemes.subDict("divSchemes").add
            (
                word(convSchemeKeyword),
                ITstream(solverOptions_.lookup("convectionScheme")),
                true
            );
        }
        else
        {
            char pSdiv[] = "Gauss limitedLinear 1";

            fvSchemes.subDict("divSchemes").add
            (
                word(convSchemeKeyword),
                pSdiv,
                false
            );

        }
    }

    if //replace existing field specific laplacian schemes if necessary
    (
        solveDiffusion_
    )
    {
        if (!fvSchemes.found(word("laplacianSchemes")))
        {
            fvSchemes.add(word("laplacianSchemes"), dictionary(), false);
        }
        string diffSchemeKeyword1 = "laplacian(nuEff,"+ fieldName_ +")";
        string diffSchemeKeyword2 = "laplacian(muEff,"+ fieldName_ +")";

        if (solverOptions_.found("diffusionScheme"))
        {
            fvSchemes.subDict("laplacianSchemes").add
            (
                word(diffSchemeKeyword1),
                ITstream(solverOptions_.lookup("diffusionScheme")),
                true
            );
            fvSchemes.subDict("laplacianSchemes").add
            (
                word(diffSchemeKeyword2),
                ITstream(solverOptions_.lookup("diffusionScheme")),
                true
            );
        }
    }

    refMesh.schemes().setLocalSchemeDict(fvSchemes);
}

tmp<volScalarField> Foam::fv::passiveScalarTransport::diffusivity()
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;
    tmp<volScalarField> alphaEff;

    if (mesh_.foundObject<cmpTurbModel>(cmpTurbModel::typeName))
    {
        const cmpTurbModel& turb =
            mesh_.lookupObject<cmpTurbModel>(cmpTurbModel::typeName);
        alphaEff = turb.muEff();
    }
    else if (mesh_.foundObject<icoTurbModel>(icoTurbModel::typeName))
    {
        const incompressible::turbulenceModel& turb =
            mesh_.lookupObject<icoTurbModel>(icoTurbModel::typeName);
        alphaEff = turb.nuEff();
    }
    else
    {
        //no turb model available, exit with FatalError
        FatalErrorInFunction
            << "A valid turbulence model could not be found in the database."
            << abort(FatalError);
    }

    return alphaEff;

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::passiveScalarTransport::passiveScalarTransport
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict),
    fieldName_(dict.lookup("fieldName")),
    phiName_(dict.lookupOrDefault(word("phiName"), word("phi"))),
    rhoName_(dict.lookupOrDefault(word("rhoName"), word("rho"))),
    solveTime_
    (
        dict.subDict("solverOptions").lookupOrDefault(word("solveTime"), true)
    ),
    solveConvection_
    (
        dict.subDict("solverOptions").lookupOrDefault
        (
            word("solveConvection"), true
        )
    ),
    solveDiffusion_
    (
        dict.subDict("solverOptions").lookupOrDefault
        (
            word("solveDiffusion"), true
        )
    ),
    bounded01_(false),
    field0_(0),
    solverOptions_(),
    sourceCellSets_(0),
    fixedValueSets_(0),
    fixedCellValues_(0),
    SuSpSets_(0),
    SuSpCoeffs_(0),
    SuSets_(0),
    SuValues_(0),
    fieldDependency_(dict.lookupOrDefault<word>("fieldDependency", "none"))
{
    Info<< "solverObject:" << type() << ":" << name_ << endl;

    read(dict);

    Info<< endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::passiveScalarTransport::~passiveScalarTransport()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::fv::passiveScalarTransport::correct(const word& solveName)
{
    fv::options& fvOptions = this->fvOptions();

    if (debug)
    {
        Info<< "    " << "Solving passive scalar transport variable: "
             << fieldName_ <<  endl;
    }

    volScalarField& field
        = const_cast<volScalarField&>
        (mesh_.lookupObject<volScalarField>(fieldName_));

    const surfaceScalarField& phi
        = mesh_.lookupObject<surfaceScalarField>(phiName_);

    //Info<< field <<endl;

    if (phi.dimensions() == dimVelocity * dimArea)
    {
        tmp<fvScalarMatrix> Eq
        (
            new fvScalarMatrix(field, field.dimensions() * dimVolume / dimTime)
        );

        //f (solveTime_)
        {
            Eq.ref() += fvScalarMatrix(fvm::ddt(field));
        }

        //if (solveConvection_)
        {
            Eq.ref() += fvScalarMatrix(fvm::div(phi,field));
        }

        //if (solveDiffusion_)
        {
            tmp<volScalarField> D = diffusivity();
            Eq.ref() -= fvScalarMatrix(fvm::laplacian(D,field));
        }

        Eq.ref() -= fvOptions(field);

        Eq->relax();

        fvOptions.constrain(Eq.ref());

        Eq->solve();
    }
    else //variable density
    {
        const volScalarField& rho
            = mesh_.lookupObject<volScalarField>(rhoName_);

        tmp<fvScalarMatrix> Eq
        (
            new fvScalarMatrix
            (
                field, rho.dimensions() * field.dimensions()
                * dimVolume / dimTime
            )
        );

        //if (solveTime_)
        {
            Eq.ref() += fvScalarMatrix(fvm::ddt(rho,field));
        }

        //if (solveConvection_)
        {
            Eq.ref() += fvScalarMatrix(fvm::div(phi,field));
        }

        //if (solveDiffusion_)
        {
            tmp<volScalarField> D = diffusivity();
            Eq.ref() -= fvScalarMatrix(fvm::laplacian(D,field));
        }

        Eq.ref() -= fvOptions(rho, field);

        Eq->relax();

        fvOptions.constrain(Eq.ref());

        Eq->solve();
    }

    if (bounded01_)
    {
        bound
        (
            field,
            dimensionedScalar(fieldName_, field.dimensions(), field0_)
        );
    }

}

void Foam::fv::passiveScalarTransport::write()
{
    //TODO: surely this is not needed - it will overwrite the functionObject's writeTime setting?
    if (mesh_.time().outputTime())
    {
        mesh_.lookupObject<volScalarField>(fieldName_).write();
    }

    Info<< endl;
}


void Foam::fv::passiveScalarTransport::read(const dictionary& dict)
{
    //solver options need to be applied every iteration to make sure they
    //are not overwritten by dictionary updates
    solverOptions_ = dict.subDict("solverOptions");

    if (solverOptions_.found("boundValue"))
    {
        bounded01_ = true;
        field0_ = readScalar(solverOptions_.lookup("boundValue"));
    }

    // create field if it does not already exist
    if (!mesh_.foundObject<volScalarField>(fieldName_))
    {

        IOobject fieldHeader
        (
            fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
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

    //check availability of fvSchemes and fvSolution settings
    setFvDicts();

    dict.readIfPresent("fieldDependency", fieldDependency_);
}



// ************************************************************************* //
