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
    (c) 2017 OpenFOAM Foundation
    (c) 2020-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "radiationSolver.H"
#include "fluidThermo/fluidThermo.H"
#include "fvMatrices/fvMatrices.H"
#include "radiationModels/fvDOM/fvDOM/fvDOM.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "derivedFvPatchFields/greyDiffusiveRadiation/greyDiffusiveRadiationMixedFvPatchScalarField.H"
#include "radiationModels/noRadiation/noRadiation.H"
#include "incompressible/transportModel/transportModel.H"


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(radiationSolver, 0);
}
}

makeFvSolverOption(radiationSolver);


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void Foam::fv::radiationSolver::updateQrConst()
{
    if (!qrConst_.valid())
    {
        qrConst_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "qrConst",
                    mesh_.time().timeName(),
                    this->db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_.V(), // Dummy internal field
                volScalarField::Boundary
                (
                    mesh_.boundary(),
                    mesh_.V(),           // Dummy internal field,
                    calculatedFvPatchScalarField::typeName
                )
            )
        );
        qrConst_().boundaryFieldRef().forceAssign(scalar(0));
    }

    radiation::radiationModel& radiation(*radiationPtr_);
    const Foam::radiation::fvDOM& dom =
        refCast<Foam::radiation::fvDOM>(radiation);

    forAll(dom.IRayLambda(0,0).boundaryField(), patchi)
    {
        const Foam::radiation::fvDOM& dom =
            refCast<Foam::radiation::fvDOM>(radiation);
        if
        (
            isA
            <Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField>
            (
                dom.IRayLambda(0,0).boundaryField()[patchi]
            )
        )
        {
            const scalarField& pT0 = TPtr_->boundaryField()[patchi];
            scalarField eri0
            (
                dom.emittedRadiantIntensity(patchi, pT0)
            );
            qrConst_().boundaryFieldRef()[patchi].forceAssign(
                dom.qr().boundaryField()[patchi] + eri0
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::radiationSolver::radiationSolver
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict),
    radiationPtr_(nullptr),
    thermoPtr_(nullptr),
    transportPtr_(nullptr),
    TPtr_(nullptr),
    rhoCpRef_
    (
        "rhoCpRef",
        dimDensity*dimEnergy/dimMass/dimTemperature,
        1.0
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::radiationSolver::initialise()
{
    // Support legacy use of transportProperties
    transportPtr_ =
        obr_.lookupObjectPtr<transportModel>("transportProperties");
    if (transportPtr_)
    {
        TPtr_ = obr_.lookupObjectPtr<volScalarField>(addPhaseName("T"));
        if (!TPtr_)
        {
            WarningInFunction
                << "Could not look up T field " << addPhaseName("T")
                << " - disabling radiation"
                << nl << endl;
            return false;
        }
    }
    else
    {
        thermoPtr_ = &basicThermo::lookupOrCreate(obr_);
        TPtr_ = &thermoPtr_->T();
    }
    radiationPtr_ =
        &radiation::radiationModel::lookupOrCreate(*TPtr_);

    if
    (
        transportPtr_
     && radiationPtr_->radiation()
     && radiationPtr_->participating()
    )
    {
        dimensionedScalar rhoRef
        (
            "rho",
            dimDensity,
            refCast<const IOdictionary>(*transportPtr_)
        );
        dimensionedScalar CpRef
        (
            "Cp",
            dimSpecificHeatCapacity,
            refCast<const IOdictionary>(*transportPtr_)
        );
        rhoCpRef_ = rhoRef*CpRef;
    }

    // Opt out if no radiation model present, but radiation model created
    // above is kept in object-registry as it may be needed for possible
    // inter-region queries
    return !isA<radiation::noRadiation>(*radiationPtr_);
}


void Foam::fv::radiationSolver::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    solveNames.append("radiation");
    correctorMembers.insert(solverObject::outerCorrectorName, solveNames);
}


void Foam::fv::radiationSolver::correct
(
    const word& solveName,
    const word& regionName
)
{
    // Returns true if recalculation was done this iteration
    radiation::radiationModel& radiation(*radiationPtr_);
    if (radiation.correct())
    {
        // Record the constant part of the radiation for use in boundary source
        if (isA<Foam::radiation::fvDOM>(radiation))
        {
            updateQrConst();
        }
    }
}


void Foam::fv::radiationSolver::getSourceGraph
(
    wordList& fields,
    HashTable<wordList>& sourceDependencies
)
{
    if (thermoPtr_ && thermoPtr_->solvesTFromhe())
    {
        fields = {thermoPtr_->he().name()};
    }
    else
    {
        fields = {TPtr_->name()};
    }
    sourceDependencies.insert(fields[0], {"radiation"});
}


void Foam::fv::radiationSolver::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    eqn += radiationPtr_->ST(rhoCpRef_, const_cast<volScalarField&>(eqn.psi()));
}


void Foam::fv::radiationSolver::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (thermoPtr_->solvesTFromhe())
    {
        // Equation is for he
        eqn += radiationPtr_->Sh(*thermoPtr_, eqn.psi());
    }
    else
    {
        // Equation is for T
        eqn +=
            radiationPtr_->ST
            (
                dimensionedScalar("one", dimless, 1),
                const_cast<volScalarField&>(eqn.psi())
            );
    }
}


void Foam::fv::radiationSolver::getBoundarySourceGraph
(
    HashTable<labelList>& fieldPatchIDs,
    HashTable<wordList>& boundarySourceDependencies
)
{
    DynamicList<label> patchIDs;
    forAll(mesh_.boundary(), i)
    {
        if (isA<wallFvPatch>(mesh_.boundary()[i]))
        {
            patchIDs.append(i);
        }
    }
    if (thermoPtr_)
    {
        fieldPatchIDs.insert(thermoPtr_->T().name(), patchIDs);
        boundarySourceDependencies.insert(thermoPtr_->T().name(), {"radiation"});
    }
    else
    {
        fieldPatchIDs.insert(TPtr_->name(), patchIDs);
        boundarySourceDependencies.insert(TPtr_->name(), {"radiation"});
    }
}


void Foam::fv::radiationSolver::addBoundarySource
(
    const word& fieldName,
    const label patchID,
    const scalarField& pf,
    scalarField& f,
    scalarField& df
)
{
    // We can provide the derivative for DOM radiation, otherwise just
    // treat it like a constant source
    radiation::radiationModel& radiation(*radiationPtr_);
    if (radiation.radiation())
    {
        if (isA<Foam::radiation::fvDOM>(radiation))
        {
            const Foam::radiation::fvDOM& dom =
                refCast<Foam::radiation::fvDOM>(radiation);
            if
            (
                isA
                <Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField>
                (
                    dom.IRayLambda(0,0).boundaryField()[patchID]
                )
            )
            {
                if (!qrConst_.valid())
                {
                    // In the legacy case, may not have been corrected yet
                    updateQrConst();
                }
                scalarField eri( dom.emittedRadiantIntensity(patchID, pf) );
                f += qrConst_().boundaryField()[patchID] - eri;
                df -= 4*eri/pf;
            }
        }
        else
        {
            if (mesh_.foundObject<volScalarField>("Qr"))
            {
                volScalarField Qr = mesh_.lookupObject<volScalarField>("Qr");
                f += Qr.boundaryField()[patchID];
            }
        }
    }
}


// ************************************************************************* //
