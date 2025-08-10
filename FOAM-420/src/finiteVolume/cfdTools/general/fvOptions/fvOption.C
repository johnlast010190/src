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
    (c) 2010-2023 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/fvOptions/fvOption.H"
#include "fields/volFields/volFields.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<> const char* Foam::NamedEnum< fv::option::execHookType, 6 >::names[] =
    {
        "constrain",
        "correct",
        "operator",
        "outerCorrect",
        "solve",
        "none"
    };

    namespace fv
    {
        defineTypeNameAndDebug(option, 0);
        defineRunTimeSelectionTable(option, dictionary);

        const NamedEnum<fv::option::execHookType, 6>
            fv::option::execHookTypeNames_;
    }
}

const Foam::fvMesh& Foam::fv::option::getMesh(const objectRegistry& obr) const
{
    if (isA<fvSolutionRegistry>(obr))
    {
        return dynamic_cast<const fvSolutionRegistry&>(obr).mesh();
    }
    return dynamic_cast<const fvMesh&>(obr);
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::option::option
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    name_(name),
    modelType_(modelType),
    mesh_(getMesh(obr)),
    obr_(obr),
    dict_(dict),
    coeffs_(dict.optionalSubDict(modelType + "Coeffs")),
    active_(dict_.lookupOrDefault<Switch>("active", true)),
    MRF_(dict_.lookupOrDefault<Switch>("MRF", false)),
    fieldNames_(),
    execHook_(execHookTypeNames_[dict_.lookupOrDefault<word>("hookOp", "none")]),
    applied_(),
    boundarySourcePatchIDs_(),
    boundaryApplied_()
{
    Info<< incrIndent << indent << "Source: " << name_ << endl << decrIndent;
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fv::option> Foam::fv::option::New
(
    const word& name,
    const dictionary& coeffs0,
    const objectRegistry& obr
)
{
    // Backward compatibility
    dictionary coeffs = coeffs0;
    if (word(coeffs.lookup("type")) == "solverOption")
    {
        WarningInFunction
            << "Using deprecated fvOption type 'solverOption'. "
            << "Please use solver type directly as an fvOption." << nl << endl;
        coeffs.merge(coeffs0.subDict("solverOptionCoeffs"));
    }

    word modelType(coeffs.lookup("type"));

    // Backward compatibility
    if (modelType == "radiation")
    {
        WarningInFunction
            << "The 'radiation' fvOption is deprecated. "
            << "Please use 'radiationSolver' instead."
            << endl;
        modelType = "radiationSolver";
    }

    Info<< indent
        << "Selecting finite volume options model type " << modelType << endl;

    const_cast<Time&>(obr.time()).libs().open
    (
        coeffs,
        "libs",
        dictionaryConstructorTable_()
    );

    const auto ctor = ctorTableLookup("Model type", dictionaryConstructorTable_(), modelType);
    return autoPtr<option>(ctor(name, modelType, coeffs, obr));
}


Foam::fv::option::~option()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::option::isActive()
{
    return active_;
}

bool Foam::fv::option::isMRF() const
{
    return MRF_;
}

Foam::label Foam::fv::option::applyToField
(
    const word& fieldName, const word& regionName
) const
{
    forAll(fieldNames_, i)
    {
        if (fieldName == fieldNames_[i])
        {
            if (regionNames_.size())
            {
                if (regionName == regionNames_[i])
                {
                    return i;
                }
            }
            else if
            (
                regionName == (obr_.name() == "" ? mesh_.name() : obr_.name())
            )
            {
                return i;
            }
        }
    }
    return -1;
}

Foam::label Foam::fv::option::applyToBoundaryFieldAndPatch
(
    const word& fieldName,
    const label patchID
) const
{
    if (!boundarySourcePatchIDs_.found(fieldName))
    {
        return -1;
    }
    else
    {
        return findIndex(boundarySourcePatchIDs_[fieldName], patchID);
    }
}


void Foam::fv::option::checkApplied() const
{
    forAll(applied_, i)
    {
        if (!applied_[i])
        {
            WarningInFunction
                << "Source " << name_ << " defined for field "
                << fieldNames_[i] << " but never used" << endl;
        }
    }
}


void Foam::fv::option::checkBoundaryApplied() const
{
    forAllConstIters(boundaryApplied_, iter)
    {
        forAll(iter(), i)
        {
            if (!iter()[i])
            {
                label patchi = boundarySourcePatchIDs_[iter.key()][i];
                WarningInFunction
                    << "Boundary source " << name_ << " defined for field "
                    << iter.key() << ", " << "patch "
                    << mesh_.boundary()[patchi].name()
                    << " but never used" << endl;
            }
        }
    }
}


void Foam::fv::option::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    fvMatrix<sphericalTensor>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    fvMatrix<symmTensor>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    fvMatrix<tensor>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    fvBlockMatrix<vector>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    const volScalarField& rho,
    fvMatrix<sphericalTensor>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    const volScalarField& rho,
    fvMatrix<symmTensor>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    const volScalarField& rho,
    fvMatrix<tensor>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    const volScalarField& rho,
    fvBlockMatrix<vector>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    addSup(alpha*rho, eqn, fieldi);
}


void Foam::fv::option::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    addSup(alpha*rho, eqn, fieldi);
}


void Foam::fv::option::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<sphericalTensor>& eqn,
    const label fieldi
)
{
    addSup(alpha*rho, eqn, fieldi);
}


void Foam::fv::option::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<symmTensor>& eqn,
    const label fieldi
)
{
    addSup(alpha*rho, eqn, fieldi);
}


void Foam::fv::option::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<tensor>& eqn,
    const label fieldi
)
{
    addSup(alpha*rho, eqn, fieldi);
}


void Foam::fv::option::constrain(fvMatrix<scalar>& eqn, const label fieldi)
{}


void Foam::fv::option::constrain(fvMatrix<vector>& eqn, const label fieldi)
{}


void Foam::fv::option::constrain
(
    fvMatrix<sphericalTensor>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::constrain
(
    fvMatrix<symmTensor>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::constrain(fvMatrix<tensor>& eqn, const label fieldi)
{}


void Foam::fv::option::correct(volScalarField& field)
{}


void Foam::fv::option::correct(volVectorField& field)
{}


void Foam::fv::option::correct(volSphericalTensorField& field)
{}


void Foam::fv::option::correct(volSymmTensorField& field)
{}


void Foam::fv::option::correct(volTensorField& field)
{}


void Foam::fv::option::correct()
{}


// ************************************************************************* //
