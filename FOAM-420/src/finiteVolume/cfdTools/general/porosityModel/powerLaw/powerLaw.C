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
    (c) 2012-2017 OpenFOAM Foundation
    (c) 2010-2016 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "cfdTools/general/porosityModel/powerLaw/powerLaw.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"
#include "fvMatrices/fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace porosityModels
    {
        defineTypeNameAndDebug(powerLaw, 0);
        addToRunTimeSelectionTable(porosityModel, powerLaw, mesh);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModels::powerLaw::powerLaw
(
    const word& name,
    const word& modelType,
    const objectRegistry& obr,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    porosityModel(name, modelType, obr, mesh, dict, cellZoneName),
    C0_(Function1<scalar>::New("C0", coeffs_)),
    C1_(Function1<scalar>::New("C1", coeffs_)),
    rhoName_(coeffs_.lookupOrDefault<word>("rho", "rho"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porosityModels::powerLaw::~powerLaw()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityModels::powerLaw::calcTransformModelData()
{
    // nothing to be transformed
}


void Foam::porosityModels::powerLaw::calcForce
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& mu,
    vectorField& force
) const
{
    scalarField Udiag(U.size(), 0.0);
    vectorField Usource(U.size(), vector::zero);
    const scalarField& V = mesh_.V();

    apply(Udiag, Usource, V, rho, U);

    force = Udiag*U;
}


void Foam::porosityModels::powerLaw::correct
(
    fvVectorMatrix& UEqn
) const
{
    const volVectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho = obr_.lookupObject<volScalarField>
        (
            IOobject::groupName(rhoName_, U.group())
        );

        apply(Udiag, Usource, V, rho, U);
    }
    else
    {
        apply(Udiag, Usource, V, geometricOneField(), U);
    }
}


void Foam::porosityModels::powerLaw::correct
(
    fvBlockMatrix<vector>& UEqn
) const
{
    const volVectorField& U = UEqn.psi();

    typename CoeffField<vector>::squareTypeField& blockDiag =
        UEqn.diag().asSquare();

    apply(blockDiag, UEqn.source(), geometricOneField(), U);
}


void Foam::porosityModels::powerLaw::correct
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    apply(Udiag, Usource, V, rho, U);
}


void Foam::porosityModels::powerLaw::correct
(
    fvBlockMatrix<vector>& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
) const
{
    const volVectorField& U = UEqn.psi();

    typename CoeffField<vector>::squareTypeField& blockDiag =
        UEqn.diag().asSquare();

    apply(blockDiag, UEqn.source(), rho, U);
}


void Foam::porosityModels::powerLaw::correct
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    const volVectorField& U = UEqn.psi();
    vectorField nullSource(AU.size(), Zero);

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho = obr_.lookupObject<volScalarField>
        (
            IOobject::groupName(rhoName_, U.group())
        );

        apply(AU, nullSource, rho, U);
    }
    else
    {
        apply(AU, nullSource, geometricOneField(), U);
    }
}


void Foam::porosityModels::powerLaw::adjointCorrect
(
    fvVectorMatrix& UaEqn,
    const volVectorField& Uprimal
) const
{
    if (cellZoneIDs_.empty())
    {
        return;
    }

    const scalarField& V = mesh_.V();
    scalarField& Udiag = UaEqn.diag();
    dimensionSet idDim =
        UaEqn.dimensions()*dimTime/dimVolume/UaEqn.psi().dimensions();

    scalar C0 = C0_->value(mesh_.time().timeOutputValue());

    if (idDim == dimDensity)
    {
        const volScalarField& rho =
            obr_.lookupObject<volScalarField>(rhoName_);

        if (C0 > VSMALL)
        {
            adjointApply
            (
                Udiag,
                V,
                rho,
                Uprimal
            );
        }
    }
    else if (idDim == dimless)
    {
        if (C0 > VSMALL)
        {
            adjointApply
            (
                Udiag,
                V,
                geometricOneField(),
                Uprimal
            );
        }
    }
    else
    {
        //unsupported dimensions
        FatalErrorInFunction
            << "Unsupported adjoint matrix dimensions: " << UaEqn.dimensions()
            << exit(FatalError);
    }
}

void Foam::porosityModels::powerLaw::adjointCorrect
(
    fvVectorMatrix& UaEqn,
    const volScalarField& rho,
    const volScalarField& mu,
    const volVectorField& Uprimal
) const
{
    adjointCorrect(UaEqn, Uprimal);
}


void Foam::porosityModels::powerLaw::adjointCorrect
(
    fvBlockMatrix<vector>& UaEqn,
    const volVectorField& U
) const
{
    if (validParentFrame())
    {
        FatalErrorInFunction
           << "Adjoint powerLaw not supported in moving "
           << "reference frame." << exit(FatalError);
    }

    const scalarField& V = mesh_.V();

    const scalar t = db().time().timeOutputValue();
    const scalar C0 = C0_->value(t);
    const scalar C1m1b2 = (C1_->value(t) - 1.0)/2.0;

    typename CoeffField<vector>::squareTypeField& diag =
        UaEqn.diag().asSquare();

    forAll(cellZoneIDs_, zoneI)
    {
        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label cellI = cells[i];
            const tensor dI = V[cellI]*C1m1b2*I*C0*
                pow(magSqr(U[cellI]), C1m1b2-1);
            diag[cellI](0, 0) += dI.xx();
            diag[cellI](0, 1) += dI.yx();
            diag[cellI](0, 2) += dI.zx();
            diag[cellI](1, 0) += dI.xy();
            diag[cellI](1, 1) += dI.yy();
            diag[cellI](1, 2) += dI.zy();
            diag[cellI](2, 0) += dI.xz();
            diag[cellI](2, 1) += dI.yz();
            diag[cellI](2, 2) += dI.zz();
        }
    }
}


bool Foam::porosityModels::powerLaw::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);

    return true;
}


// ************************************************************************* //
