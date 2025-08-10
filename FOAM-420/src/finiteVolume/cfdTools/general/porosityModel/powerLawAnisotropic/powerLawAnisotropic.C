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
    (c) 2018 Esi Ltd.
    (c) 2012 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "cfdTools/general/porosityModel/powerLawAnisotropic/powerLawAnisotropic.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"
#include "fvMatrices/fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace porosityModels
    {
        defineTypeNameAndDebug(powerLawAnisotropic, 0);
        addToRunTimeSelectionTable(porosityModel, powerLawAnisotropic, mesh);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModels::powerLawAnisotropic::powerLawAnisotropic
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
    B_(readScalar(coeffs_.lookup("B"))),
    C_(Zero),
    rhoName_(coeffs_.lookupOrDefault<word>("rho", "rho"))
{
    // local-to-global transformation tensor
    const tensor& E = csys().R();

    vector c(coeffs_.lookup("C"));


    C_.xx() = c.x();
    C_.yy() = c.y();
    C_.zz() = c.z();

    C_ = (E & C_ & E.T());

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porosityModels::powerLawAnisotropic::~powerLawAnisotropic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityModels::powerLawAnisotropic::calcTransformModelData()
{
    // nothing to be transformed
}


void Foam::porosityModels::powerLawAnisotropic::calcForce
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

    force = Udiag*U - Usource;
}


void Foam::porosityModels::powerLawAnisotropic::correct
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


void Foam::porosityModels::powerLawAnisotropic::correct
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


void Foam::porosityModels::powerLawAnisotropic::correct
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    const volVectorField& U = UEqn.psi();

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho = obr_.lookupObject<volScalarField>
        (
            IOobject::groupName(rhoName_, U.group())
        );

        apply(AU, rho, U);
    }
    else
    {
        apply(AU, geometricOneField(), U);
    }
}


void Foam::porosityModels::powerLawAnisotropic::adjointCorrect
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
    vectorField& Usource = UaEqn.source();
    const vectorField& U = UaEqn.psi();
    dimensionSet idDim
        = UaEqn.dimensions()*dimTime/dimVolume/UaEqn.psi().dimensions();

    const scalar& B = B_;
    const tensor& C = C_;

    if (magSqr(B) > VSMALL || magSqr(C) > VSMALL)
    {
        if (idDim == dimDensity)
        {
            const volScalarField& rho =
                obr_.lookupObject<volScalarField>(rhoName_);

            adjointApply
            (
                Udiag,
                Usource,
                V,
                rho,
                Uprimal,
                U
            );
        }
        else if (idDim == dimless)
        {
            adjointApply
            (
                Udiag,
                Usource,
                V,
                geometricOneField(),
                Uprimal,
                U
            );
        }
        else
        {
            //unsupported dimensions
            FatalErrorInFunction
                << "Unsupported adjoint matrix dimensions: " << UaEqn.dimensions()
                << exit(FatalError);
        }
    }
}

void Foam::porosityModels::powerLawAnisotropic::adjointCorrect
(
    fvVectorMatrix& UaEqn,
    const volScalarField& rho,
    const volScalarField& mu,
    const volVectorField& Uprimal
) const
{
    adjointCorrect(UaEqn, Uprimal);
}


bool Foam::porosityModels::powerLawAnisotropic::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);

    return true;
}


// ************************************************************************* //
