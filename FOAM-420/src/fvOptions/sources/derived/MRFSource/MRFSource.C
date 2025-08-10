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
    (c) 2012-2013 OpenFOAM Foundation
    (c) 2017-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "MRFSource.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrices.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(MRFSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        MRFSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::scalar Foam::fv::MRFSource::getRotationalFlux
(
    const label facei
) const
{
    if (facei >= mesh_.nInternalFaces())
    {
        const label patchi = mesh_.boundaryMesh().whichPatch(facei);
        const label patchFacei = facei - mesh_.boundaryMesh()[patchi].start();
        return getBoundaryRotationalFlux(patchi, patchFacei);
    }

    const face& f = mesh_.faces()[facei];
    if (f.size() > 3 && getFrameSwitches().conservative())
    {
        return getTriangularFlux(f);
    }
    else
    {
        return
            coorFrame().frameVelocity(mesh_.Cf()[facei], false)
          & mesh_.Sf()[facei];
    }
}


Foam::scalar Foam::fv::MRFSource::getBoundaryRotationalFlux
(
    const label patchi,
    const label patchFacei
) const
{
    const label facei = patchFacei + mesh_.boundaryMesh()[patchi].start();
    const face& f = mesh_.faces()[facei];
    if (f.size() > 3 && getFrameSwitches().conservative())
    {
        return getTriangularFlux(f);
    }
    else
    {
        const vector& Cf = mesh_.Cf().boundaryField()[patchi][patchFacei];
        const vector& Sf = mesh_.Sf().boundaryField()[patchi][patchFacei];
        return (coorFrame().frameVelocity(Cf, false)) & Sf;
    }
}


Foam::scalar Foam::fv::MRFSource::getTriangularFlux
(
    const face& f
) const
{
    const pointField& p = mesh_.points();
    scalar totalRotFlux = 0;
    point fCentre = p[f[0]];
    label nPoints = f.size();
    for (label pi = 1; pi < nPoints; pi++)
    {
        fCentre += p[f[pi]];
    }
    for (label pi = 0; pi < nPoints; pi++)
    {
        const point& nextPoint = p[f[(pi + 1) % nPoints]];
        vector c = p[f[pi]] + nextPoint + fCentre;
        c /= 3.;
        vector Sf = 0.5*(nextPoint - p[f[pi]])^(fCentre - p[f[pi]]);
        totalRotFlux += (coorFrame().frameVelocity(c, false) & Sf);
    }
    return totalRotFlux;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::MRFSource::MRFSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    frameSources(name, modelType, dict, obr),
    rhoName_(coeffs_.lookupOrDefault<word>("rhoName", "rho")),
    UName_(coeffs_.lookupOrDefault<word>("UName", "U")),
    isInitialised_(false)
{
    if (!mesh_.conformal())
    {
        FatalErrorInFunction
            << "The " << type() << " fvOption is not compatible with "
            << "non-conformal coupling." << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::MRFSource::initialise()
{
    if (!isInitialised_)
    {
        coorFramePtr_ = coordinateFrame::lookupNew(mesh_, coeffs_);
        coorFrame().registry().attachObject
        (
            this->name(),
            this->frameSourceFaces_.typeName
        );
        volVectorField& U = obr_.lookupObjectRef<volVectorField>(UName_);
        forAll(U.boundaryField(), patchi)
        {
            if (isA<refFvPatch>(U.boundaryField()[patchi]))
            {
                const refFvPatch& frameFvPatch =
                    dynamic_cast<const refFvPatch&>(U.boundaryField()[patchi]);
                frameFvPatch.updateCordinateFrameRegistry();
                U.boundaryFieldRef()[patchi].resetUpdate();
            }
        }
        U.correctBoundaryConditions();
        MRF_ = true;
        fieldNames_.setSize(0);
        applied_.setSize(0);
        isInitialised_ = true;
    }

    return true;
}


void Foam::fv::MRFSource::addSup
(
    fvVectorMatrix& eqn,
    const label fieldI
)
{
    if (eqn.dimensions() == dimForce)
    {
        addAcceleration
        (
            obr_.lookupObject<volScalarField>(rhoName_),
            eqn,
            false
        );
    }
    else
    {
        addAcceleration(eqn, false);
    }
}


void Foam::fv::MRFSource::addSup
(
    const volScalarField& rho,
    fvVectorMatrix& eqn,
    const label fieldI
)
{
    addAcceleration(rho, eqn, false);
}


void Foam::fv::MRFSource::makeRelative(volVectorField& U) const
{
    if (!active_)
    {
        return;
    }

    const volVectorField& C = mesh_.C();

    forAll(cells_, i)
    {
        const label celli = cells_[i];
        U[celli] -= coorFrame().frameVelocity(C[celli], false);
    }

    volVectorField::Boundary& Ubf = U.boundaryFieldRef();
    forAll(frameSourceFaces_.includedFaces(), patchi)
    {
        const labelList& includedFaces =
            frameSourceFaces_.includedFaces()[patchi];
        forAll(includedFaces, i)
        {
            const label patchFacei = includedFaces[i];
            Ubf[patchi][patchFacei] -=
                coorFrame().frameVelocity
                (
                    C.boundaryField()[patchi][patchFacei],
                    false
                );
        }
    }

    forAll(frameSourceFaces_.excludedFaces(), patchi)
    {
        const labelList& excludedFaces =
            frameSourceFaces_.excludedFaces()[patchi];
        forAll(excludedFaces, i)
        {
            const label patchFacei = excludedFaces[i];
            Ubf[patchi][patchFacei] -=
                coorFrame().frameVelocity
                (
                    C.boundaryField()[patchi][patchFacei],
                    false
                );
        }
    }
}


void Foam::fv::MRFSource::makeRelative(surfaceScalarField& phi) const
{
    makeRelativeRhoFlux(geometricOneField(), phi);
}


void Foam::fv::MRFSource::makeRelative
(
    FieldField<fvsPatchField, scalar>& phi
) const
{
    makeRelativeRhoFlux(oneFieldField(), phi);
}


void Foam::fv::MRFSource::makeRelative
(
    Field<scalar>& phi,
    const label patchi
) const
{
    makeRelativeRhoFlux(oneField(), phi, patchi);
}


void Foam::fv::MRFSource::makeRelative
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    makeRelativeRhoFlux(rho, phi);
}


void Foam::fv::MRFSource::makeAbsolute(volVectorField& U) const
{
    if (!active_)
    {
        return;
    }

    const volVectorField& C = mesh_.C();

    forAll(cells_, i)
    {
        const label celli = cells_[i];
        U[celli] += coorFrame().frameVelocity(C[celli], false);
    }

    volVectorField::Boundary& Ubf = U.boundaryFieldRef();

    forAll(frameSourceFaces_.includedFaces(), patchi)
    {
        const labelList& includedFaces =
            frameSourceFaces_.includedFaces()[patchi];
        forAll(includedFaces, i)
        {
            const label patchFacei = includedFaces[i];
            Ubf[patchi][patchFacei] +=
                coorFrame().frameVelocity
                (
                    C.boundaryField()[patchi][patchFacei],
                    false
                );
        }
    }

    forAll(frameSourceFaces_.excludedFaces(), patchi)
    {
        const labelList& excludedFaces =
            frameSourceFaces_.excludedFaces()[patchi];
        forAll(excludedFaces, i)
        {
            const label patchFacei = excludedFaces[i];
            Ubf[patchi][patchFacei] +=
                coorFrame().frameVelocity
                (
                    C.boundaryField()[patchi][patchFacei],
                    false
                );
        }
    }
}


void Foam::fv::MRFSource::makeAbsolute(surfaceScalarField& phi) const
{
    makeAbsoluteRhoFlux(geometricOneField(), phi);
}


void Foam::fv::MRFSource::makeAbsolute
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    makeAbsoluteRhoFlux(rho, phi);
}


void Foam::fv::MRFSource::zero
(
    surfaceScalarField& phi
) const
{
    if (!active_ || getFrameSwitches().ddtPhiCorr())
    {
        return;
    }

    UIndirectList<scalar>(phi, frameSourceFaces_.internalFaces()) = 0.0;

    forAll(frameSourceFaces_.includedFaces(), patchi)
    {
        const labelList& patchFaces =
            frameSourceFaces_.includedFaces()[patchi];
        const scalarField& pphi = phi.boundaryFieldRef()[patchi];
        UIndirectList<scalar>(pphi, patchFaces) = 0.0;
    }
    forAll(frameSourceFaces_.excludedFaces(), patchi)
    {
        const labelList& excludedFaces =
            frameSourceFaces_.excludedFaces()[patchi];
        const scalarField& pphi = phi.boundaryFieldRef()[patchi];
        UIndirectList<scalar>(pphi, excludedFaces) = 0.0;
    }
}


void Foam::fv::MRFSource::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


bool Foam::fv::MRFSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.readIfPresent("rhoName", rhoName_);

        initialise();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
