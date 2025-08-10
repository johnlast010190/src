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
    (c) 2022-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "flowSolver.H"
#include "fields/fvPatchFields/constraint/empty/emptyFvPatchField.H"
#include "fields/fvPatchFields/derived/fixedValueZone/fixedValueZoneFvPatchFields.H"
#include "fields/fvsPatchFields/constraint/empty/emptyFvsPatchField.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::fv::flowSolver::correctInactiveGIBZoneFaces
(
    surfaceScalarField& fb
) const
{
    const volScalarField& p = this->p();
    const fvMesh& mesh = this->mesh();

    //- loop in the p field
    //  Check if there is an inactive GIB zone from BC
    //  If yes, grab cells and assign the at their faces fb = 0

    forAll(p.boundaryField(), pI)
    {
        const fvPatchScalarField& pfb = p.boundaryField()[pI];
        if (isA<fixedValueZoneFvPatchField<scalar>>(pfb))
        {
            const indirectPolyPatch& gibPolyPatch =
                refCast<const indirectPolyPatch>(pfb.patch().patch());

            const word czName =
                mesh.cellZones()[gibPolyPatch.zoneId()].name();

            const label& zoneId = mesh.cellZones().findZoneID(czName);
            const labelList& cz = mesh.cellZones()[zoneId];

            const cellList& cells = mesh.cells();
            forAll(cz, czI)
            {
                const label& gcI = cz[czI];
                forAll(cells[gcI], fI)
                {
                    const label gfI = cells[gcI][fI];
                    if (gfI< mesh.nInternalFaces())
                    {
                        fb[gfI] = 0;
                    }
                    else
                    {
                        label patchi = mesh.boundaryMesh().whichPatch(gfI);
                        label lfI = gfI - mesh.boundaryMesh()[patchi].start();
                        if
                        (
                            !isA<emptyFvsPatchField<scalar>>
                            (
                                fb.boundaryField()[patchi]
                            )
                        )
                        {
                            fb.boundaryFieldRef()[patchi][lfI] = 0;
                        }
                    }
                }
            }
            fb.boundaryFieldRef()[pI] = 0;
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::flowSolver::flowSolver
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::flowSolver::initializeGravityHref()
{
    const Time& runTime = mesh().time();
    Info<< "Buoyancy active\n" << endl;
    Info<< "Reading g\n" << endl;
    g_.set
    (
        new uniformDimensionedVectorField
        (
            IOobject
            (
                "g",
                runTime.constant(),
                obr(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    Info<< "Reading hRef\n" << endl;
    hRef_.set
    (
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                obr(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar("hRef", dimLength, 0)
        )
    );
    calculateghFields(true);
}


void Foam::fv::flowSolver::calculateghFields(bool buoyant)
{
    if (!buoyant)
    {
        return;
    }

    dimensionedScalar ghRef
    (
        mag(g_->value()) > SMALL
      ? g_() & (cmptMag(g_->value())/mag(g_->value()))*hRef_()
      : dimensionedScalar("ghRef", g_->dimensions()*dimLength, 0)
    );
    if (!gh_.valid())
    {
        Info<< "Creating field gh\n" << endl;
        gh_.set
        (
            new volScalarField("gh", (g_() & mesh().C()) - ghRef)
        );
        ghf_.set
        (
            new surfaceScalarField("ghf", (g_() & mesh().Cf()) - ghRef)
        );
    }
    else
    {
        gh_() = (g_() & mesh().C()) - ghRef;
        ghf_() = (g_() & mesh().Cf()) - ghRef;
    }

    // Non-orthogonal correction for boundary gradient
    forAll(gh_->boundaryField(), patchi)
    {
        const vectorField pd(mesh().boundary()[patchi].delta());
        const labelUList& fc = mesh().boundary()[patchi].faceCells();
        if
        (
            !isA<emptyFvPatchField<scalar>>(gh_->boundaryField()[patchi])
         && !gh_->boundaryField()[patchi].coupled()
        )
        {
            fvPatchScalarField& pgh = gh_->boundaryFieldRef()[patchi];
            forAll(pgh, bfi)
            {
                pgh[bfi] =
                    (g_->value() & (mesh().C()[fc[bfi]] + pd[bfi]))
                  - ghRef.value();
            }
            ghf_->boundaryFieldRef()[patchi].forceAssign(gh_->boundaryField()[patchi]);
        }
    }
}


Foam::tmp<Foam::surfaceScalarField>
Foam::fv::flowSolver::faceBuoyancyForce(bool includeSnGradP) const
{
    const volScalarField& bRho = buoyantRho();

    // Add and subtract snGrad p so that (p-rho*gh) is lumped into a single
    // snGrad, to minimise clipping by limited scheme
    // Use the same p schemes for all for consistency
    tmp<surfaceScalarField> fbSource =
        - ghf_()
         *fvc::snGrad
          (
              bRho,
              "snGrad("+p().name()+')', "grad("+p().name()+')'
          )
        - fvc::snGrad
          (
              p()-bRho*gh_(),
              "snGrad("+p().name()+')',
              "grad("+p().name()+')'
          );
    // Note: if includeSnGradP==true, this leaves it unbalanced, thereby
    // including it
    if (!includeSnGradP)
    {
        fbSource.ref() += fvc::snGrad(p());
    }
    correctInactiveGIBZoneFaces(fbSource.ref());

    return fbSource;
}


void Foam::fv::flowSolver::setOrComputeRhof()
{
    if (!thermo().isochoric() || !rhof_.valid() || !isStatic())
    {
        rhof_.clear();
        rhof_.reset
        (
            new surfaceScalarField
            (
                "rhof", rhoInterpolation()
            )
        );
    }
}


const Foam::tmp<Foam::surfaceScalarField>
Foam::fv::flowSolver::rhoInterpolation(const volScalarField& rho) const
{
    if (solnControl().transonic())
    {
        return
            fvc::interpolate
            (
                rho,
                phiv(),
                "interpolate("+phiv().name()+","+rho.name()+")"
            );
    }
    else
    {
        return fvc::interpolate(rho);
    }
}


// ************************************************************************* //
