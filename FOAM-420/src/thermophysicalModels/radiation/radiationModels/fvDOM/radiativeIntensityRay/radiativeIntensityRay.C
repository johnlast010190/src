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
    (c) 2010-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "radiationModels/fvDOM/radiativeIntensityRay/radiativeIntensityRay.H"
#include "finiteVolume/fvm/fvm.H"
#include "radiationModels/fvDOM/fvDOM/fvDOM.H"
#include "global/constants/constants.H"
#include "derivedFvPatchFields/greyDiffusiveRadiation/greyDiffusiveRadiationMixedFvPatchScalarField.H"
#include "containers/Lists/SortableList/SortableList.H"
#include "cfdTools/general/explicitScalarWaveSolve/explicitScalarWaveSolve.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"

using namespace Foam::constant;

const Foam::word
Foam::radiation::radiativeIntensityRay::intensityPrefix("ILambda");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::radiativeIntensityRay::radiativeIntensityRay
(
    const fvDOM& dom,
    const fvMesh& mesh,
    const scalar phi,
    const scalar theta,
    const scalar deltaPhi,
    const scalar deltaTheta,
    const label nLambda,
    const absorptionEmissionModel& absorptionEmission,
    const blackBodyEmission& blackBody,
    const label rayId
)
:
    dom_(dom),
    mesh_(mesh),
    absorptionEmission_(absorptionEmission),
    blackBody_(blackBody),
    d_(Zero),
    dAve_(Zero),
    theta_(theta),
    phi_(phi),
    omega_(0.0),
    nLambda_(nLambda),
    ILambda_(nLambda),
    solutionOrderPtr_(),
    myRayId_(rayId)
{
    scalar sinTheta = Foam::sin(theta);
    scalar cosTheta = Foam::cos(theta);
    scalar sinPhi = Foam::sin(phi);
    scalar cosPhi = Foam::cos(phi);

    omega_ = 2.0*sinTheta*Foam::sin(deltaTheta/2.0)*deltaPhi;
    d_ = vector(sinTheta*sinPhi, sinTheta*cosPhi, cosTheta);
    dAve_ = vector
    (
        sinPhi
       *Foam::sin(0.5*deltaPhi)
       *(deltaTheta - Foam::cos(2.0*theta)
       *Foam::sin(deltaTheta)),
        cosPhi
       *Foam::sin(0.5*deltaPhi)
       *(deltaTheta - Foam::cos(2.0*theta)
       *Foam::sin(deltaTheta)),
        0.5*deltaPhi*Foam::sin(2.0*theta)*Foam::sin(deltaTheta)
    );

    if (mesh_.nSolutionD() == 2)
    {
        vector meshDir(Zero);
        if (dom_.meshOrientation() != vector::zero)
        {
            meshDir = dom_.meshOrientation();
        }
        else
        {
            for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
            {
                if (mesh_.geometricD()[cmpt] == -1)
                {
                    meshDir[cmpt] = 1;
                }
            }
        }
        const vector normal(vector(0, 0, 1));

        const tensor coordRot = rotationTensor(normal, meshDir);

        dAve_ = coordRot & dAve_;
        d_ = coordRot & d_;
    }
    else if (mesh_.nSolutionD() == 1)
    {
        vector meshDir(Zero);
        if (dom_.meshOrientation() != vector::zero)
        {
            meshDir = dom_.meshOrientation();
        }
        else
        {
            for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
            {
                if (mesh_.geometricD()[cmpt] == 1)
                {
                    meshDir[cmpt] = 1;
                }
            }
        }
        const vector normal(vector(1, 0, 0));

        dAve_ = (dAve_ & normal)*meshDir;
        d_ = (d_ & normal)*meshDir;
    }

    autoPtr<volScalarField> IDefaultPtr;

    forAll(ILambda_, lambdaI)
    {
        IOobject IHeader
        (
            intensityPrefix + "_" + name(rayId) + "_" + name(lambdaI),
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        // Check if field exists and can be read
        if (IHeader.typeHeaderOk<volScalarField>(true))
        {
            ILambda_.set
            (
                lambdaI,
                new volScalarField(IHeader, mesh_)
            );
        }
        else
        {
            // Demand driven load the IDefault field
            if (!IDefaultPtr.valid())
            {
                IDefaultPtr.reset
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            "IDefault",
                            mesh_.time().timeName(),
                            mesh_,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh_
                    )
                );
            }

            // Reset the MUST_READ flag
            IOobject noReadHeader(IHeader);
            noReadHeader.readOpt() = IOobject::NO_READ;

            ILambda_.set
            (
                lambdaI,
                new volScalarField(noReadHeader, IDefaultPtr())
            );
        }
    }

    // build solutionOrder
    if (!dom_.participating())
    {
        SortableList<scalar> dist((d_ & (mesh_.C().internalField()))());

        solutionOrderPtr_.reset(new labelList(dist.indices()));
    }


}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::radiativeIntensityRay::~radiativeIntensityRay()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



Foam::scalar Foam::radiation::radiativeIntensityRay::correct()
{
    // Reset boundary heat flux to zero
    qr_.clear();
    qin_.clear();
    qem_.clear();
    qg_.clear();
    qr_.set
    (
        new volScalarField::Boundary
        (
            mesh_.boundary(),
            mesh_.V(),           // Dummy internal field,
            calculatedFvPatchScalarField::typeName
        )
    );
    qin_.set
    (
        new volScalarField::Boundary
        (
            mesh_.boundary(),
            mesh_.V(),           // Dummy internal field,
            calculatedFvPatchScalarField::typeName
        )
    );
    qem_.set
    (
        new volScalarField::Boundary
        (
            mesh_.boundary(),
            mesh_.V(),           // Dummy internal field,
            calculatedFvPatchScalarField::typeName
        )
    );
    qg_.set
    (
        new volScalarField::Boundary
        (
            mesh_.boundary(),
            mesh_.V(),           // Dummy internal field,
            calculatedFvPatchScalarField::typeName
        )
    );
    qr_() = 0.0;
    qin_() = 0.0;
    qem_() = 0.0;
    qg_() = 0.0;

    scalar maxResidual = -GREAT;

    forAll(ILambda_, lambdaI)
    {
        const surfaceScalarField Ji(dAve_ & mesh_.Sf());

        if (dom_.participating())
        {
            const volScalarField& k = dom_.aLambda(lambdaI);

            fvScalarMatrix IiEq
            (
                fvm::div(Ji, ILambda_[lambdaI], "div(Ji,Ii_h)")
              + fvm::Sp(k*omega_, ILambda_[lambdaI])
            ==
                1.0/constant::mathematical::pi*omega_
               *(
                    (k - absorptionEmission_.aDisp(lambdaI))
                   *blackBody_.bLambda(lambdaI)

                  + absorptionEmission_.E(lambdaI)/4

                )
            );

            //look up to see if there is a fvOptions in the registry
            if (mesh_.thisDb().foundObject<fv::options>(fv::options::typeName))
            {
                fv::options& fvOptions = const_cast<fv::options&>
                (
                    mesh_.lookupObject<fv::options>(fv::options::typeName)
                );

                //dummy to pass dimensions check.
                //this is not elegant but I could not find a better way
                dimensionedScalar time_("time_", dimTime, 1.0);
                IiEq -= fvOptions(k*omega_*time_, ILambda_[lambdaI]);
                // IiEq -= fvOptions(ILambda_[lambdaI]);
            }

            IiEq.relax();

            const solverPerformance ILambdaSol = solve
            (
                IiEq,
                mesh_.solution().solver("Ii")
            );

            const scalar initialRes =
                ILambdaSol.initialResidual()*omega_/dom_.omegaMax();

            maxResidual = max(initialRes, maxResidual);
        }
        else //explicit solve for non-participating media
        {
            ILambda_[lambdaI].correctBoundaryConditions();

            surfaceScalarField ILambdaf( fvc::interpolate(ILambda_[lambdaI]) );

            // ordered solution loop
            const labelList& solor(solutionOrderPtr_());

            // call explicit wave solver
            maxResidual = max
            (
                maxResidual,
                Foam::explicitScalarWaveSolve
                (
                    ILambda_[lambdaI], ILambdaf, Ji, solor
                )
            );
        }
    }

    return maxResidual;
}


Foam::tmp<Foam::scalarField>
Foam::radiation::radiativeIntensityRay::emittedRadiantIntensity
(
    label patchI,
    const scalarField& Tp
) const
{
    tmp<scalarField> erfTmp(new scalarField(Tp.size(), 0.0));
    scalarField& erf = erfTmp.ref();

    forAll(ILambda_, lambdaI)
    {
        if
        (
            isA<greyDiffusiveRadiationMixedFvPatchScalarField>
            (ILambda_[lambdaI].boundaryField()[patchI])
        )
        {
            const greyDiffusiveRadiationMixedFvPatchScalarField& gd
                = dynamic_cast<const greyDiffusiveRadiationMixedFvPatchScalarField&>
                (ILambda_[lambdaI].boundaryField()[patchI]);

            erf += gd.emittedRadiantIntensity(Tp);
        }
    }

    return erfTmp;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::radiativeIntensityRay::I() const
{
    tmp<volScalarField> Itmp
    (
        new volScalarField
        (
            IOobject
            (
                rayName(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("I", dimMass/pow3(dimTime), 0.0)
        )
    );

    volScalarField& I = Itmp.ref();

    forAll(ILambda_, lambdaI)
    {
        I += ILambda_[lambdaI];
    }

    return Itmp;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::radiativeIntensityRay::cellQr() const
{
    tmp<volScalarField> Itmp
    (
        new volScalarField
        (
            IOobject
            (
                "cellQr",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("cellQr", dimMass/pow3(dimTime), 0.0)
        )
    );

    volScalarField& I = Itmp.ref();

    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    scalarField cellSurfaceArea(mesh_.nCells(), 0);

    forAll(mesh_.faces(), faceI)
    {
        scalar faceArea = mesh_.magFaceAreas()[faceI];

        if (mesh_.isInternalFace(faceI))
        {
            cellSurfaceArea[owner[faceI]] += faceArea;
            cellSurfaceArea[neighbour[faceI]] += faceArea;
        }
        else
        {
            cellSurfaceArea[owner[faceI]] += faceArea;
        }
    }

    forAll(ILambda_, lambdaI)
    {
        surfaceScalarField ILambdaf( fvc::interpolate(ILambda_[lambdaI]) );

        forAll(mesh_.faces(), faceI)
        {
            if (mesh_.isInternalFace(faceI))
            {
                scalar faceLam = (-dAve() & mesh_.faceAreas()[faceI])
                    * ILambdaf[faceI];
                I[owner[faceI]] += pos0(faceLam) * faceLam
                    /cellSurfaceArea[owner[faceI]];
                I[neighbour[faceI]] += pos0(-faceLam) * -faceLam
                    /cellSurfaceArea[neighbour[faceI]];
            }
            else
            {
                label patchID = mesh_.boundaryMesh( ).whichPatch( faceI ) ;
                const polyPatch& patch = mesh_.boundaryMesh( )[patchID] ;
                scalar faceLam = (-dAve() & mesh_.faceAreas()[faceI]) *
                    ILambda_[lambdaI].boundaryField()
                    [patchID][faceI - patch.start()];

                I[owner[faceI]] += pos0(faceLam) * faceLam
                    /cellSurfaceArea[owner[faceI]];
            }
        }
    }

    return Itmp;
}


// ************************************************************************* //
