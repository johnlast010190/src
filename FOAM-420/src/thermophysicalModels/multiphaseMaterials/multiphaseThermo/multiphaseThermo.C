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
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "multiphaseThermo.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/unitConversion/unitConversion.H"

//TODO: Most of these includes should be removed in future
#include "mixtures/phaseVolumeFractions/alphasContactAngle/alphasContactAngleFvPatchScalarField.H"
#include "db/Time/Time.H"
#include "finiteVolume/fvc/fvcDiv.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "finiteVolume/fvc/fvcSnGrad.H"
#include "finiteVolume/fvc/fvcFlux.H"
#include "finiteVolume/fvc/fvcMeshPhi.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"

#include "materials/matHeRhoThermo/matHeRhoThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiphaseThermo, 0);
}

makeMatThermo
(
    multiphaseThermo,
    phaseVolumeFractions,
    matHeRhoThermo,
    multiphase
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseThermo::multiphaseThermo
(
    const objectRegistry& obr,
    const word& phaseName
)
:
    rhoThermo(obr, word::null),
    materials_(matLookupOrConstruct(obr, *this)),
    phases_(IOdictionary::lookup<wordList>("phases")),
    fractions_(nullptr),
    mesh_(getMesh(obr)),
    deltaN_
    (
        "deltaN",
        1e-8/pow(average(mesh_.V()), 1.0/3.0)
    )
{
    // We only need sigmas between fluid phases, so allow it to be absent
    if (this->found("sigmas"))
    {
        sigmas_ = sigmaTable(lookup("sigmas"), iNewSigma());
    }

    thermos_.setSize(phases_.size());
    forAll(phases_, phasei)
    {
        thermos_.set(phasei, basicThermo::New(obr, phases_[phasei]).ptr());
        thermos_[phasei].validate(phases_[phasei], "e");
    }

    // Reporting information about loaded models in debug mode
    materials_.reportModelsInfo();
}


Foam::basicThermo& Foam::multiphaseThermo::lookupOrCreate
(
    const objectRegistry& obr,
    const word& phaseName
)
{
    basicThermo& thermo =
        basicThermo::lookupOrCreate(obr, phaseName);

    // If phaseName is specified, try to look up and return the thermo
    // sub-object for the given phase
    if (phaseName != word::null && isA<multiphaseThermo>(thermo))
    {
        multiphaseThermo& mpThermo = refCast<multiphaseThermo>(thermo);

        forAll(mpThermo.thermos(), ti)
        {
            if (mpThermo.phases()[ti] == phaseName)
            {
                return mpThermo.thermos()[ti];
            }
        }
        FatalIOErrorInFunction(mpThermo)
            << "Phase '" << phaseName << "' not found in materials."
            << exit(FatalIOError);
        return mpThermo;
    }
    else
    {
        return thermo;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::multiphaseThermo::init()
{
    // Additional initialisation called from parent's init
    fractions_ =
        &refCast<phaseVolumeFractions>
        (
            this->p_.db().lookupObjectRef<volumeMassFractions>
            (
                "massVolumeFractions"
            )
        );
    multiphaseThermo::correct();
}

const Foam::fvMesh& Foam::multiphaseThermo::getMesh
(
    const objectRegistry& obr
) const
{
    if (isA<fvSolutionRegistry>(obr))
    {
        return dynamic_cast<const fvSolutionRegistry&>(obr).mesh();
    }
    return dynamic_cast<const fvMesh&>(obr);
}


Foam::materialTables& Foam::multiphaseThermo::matLookupOrConstruct
(
    const objectRegistry& obr,
    const dictionary& dict
)
{
    if (!obr.foundObject<objectRegistry>("materialModels"))
    {
        obr.store(new materialTables(obr, dict));
    }
    else if
    (
        !obr.subRegistry("materialModels").foundObject<materialTables>
        (
            "materialTables"
        )
    )
    {
        obr.store(new materialTables(obr, dict));
    }

    return
        obr.subRegistry
        (
            "materialModels"
        ).lookupObjectRef<materialTables>("materialTables");
}


Foam::tmp<Foam::scalarField> Foam::multiphaseThermo::tmpBoundaryField
(
    const scalarField& p,
    const scalarField& T,
    const label patchi,
    const word& modelName
) const
{
    scalarField Told, pOld;
    const bool isTsame = (&T == &this->T_.boundaryField()[patchi]);
    const bool isPsame = (&p == &this->p_.boundaryField()[patchi]);
    if (!isTsame)
    {
        Told = this->T_.boundaryField()[patchi];
        this->T_.boundaryFieldRef()[patchi].forceAssign(T);
    }
    if (!isPsame)
    {
        pOld = this->p_.boundaryField()[patchi];
        this->p_.boundaryFieldRef()[patchi].forceAssign(T);
    }

    tmp<scalarField> tField(materials_(modelName).boundaryField()[patchi]);

    // Reset T/p back to the original
    if (!isTsame)
    {
        this->T_.boundaryFieldRef()[patchi].forceAssign(Told);
    }
    if (!isPsame)
    {
        this->p_.boundaryFieldRef()[patchi].forceAssign(pOld);
    }

    return tField;
}


void Foam::multiphaseThermo::correct()
{
    forAll(thermos_, thermoi)
    {
        thermos_[thermoi].correct();
    }
}


void Foam::multiphaseThermo::correctRho(const volScalarField& dp)
{
    forAll(phases_, phasei)
    {
        if (isA<rhoThermo>(thermos_[phasei]))
        {
            rhoThermo& thermo(refCast<rhoThermo>(thermos_[phasei]));
            thermo.rho() += thermo.psi()*dp;
        }
    }
    rho_ += psi_*dp;
}


bool Foam::multiphaseThermo::incompressible() const
{
    return materials_.incompressible();
}


bool Foam::multiphaseThermo::isochoric() const
{
    return materials_.isochoric();
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseThermo::rCv() const
{
    tmp<volScalarField> trCv
    (
        fractions_->alphas().first()/thermos_.first().Cv()
    );
    for (label phasei = 1; phasei != phases_.size(); ++phasei)
    {
        trCv.ref() += fractions_->alphas()[phasei]/thermos_[phasei].Cv();
    }

    return trCv;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::multiphaseThermo::surfaceTensionForce(const volVectorField& U) const
{
    tmp<surfaceScalarField> tstf
    (
        new surfaceScalarField
        (
            IOobject
            (
                "surfaceTensionForce",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar
            (
                "surfaceTensionForce",
                dimensionSet(1, -2, -2, 0, 0),
                0.0
            )
        )
    );

    surfaceScalarField& stf = tstf.ref();
    stf.setOriented();

    forAll(phases_, phasei)
    {
        const volScalarField& alpha1 = fractions_->alphas()[phasei];

        for (label phase2i = phasei + 1; phase2i != phases_.size(); ++phase2i)
        {
            const volScalarField& alpha2 = fractions_->alphas()[phase2i];

            sigmaTable::const_iterator sigmas =
                sigmas_.find(interfacePair(alpha1, alpha2));

            if (sigmas == sigmas_.end())
            {
                FatalErrorInFunction
                    << "Cannot find interface "
                    << interfacePair(alpha1, alpha2)
                    << " in list of sigma values"
                    << exit(FatalError);
            }

            volScalarField sigma
            (
                "sigma",
                U.db(),
                U.mesh(),
                dimensionSet(1, 0, -2, 0, 0)
            );
            sigma.primitiveFieldRef() = sigmas()->value(this->T());
            sigma.correctBoundaryConditions();

            stf +=
                fvc::interpolate(sigma*K(U, alpha1, alpha2))
               *(
                    fvc::interpolate(alpha2)*fvc::snGrad(alpha1)
                  - fvc::interpolate(alpha1)*fvc::snGrad(alpha2)
                );
        }
    }

    return tstf;
}


Foam::tmp<Foam::surfaceVectorField> Foam::multiphaseThermo::nHatfv
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    surfaceVectorField gradAlphaf
    (
        fvc::interpolate(alpha2)*fvc::interpolate(fvc::grad(alpha1))
      - fvc::interpolate(alpha1)*fvc::interpolate(fvc::grad(alpha2))
    );

    // Face unit interface normal
    return gradAlphaf/(mag(gradAlphaf) + deltaN_);
}


Foam::tmp<Foam::surfaceScalarField> Foam::multiphaseThermo::nHatf
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    return nHatfv(alpha1, alpha2) & mesh_.Sf();
}


void Foam::multiphaseThermo::correctContactAngle
(
    const volVectorField& U,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    surfaceVectorField::Boundary& nHatb
) const
{
    const volScalarField::Boundary& gbf = alpha1.boundaryField();
    const fvBoundaryMesh& boundary = mesh_.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphasContactAngleFvPatchScalarField>(gbf[patchi]))
        {
            const alphasContactAngleFvPatchScalarField& acap =
                refCast<const alphasContactAngleFvPatchScalarField>
                (
                    gbf[patchi]
                );

            vectorField& nHatPatch = nHatb[patchi];

            vectorField AfHatPatch
            (
                mesh_.Sf().boundaryField()[patchi]
               /mesh_.magSf().boundaryField()[patchi]
            );

            alphasContactAngleFvPatchScalarField::thetaPropsTable::
                const_iterator tp =
                acap.thetaProps().find(interfacePair(alpha1, alpha2));

            if (tp == acap.thetaProps().end())
            {
                FatalErrorInFunction
                    << "Cannot find interface " << interfacePair(alpha1, alpha2)
                    << "\n    in table of theta properties for patch "
                    << acap.patch().name()
                    << exit(FatalError);
            }

            bool matched = (tp.key().first() == alpha1.name());

            scalar theta0 = degToRad(tp().theta0(matched));
            scalarField theta(boundary[patchi].size(), theta0);

            scalar uTheta = tp().uTheta();

            // Calculate the dynamic contact angle if required
            if (uTheta > SMALL)
            {
                scalar thetaA = degToRad(tp().thetaA(matched));
                scalar thetaR = degToRad(tp().thetaR(matched));

                // Calculated the component of the velocity parallel to the wall
                vectorField Uwall
                (
                    U.boundaryField()[patchi].patchInternalField()
                  - U.boundaryField()[patchi]
                );
                Uwall -= (AfHatPatch & Uwall)*AfHatPatch;

                // Find the direction of the interface parallel to the wall
                vectorField nWall
                (
                    nHatPatch - (AfHatPatch & nHatPatch)*AfHatPatch
                );

                // Normalise nWall
                nWall /= (mag(nWall) + SMALL);

                // Calculate Uwall resolved normal to the interface parallel to
                // the interface
                scalarField uwall(nWall & Uwall);

                theta += (thetaA - thetaR)*tanh(uwall/uTheta);
            }

            // Reset nHatPatch to correspond to the contact angle
            const scalarField a12(nHatPatch & AfHatPatch);
            const scalarField b1(cos(theta));
            scalarField b2(nHatPatch.size());

            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);
            const scalarField a((b1 - a12*b2)/det);
            const scalarField b((b2 - a12*b1)/det);

            nHatPatch = a*AfHatPatch + b*nHatPatch;

            nHatPatch /= (mag(nHatPatch) + deltaN_.value());
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseThermo::K
(
    const volVectorField& U,
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    tmp<surfaceVectorField> tnHatfv = nHatfv(alpha1, alpha2);

    // Correcting contact angle for boundary alphaContactAngle
    correctContactAngle(U, alpha1, alpha2, tnHatfv.ref().boundaryFieldRef());

    // Simple expression for curvature
    return -fvc::div(tnHatfv & mesh_.Sf());
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseThermo::nearInterface() const
{
    tmp<volScalarField> tnearInt
    (
        new volScalarField
        (
            IOobject
            (
                "nearInterface",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("nearInterface", dimless, 0.0)
        )
    );

    forAll(phases_, phasei)
    {
        tnearInt.ref() =
            max
            (
                tnearInt(),
                pos0(fractions_->alphas()[phasei] - 0.01)
               *pos0(0.99 - fractions_->alphas()[phasei])
            );
    }

    return tnearInt;
}


// ************************************************************************* //
