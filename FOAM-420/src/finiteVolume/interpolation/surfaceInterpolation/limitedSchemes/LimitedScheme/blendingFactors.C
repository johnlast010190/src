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
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace blendingFactors
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline tmp<surfaceScalarField> calculateTheta
(
    const volScalarField& vsf
)
{
    const fvMesh& mesh = vsf.mesh();

    tmp<surfaceScalarField> tTheta
    (
        new surfaceScalarField
        (
            IOobject
            (
                vsf.name() + "Theta",
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimless
        )
    );

    surfaceScalarField& Theta = tTheta.ref();

    const dimensionedScalar deltaN
    (
        "deltaN",
        1e-8/pow(average(mesh.V()), 1.0/3.0)
    );

    // smoothen the field to get better normals
    volScalarField vsfSmooth ( fvc::average(linearInterpolate(vsf)) );
    vsfSmooth = fvc::average(linearInterpolate(vsfSmooth));

    // Calculate interface/cell-face orientation in narrow band around interface
    volVectorField gradcTheta ( fvc::grad(vsfSmooth,"grad(alphaSmooth)") );
    surfaceVectorField gradcThetaf ( linearInterpolate(gradcTheta) );
    surfaceVectorField nHatf ( gradcThetaf/(mag(gradcThetaf) + deltaN) );

    Theta = (nHatf & mesh.Sf())/mesh.magSf();

    return tTheta;
}


// Blending factor for CICSAM:
inline tmp<surfaceScalarField> cosine::operator()
(
    const volScalarField& vsf
) const
{
    const fvMesh& mesh = vsf.mesh();
    tmp<surfaceScalarField> tblendingFactor
    (
        new surfaceScalarField
        (
            IOobject
            (
                vsf.name() + "BlendingFactor",
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimless
        )
    );

    surfaceScalarField& blending = tblendingFactor.ref();

    blending = calculateTheta(vsf);

    scalar kf = 2.0;
    blending = cos(2*(acos(blending)));
    blending = kf * (blending + 1)/2;
    blending = max(min(blending, scalar(1.)), scalar(0.));

    return tblendingFactor;
}


// Blending factor for HRIC:
inline tmp<surfaceScalarField> root2::operator()
(
    const volScalarField& vsf
) const
{
    const fvMesh& mesh = vsf.mesh();
    tmp<surfaceScalarField> tblendingFactor
    (
        new surfaceScalarField
        (
            IOobject
            (
                vsf.name() + "BlendingFactor",
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimless
        )
    );

    surfaceScalarField& blending = tblendingFactor.ref();

    blending = calculateTheta(vsf);

    blending = Foam::pow(mag(blending), 0.5);
    blending = max(min(blending, scalar(1.)), scalar(0.));

    return tblendingFactor;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blendingFactors
} // End namespace Foam

// ************************************************************************* //
