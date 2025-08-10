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
    (c) 2011 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Description
    Simple central-difference snGrad scheme with implicit bounded
    non-orthogonal correction

\*---------------------------------------------------------------------------*/

#include "finiteVolume/snGradSchemes/SINCSnGrad/SINCSnGrad.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
SINCSnGrad<Type>::~SINCSnGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
SINCSnGrad<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>&
) const
{
    NotImplemented;
}


template<class Type>
tmp<surfaceScalarField> SINCSnGrad<Type>::deltaCoeffs
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Implementation for implicitly limited 'Over-relaxed approach'

    const fvMesh& mesh = this->mesh();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tdeltaCoeffs
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "deltaCoeffs",
                mesh.pointsInstance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false // Do not register
            ),
            mesh.nonOrthDeltaCoeffs()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& deltaCoeffs
        = tdeltaCoeffs.ref();

    const word& gradSchemeName
    (
        defaultGradSchemeName_ == word::null
      ? word("grad(" + vf.name() + ')')
      : defaultGradSchemeName_
    );

    // non-orthogonal correction part
    GeometricField<Type, fvsPatchField, surfaceMesh> nonOrthCorr
    (
        mesh.nonOrthCorrectionVectors()
      & linear<typename outerProduct<vector, Type>::type>(mesh).interpolate
        (
            gradScheme<Type>::New
            (
                mesh,
                vf.db(),
                mesh.schemes().gradScheme(gradSchemeName)
            )().grad(vf, gradSchemeName)
        )
    );

    forAll(owner, facei)
    {
        scalar delta = vf[neighbour[facei]] - vf[owner[facei]];
        if (mag(delta) > SMALL)
        {
            nonOrthCorr[facei] /= delta;
        }
        else
        {
            // If delta is < SMALL nonOrthCorr has to be GREAT.
            // The sign has to be calculated.
            nonOrthCorr[facei] =
                GREAT*(pos0(delta*nonOrthCorr[facei])*2 - 1);
        }
    }

    nonOrthCorr /= dimensionedScalar("p", vf.dimensions(), 1.0);

    // Add orthogonal correction part
    deltaCoeffs += nonOrthCorr;

    // Limiting
    forAll(owner, facei)
    {
        deltaCoeffs[facei] = max(deltaCoeffs[facei], 0);
    }
    return tdeltaCoeffs;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
