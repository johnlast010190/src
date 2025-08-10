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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016 Esi Ltd.

Description

\*---------------------------------------------------------------------------*/

#include "finiteArea/convectionSchemes/gaussFaConvectionScheme/gaussFaConvectionScheme.H"
#include "finiteArea/fac/facEdgeIntegrate.H"
#include "faMatrices/faMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
template<class Type>
const edgeInterpolationScheme<Type>&
gaussConvectionScheme<Type>::interpScheme() const
{
    return tinterpScheme_();
}

template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
gaussConvectionScheme<Type>::flux
(
    const edgeScalarField& faceFlux,
    const GeometricField<Type, faPatchField, areaMesh>& vf
) const
{
    return faceFlux*tinterpScheme_().interpolate(vf);
}


template<class Type>
tmp<faMatrix<Type>>
gaussConvectionScheme<Type>::famDiv
(
    const edgeScalarField& faceFlux,
    const GeometricField<Type, faPatchField, areaMesh>& vf
) const
{
    tmp<edgeScalarField> tweights = tinterpScheme_().weights(vf);
    const edgeScalarField& weights = tweights();

    tmp<faMatrix<Type>> tfam
    (
        new faMatrix<Type>
        (
            vf,
            faceFlux.dimensions()*vf.dimensions()
        )
    );
    faMatrix<Type>& fam = tfam.ref();

    fam.lower() = -weights.internalField()*faceFlux.internalField();
    fam.upper() = fam.lower() + faceFlux.internalField();
    fam.negSumDiag();

    forAll(fam.psi().boundaryField(), patchI)
    {
        const faPatchField<Type>& psf = fam.psi().boundaryField()[patchI];
        const faePatchScalarField& patchFlux = faceFlux.boundaryField()[patchI];
        const faePatchScalarField& pw = weights.boundaryField()[patchI];

        fam.internalCoeffs()[patchI] = patchFlux*psf.valueInternalCoeffs(pw);
        fam.boundaryCoeffs()[patchI] = -patchFlux*psf.valueBoundaryCoeffs(pw);
    }

//     if (tinterpScheme_().corrected())
//     {
//         fam += fac::edgeIntegrate(faceFlux*tinterpScheme_().correction(vf));
//     }

    // Non-euclidian and other corrections
    GeometricField<Type, faePatchField, faEdgeMesh> convFluxCorr
    (
        flux(faceFlux, vf)
      - faceFlux*tinterpScheme_().euclidianInterpolate(vf)
    );

    fam += fac::edgeIntegrate(convFluxCorr);

    return tfam;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
gaussConvectionScheme<Type>::facDiv
(
    const edgeScalarField& faceFlux,
    const GeometricField<Type, faPatchField, areaMesh>& vf
) const
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> tConvection
    (
        fac::edgeIntegrate(flux(faceFlux, vf))
    );

    tConvection.ref().rename
    (
        "convection(" + faceFlux.name() + ',' + vf.name() + ')'
    );

    return tConvection;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
