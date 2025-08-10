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

#include "finiteArea/laplacianSchemes/gaussFaLaplacianScheme/gaussFaLaplacianScheme.H"
#include "faMesh/faMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeFaLaplacianScheme(gaussLaplacianScheme)

#define declareFamLaplacianScalarGamma(Type)                                 \
                                                                             \
template<>                                                                   \
Foam::tmp<Foam::faMatrix<Foam::Type>>                                       \
Foam::fa::gaussLaplacianScheme<Foam::Type, Foam::scalar>::famLaplacian       \
(                                                                            \
    const GeometricField<scalar, faePatchField, faEdgeMesh>& gamma,         \
    const GeometricField<Type, faPatchField, areaMesh>& vf                    \
)                                                                            \
{                                                                            \
    const faMesh& mesh = this->mesh();                                       \
                                                                             \
    GeometricField<scalar, faePatchField, faEdgeMesh> gammaMagSf            \
    (                                                                        \
        gamma*mesh.magLe()                                                   \
    );                                                                       \
                                                                             \
    tmp<faMatrix<Type>> tfam = famLaplacianUncorrected                      \
    (                                                                        \
        gammaMagSf,                                                          \
        this->tlnGradScheme_().deltaCoeffs(vf),                              \
        vf                                                                   \
    );                                                                       \
    faMatrix<Type>& fam = tfam.ref();                                            \
                                                                             \
    if (this->tlnGradScheme_().corrected())                                  \
    {                                                                        \
        if (mesh.schemes().fluxRequired(vf.name()))                                    \
        {                                                                    \
            fam.faceFluxCorrectionPtr() = new                                \
            GeometricField<Type, faePatchField, faEdgeMesh>                 \
            (                                                                \
                gammaMagSf*this->tlnGradScheme_().correction(vf)             \
            );                                                               \
                                                                             \
            fam.source() -=                                                  \
                mesh.S()*                                                    \
                fac::div                                                     \
                (                                                            \
                    *fam.faceFluxCorrectionPtr()                             \
                )().internalField();                                         \
        }                                                                    \
        else                                                                 \
        {                                                                    \
            fam.source() -=                                                  \
                mesh.S()*                                                    \
                fac::div                                                     \
                (                                                            \
                    gammaMagSf*this->tlnGradScheme_().correction(vf)         \
                )().internalField();                                         \
        }                                                                    \
    }                                                                        \
                                                                             \
    return tfam;                                                             \
}                                                                            \
                                                                             \
                                                                             \
template<>                                                                   \
Foam::tmp<Foam::GeometricField<Foam::Type, Foam::faPatchField, Foam::areaMesh>>\
Foam::fa::gaussLaplacianScheme<Foam::Type, Foam::scalar>::facLaplacian       \
(                                                                            \
    const GeometricField<scalar, faePatchField, faEdgeMesh>& gamma,         \
    const GeometricField<Type, faPatchField, areaMesh>& vf                    \
)                                                                            \
{                                                                            \
    const faMesh& mesh = this->mesh();                                       \
                                                                             \
    tmp<GeometricField<Type, faPatchField, areaMesh>> tLaplacian             \
    (                                                                        \
        fac::div(gamma*this->tlnGradScheme_().lnGrad(vf)*mesh.magLe())       \
    );                                                                       \
                                                                             \
    tLaplacian.ref().rename("laplacian(" + gamma.name() + ',' + vf.name() + ')');\
                                                                             \
    return tLaplacian;                                                       \
}


declareFamLaplacianScalarGamma(scalar);
declareFamLaplacianScalarGamma(vector);
declareFamLaplacianScalarGamma(sphericalTensor);
declareFamLaplacianScalarGamma(symmTensor);
declareFamLaplacianScalarGamma(tensor);


// ************************************************************************* //
