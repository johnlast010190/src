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
    Fourth-order snGrad scheme with non-orthogonal correction.

\*---------------------------------------------------------------------------*/

#include "finiteArea/lnGradSchemes/fourthLnGrad/fourthLnGrad.H"
#include "fields/areaFields/areaFields.H"
#include "fields/edgeFields/edgeFields.H"
#include "finiteArea/lnGradSchemes/correctedLnGrad/correctedLnGrad.H"
#include "interpolation/edgeInterpolation/schemes/linear/linearEdgeInterpolation.H"
#include "finiteArea/gradSchemes/gaussFaGrad/gaussFaGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
fourthLnGrad<Type>::~fourthLnGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
fourthLnGrad<Type>::correction
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
) const
{
    const faMesh& mesh = this->mesh();

    tmp<GeometricField<Type, faePatchField, faEdgeMesh>> tcorr
    (
        new GeometricField<Type, faePatchField, faEdgeMesh>
        (
            IOobject
            (
                "lnGradCorr("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            vf.dimensions()*this->mesh().deltaCoeffs().dimensions()
        )
    );
    GeometricField<Type, faePatchField, faEdgeMesh>& corr = tcorr.ref();

    edgeVectorField m ( mesh.Le()/mesh.magLe() );

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        corr.replace
        (
            cmpt,
          - (1.0/15.0)*m
          & linearEdgeInterpolation
            <
                typename
                outerProduct<vector, typename pTraits<Type>::cmptType>::type
            >(mesh).interpolate
            (
                gaussGrad<typename pTraits<Type>::cmptType>(mesh)
               .grad(vf.component(cmpt))
            )
        );
    }

    corr += (1.0/15.0)*correctedLnGrad<Type>(mesh).lnGrad(vf);


//     tmp<GeometricField<Type, faePatchField, faEdgeMesh>> tcorr
//     (
//         (1.0/15.0)
//        *(
//             correctedLnGrad<Type>(mesh).lnGrad(vf)
//           - (
//               linearEdgeInterpolate(gaussGrad<Type>(mesh).grad(vf)) & mesh.Le()
//             )/mesh.magLe()
//         )
//     );

    if (correctedLnGrad<Type>(mesh).corrected())
    {
        tcorr.ref() += correctedLnGrad<Type>(mesh).correction(vf);
    }

    return tcorr;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
