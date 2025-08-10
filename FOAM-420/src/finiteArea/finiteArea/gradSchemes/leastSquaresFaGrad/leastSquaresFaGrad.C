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
    (c) 2017 Wikki Ltd.
    (c) 2016 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteArea/gradSchemes/leastSquaresFaGrad/leastSquaresFaGrad.H"
#include "finiteArea/gradSchemes/leastSquaresFaGrad/leastSquaresFaVectors.H"
#include "finiteArea/gradSchemes/gaussFaGrad/gaussFaGrad.H"
#include "faMesh/faMesh.H"
#include "areaMesh/areaMesh.H"
#include "faEdgeMesh/faEdgeMesh.H"
#include "fields/GeometricFields/GeometricField/GeometricField.H"
#include "fields/faPatchFields/basic/zeroGradient/zeroGradientFaPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::faPatchField,
        Foam::areaMesh
    >
>
Foam::fa::leastSquaresFaGrad<Type>::grad
(
    const GeometricField<Type, faPatchField, areaMesh>& vsf
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const faMesh& mesh = vsf.mesh();

    tmp<GeometricField<GradType, faPatchField, areaMesh>> tlsGrad
    (
        new GeometricField<GradType, faPatchField, areaMesh>
        (
            IOobject
            (
                "grad(" + vsf.name() + ')',
                vsf.instance(),
                vsf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<GradType>
            (
                "zero",
                vsf.dimensions()/dimLength,
                pTraits<GradType>::zero
            ),
            zeroGradientFaPatchField<GradType>::typeName
        )
    );
    GeometricField<GradType, faPatchField, areaMesh>& lsGrad = tlsGrad.ref();

    // Get reference to least square vectors
    const leastSquaresFaVectors& lsv = leastSquaresFaVectors::New(mesh);

    const edgeVectorField& ownLs = lsv.pVectors();
    const edgeVectorField& neiLs = lsv.nVectors();

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    forAll(own, edgei)
    {
        label ownEdgeI = own[edgei];
        label neiEdgeI = nei[edgei];

        Type deltaVsf = vsf[neiEdgeI] - vsf[ownEdgeI];

        lsGrad[ownEdgeI] += ownLs[edgei]*deltaVsf;
        lsGrad[neiEdgeI] -= neiLs[edgei]*deltaVsf;
    }

    // Boundary edges
    forAll(vsf.boundaryField(), patchi)
    {
        const faePatchVectorField& patchOwnLs = ownLs.boundaryField()[patchi];

        const labelUList& edgeFaces =
            lsGrad.boundaryField()[patchi].patch().edgeFaces();

        if (vsf.boundaryField()[patchi].coupled())
        {
            Field<Type> neiVsf ( vsf.boundaryField()[patchi].patchNeighbourField() );

            forAll(neiVsf, patchEdgeI)
            {
                lsGrad[edgeFaces[patchEdgeI]] +=
                    patchOwnLs[patchEdgeI]
                   *(neiVsf[patchEdgeI] - vsf[edgeFaces[patchEdgeI]]);
            }
        }
        else
        {
            const faPatchField<Type>& patchVsf = vsf.boundaryField()[patchi];

            forAll(patchVsf, patchEdgeI)
            {
                lsGrad[edgeFaces[patchEdgeI]] +=
                     patchOwnLs[patchEdgeI]
                    *(patchVsf[patchEdgeI] - vsf[edgeFaces[patchEdgeI]]);
            }
        }
    }


    // Remove component of gradient normal to surface (area)
    const areaVectorField& n = vsf.mesh().faceAreaNormals();

    lsGrad -= n*(n & lsGrad);
    lsGrad.correctBoundaryConditions();

    gaussGrad<Type>::correctBoundaryConditions(vsf, lsGrad);

    return tlsGrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
