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
    (c) 2015 Esi Ltd.

\*---------------------------------------------------------------------------*/
#include "dynamicGIBFvMesh/deformingBodyGIBFvMesh/deformingBodyGIBFvMesh.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh
>
void deformingBodyGIBFvMesh::MapGIBField
(
    const fvMesh& mesh,
    const pointField& pF1,
    const pointField& pF2,
    const pointField& pF3,
    const pointField& pFfc,
    const pointField& pFhp,
    const labelList& faceAdd,
    const faceList& fList,
    PrimitivePatchInterpolation<indirectPrimitivePatch>& ppInter
)
{
    HashTable<const GeometricField<Type, PatchField, GeoMesh>*> fields
    (
        mesh.thisDb().objectRegistry::template
            lookupClass<GeometricField<Type, PatchField, GeoMesh>>()
    );

    for
    (
        typename HashTable<const GeometricField<Type, PatchField, GeoMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        GeometricField<Type, PatchField, GeoMesh>& field =
            const_cast<GeometricField<Type, PatchField, GeoMesh>&>
            (*fieldIter());

        typename GeometricField<Type, PatchField, GeoMesh>
        ::GeometricBoundaryField& bfield = field.boundaryField();

        forAll(mesh.boundary(), pI)
        {
            const polyPatch& poly = mesh.boundary()[pI].patch();
            if (isA<indirectPolyPatch>(poly))
            {
                Field<Type>& f = bfield[pI];
                const Field<Type>& fold = bfield[pI].oldTimeField();
                tmp<Field<Type>> pFt
                (
                    new Field<Type>
                    (
                        ppInter.faceToPointInterpolate(fold)
                    )
                );
                Field<Type>& pF = pFt();

                forAll(f, fI)
                {
                    scalar invDis1 = 1/(mag(pF1[fI]-pFhp[fI])+SMALL);
                    scalar invDis2 = 1/(mag(pF2[fI]-pFhp[fI])+SMALL);
                    scalar invDis3 = 1/(mag(pF3[fI]-pFhp[fI])+SMALL);
                    Type invDist =
                        invDis1*pF[fList[fI][0]] +
                        invDis2*pF[fList[fI][1]] +
                        invDis3*pF[fList[fI][2]];
                    invDist /= invDis1 + invDis2 + invDis3;
                    f[fI] = invDist;
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// end namespace FOAM
}
// ************************************************************************* //
