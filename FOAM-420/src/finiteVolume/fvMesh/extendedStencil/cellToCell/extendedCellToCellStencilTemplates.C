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
    (c) 2013-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fvMesh/extendedStencil/cellToCell/extendedCellToCellStencil.H"
#include "fvMesh/extendedStencil/cellToFace/extendedCellToFaceStencil.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class WeightType>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<WeightType, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
> Foam::extendedCellToCellStencil::weightedSum
(
    const mapDistribute& map,
    const labelListList& stencil,
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    const List<List<WeightType>>& stencilWeights
)
{
    typedef typename outerProduct<WeightType, Type>::type WeightedType;
    typedef GeometricField<WeightedType, fvPatchField, volMesh>
        WeightedFieldType;

    const fvMesh& mesh = fld.mesh();

    // Collect internal and boundary values
    List<List<Type>> stencilFld;
    extendedCellToFaceStencil::collectData(map, stencil, fld, stencilFld);

    tmp<WeightedFieldType> twf
    (
        new WeightedFieldType
        (
            IOobject
            (
                fld.name(),
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensioned<WeightedType>
            (
                fld.name(),
                fld.dimensions(),
                Zero
            )
        )
    );
    WeightedFieldType& wf = twf();

    forAll(wf, celli)
    {
        const List<Type>& stField = stencilFld[celli];
        const List<WeightType>& stWeight = stencilWeights[celli];

        forAll(stField, i)
        {
            wf[celli] += stWeight[i]*stField[i];
        }
    }

    // Boundaries values?

    return twf;
}


// ************************************************************************* //
