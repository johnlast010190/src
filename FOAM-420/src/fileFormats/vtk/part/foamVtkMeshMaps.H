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
    (c) 2017 OpenCFD Ltd.

Class
    Foam::foamVtkMeshMaps

Description
    Bookkeeping for mesh subsetting and/or polyhedral cell decomposition.
    Although the main use case is for handling vtk meshes, it is not specific
    to VTK alone.

    The cellMap is a local-to-global lookup for normal and decomposed cells.
    The pointMap is an optional local-to-global lookup for point ids.
    The additional ids is typically used to store the cell-centre labels
    for additional points of decomposed cells

SourceFiles
    foamVtkMeshMapsI.H

\*---------------------------------------------------------------------------*/

#ifndef foamVtkMeshMaps_H
#define foamVtkMeshMaps_H

#include "containers/Lists/DynamicList/DynamicList.H"
#include "primitives/ints/lists/labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class foamVtkMeshMaps Declaration
\*---------------------------------------------------------------------------*/

class foamVtkMeshMaps
{
    // Private Data

        //- Original cell ids for all cells (regular and decomposed)
        DynamicList<label> cellMap_;

        //- Point labels for subsetted meshes
        DynamicList<label> pointMap_;

        //- Any additional (user) labels.
        //  Eg, cell-centre labels for additional points of decomposed cells
        DynamicList<label> additionalIds_;


public:

    // Constructors

        //- Construct null
        inline explicit foamVtkMeshMaps(const label size = 0);


    //- Destructor
    ~foamVtkMeshMaps() = default;


    // Member Functions

    // Access

        //- Original cell ids for all cells (regular and decomposed).
        //  A regular mesh comprising only primitive cell types, this will just
        //  be an identity list. However, for subsetted meshes and decomposed
        //  cells this becomes a useful means of mapping from the original mesh.
        inline const labelList& cellMap() const;

        //- Point labels for subsetted meshes
        inline const labelList& pointMap() const;

        //- Any additional (user) labels.
        //  Eg, cell-centre labels for additional points of decomposed cells
        inline const labelList& additionalIds() const;


    // Edit

        //- Clear
        inline void clear();

        //- Renumber cell ids (cellMap and additionalIds) to account for
        //  subset meshes
        void renumberCells(const labelUList& mapping);

        //- Renumber point ids (pointMap) to account for subset meshes
        void renumberPoints(const labelUList& mapping);


        //- Original cell ids for all cells (regular and decomposed).
        //  For a regular mesh comprising only primitive cell types, this
        //  will simply be an identity list. However, for subsetted meshes
        //  and decomposed cells this becomes a useful means of mapping from
        //  the original mesh.
        inline DynamicList<label>& cellMap();

        //- Point labels for subsetted meshes
        inline DynamicList<label>& pointMap();

        //- Any additional (user) labels.
        //  Eg, cell-centre labels for additional points of decomposed cells
        inline DynamicList<label>& additionalIds();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "vtk/part/foamVtkMeshMapsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
