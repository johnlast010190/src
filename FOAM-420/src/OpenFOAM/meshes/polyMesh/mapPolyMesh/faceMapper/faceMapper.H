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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::faceMapper

Description
    This object provides mapping and fill-in information for face data
    between the two meshes after the topological change.  It is
    constructed from mapPolyMesh.

SourceFiles
    faceMapper.C

\*---------------------------------------------------------------------------*/

#ifndef faceMapper_H
#define faceMapper_H

#include "meshes/polyMesh/mapPolyMesh/morphFieldMapper.H"
#include "containers/HashTables/HashSet/HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class faceMapper Declaration
\*---------------------------------------------------------------------------*/

class faceMapper
:
    public morphFieldMapper
{
    // Private data

        //- Reference to polyMesh
        const polyMesh& mesh_;

        //- Reference to mapPolyMesh
        const mapPolyMesh& mpm_;

        //- Are there any inserted (unmapped) faces
        bool insertedFaces_;

        //- Is the mapping direct
        bool direct_;


    // Demand-driven private data

        //- Direct addressing (only one for of addressing is used)
        mutable labelList* directAddrPtr_;

        //- Interpolated addressing (only one for of addressing is used)
        mutable labelListList* interpolationAddrPtr_;

        //- Interpolation weights
        mutable scalarListList* weightsPtr_;

        //- Inserted faces
        mutable labelList* insertedFaceLabelsPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        faceMapper(const faceMapper&);

        //- Disallow default bitwise assignment
        void operator=(const faceMapper&);


        //- Calculate addressing for mapping with inserted faces
        void calcAddressing() const;

        //- Clear out local storage
        void clearOut();


public:

    // Constructors

        //- Construct from mapPolyMesh
        faceMapper(const mapPolyMesh& mpm);


    //- Destructor
    virtual ~faceMapper();


    // Member Functions

        //- Return size of field before mapping
        virtual label sizeBeforeMapping() const override;

        //- Return number of internal faces before mapping
        virtual label internalSizeBeforeMapping() const;

        //- Is the mapping direct
        virtual bool direct() const override
        {
            return direct_;
        }

        //- Are there unmapped values? I.e. do all size() elements get value?
        virtual bool hasUnmapped() const override
        {
            return insertedObjects();
        }

        //- Return the direct map addressing
        virtual const labelUList& directAddressing() const override;

        //- Return the interpolated map addressing
        virtual const labelListList& addressing() const override;

        //- Return interpolation weights
        virtual const scalarListList& weights() const override;

        //- Return flux flip map
        virtual const labelHashSet& flipFaceFlux() const;

        //- Return number of old internalFaces
        virtual label nOldInternalFaces() const;

        //- Return old patch starts
        virtual const labelList& oldPatchStarts() const;

        //- Return old patch sizes
        virtual const labelList& oldPatchSizes() const;

        //- Are there any inserted faces
        virtual bool insertedObjects() const override
        {
            return insertedFaces_;
        }

        //- Return list of inserted faces
        virtual const labelList& insertedObjectLabels() const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
