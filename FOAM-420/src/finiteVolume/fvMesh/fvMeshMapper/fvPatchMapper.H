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
    Foam::fvPatchMapper

Description
    Mapping class for a fvPatchField.

SourceFiles
    fvPatchMapper.C

\*---------------------------------------------------------------------------*/

#ifndef fvPatchMapper_H
#define fvPatchMapper_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "meshes/polyMesh/mapPolyMesh/faceMapper/faceMapper.H"
#include "fvMesh/fvPatches/fvPatch/fvPatch.H"
#include "fields/Fields/primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvPatch;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                        Class fvPatchMapper Declaration
\*---------------------------------------------------------------------------*/

class fvPatchMapper
:
    public generalFvPatchFieldMapper
{
    // Private data

        //- Reference to patch
        const fvPatch& patch_;

        //- Reference to face mapper
        const faceMapper& faceMap_;

        //- Size before mapping
        const label sizeBeforeMapping_;


    // Demand-driven private data

        mutable bool hasUnmapped_;

        //- Direct addressing (only one for of addressing is used)
        mutable labelList* directAddrPtr_;

        //- Interpolated addressing (only one for of addressing is used)
        mutable labelListList* interpolationAddrPtr_;

        //- Interpolation weights
        mutable scalarListList* weightsPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        fvPatchMapper
        (
            const fvPatchMapper&
        );

        //- Disallow default bitwise assignment
        void operator=(const fvPatchMapper&);


        //- Calculate addressing for mapping with inserted cells
        void calcAddressing() const;

        //- Clear out local storage
        void clearOut();


public:

    // Static data members

    // Constructors

        //- Construct from mappers
        fvPatchMapper
        (
            const fvPatch& patch,
            const faceMapper& faceMap
        );


    //- Destructor
    virtual ~fvPatchMapper();


    // Member Functions

        //- Return size of field before mapping
        virtual label sizeBeforeMapping() const
        {
            return sizeBeforeMapping_;
        }

        //- Is the mapping direct
        virtual bool direct() const override
        {
            return faceMap_.direct();
        }

        //- No indirect mapping
        virtual bool indirect() const override
        {
            return false;
        }

        //- Are there unmapped values? I.e. do all size() elements get value?
        virtual bool hasUnmapped() const override
        {
            return hasUnmapped_;
        }

        //- Return the direct map addressing
        virtual const labelUList& directAddressing() const override;

        //- Return the interpolated map addressing
        virtual const labelListList& addressing() const override;

        //- Return interpolation weights
        virtual const scalarListList& weights() const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
