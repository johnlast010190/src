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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::pointFieldDecomposer

Description
    Point field decomposer.

SourceFiles
    pointFieldDecomposer.C
    pointFieldDecomposerDecomposeFields.C

\*---------------------------------------------------------------------------*/

#ifndef pointFieldDecomposer_H
#define pointFieldDecomposer_H

#include "meshes/pointMesh/pointMesh.H"
#include "fields/pointPatchFields/pointPatchField/pointPatchFieldMapperPatchRef.H"
#include "fields/GeometricFields/pointFields/pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class pointFieldDecomposer Declaration
\*---------------------------------------------------------------------------*/

class pointFieldDecomposer
{
public:

    // Public Classes

        //- Point patch field decomposer class
        class patchFieldDecomposer
        :
            public pointPatchFieldMapperPatchRef
        {
            // Private data

                //- Direct map addressing
                labelList directAddressing_;

                //- Does map contain any unmapped values
                bool hasUnmapped_;


        public:

            // Constructors

                //- Construct given addressing
                patchFieldDecomposer
                (
                    const pointPatch& completeMeshPatch,
                    const pointPatch& procMeshPatch,
                    const labelList& directAddr
                );


            // Member functions

                //- Is the mapping direct
                virtual bool direct() const override
                {
                    return true;
                }

                //- No indirect mapping
                virtual bool indirect() const override
                {
                    return false;
                }

                //- Are there unmapped values?
                //  I.e. do all size() elements get value?
                virtual bool hasUnmapped() const override
                {
                    return hasUnmapped_;
                }

                //- Return the direct map addressing
                virtual const labelUList& directAddressing() const override
                {
                    return directAddressing_;
                }
        };


private:

    // Private data

        //- Reference to complete mesh
        const pointMesh& completeMesh_;

        //- Reference to processor mesh
        const pointMesh& procMesh_;

        //- Reference to point addressing
        const labelList& pointAddressing_;

        //- List of patch field decomposers
        PtrList<patchFieldDecomposer> patchFieldDecomposers_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        pointFieldDecomposer(const pointFieldDecomposer&);

        //- Disallow default bitwise assignment
        void operator=(const pointFieldDecomposer&);


public:

    // Constructors

        //- Construct from components
        pointFieldDecomposer
        (
            const pointMesh& completeMesh,
            const pointMesh& procMesh,
            const labelList& pointAddressing
        );


    //- Destructor
    ~pointFieldDecomposer();


    // Member Functions

        //- Decompose point field
        template<class Type>
        tmp<GeometricField<Type, pointPatchField, pointMesh>>
        decomposeField
        (
            const GeometricField<Type, pointPatchField, pointMesh>&
        ) const;

        //- Decompose a list of fields
        template<class GeoField>
        void decomposeFields(const PtrList<GeoField>& fields) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "pointFieldDecomposerDecomposeFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
