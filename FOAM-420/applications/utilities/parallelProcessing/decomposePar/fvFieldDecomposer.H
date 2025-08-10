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
    Foam::fvFieldDecomposer

Description
    Finite Volume volume and surface field decomposer.

SourceFiles
    fvFieldDecomposer.C
    fvFieldDecomposerDecomposeFields.C

\*---------------------------------------------------------------------------*/

#ifndef fvFieldDecomposer_H
#define fvFieldDecomposer_H

#include "fvMesh/fvMesh.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/directFvPatchFieldMapper.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class IOobjectList;

/*---------------------------------------------------------------------------*\
                      Class fvFieldDecomposer Declaration
\*---------------------------------------------------------------------------*/

class fvFieldDecomposer
{
public:

    // Public Classes

        //- Patch field decomposer class
        class patchFieldDecomposer
        :
            public labelList,
            public directFvPatchFieldMapper
        {
        public:

            // Constructors

                //- Construct given addressing
                patchFieldDecomposer(const labelUList& addressing);
        };


        //- Patch field decomposer for indirect patches (GIB)
        class indirectPatchFieldDecomposer
        :
            public generalFvPatchFieldMapper
        {
            // Private Data

                labelList indirectAddressing_;


        public:

            // Constructors

                //- Construct given addressing
                indirectPatchFieldDecomposer
                (
                    const labelUList&,
                    const labelUList&,
                    const labelList&
                );


            // Member Functions

                //- Is the mapping direct
                virtual bool direct() const override
                {
                    return false;
                }

                //- Is the mapping indirect
                virtual bool indirect() const override
                {
                    return true;
                }

                //- Are there unmapped values?
                //  I.e. do all size() elements get value?
                virtual bool hasUnmapped() const override
                {
                    return false;
                }

                //- Return the indirect map addressing
                virtual const labelUList& indirectAddressing() const override
                {
                    return indirectAddressing_;
                }
        };


private:

    // Private Data

        //- Reference to complete mesh
        const fvMesh& completeMesh_;

        //- Reference to processor mesh
        const fvMesh& procMesh_;

        //- Reference to face addressing
        const labelList& faceAddressing_;

        //- Reference to cell addressing
        const labelList& cellAddressing_;

        //- Reference to face addressing boundary field
        const surfaceLabelField::Boundary& faceAddressingBf_;

        //- List of patch field decomposers
        PtrList<patchFieldDecomposer> patchFieldDecomposers_;

        //- List of indirect patch field decomposers
        PtrList<indirectPatchFieldDecomposer> indirectPatchFieldDecomposers_;


    // Private Member Functions

        //- Convert a processor patch to the corresponding complete patch index
        label completePatchID(const label procPatchi) const;

        //- Map cell values to faces
        template<class Type>
        static tmp<Field<Type>> mapCellToFace
        (
            const labelUList& owner,
            const labelUList& neighbour,
            const Field<Type>& field,
            const labelUList& addressing
        );

        //- Map face values to faces
        template<class Type>
        static tmp<Field<Type>> mapFaceToFace
        (
            const Field<Type>& field,
            const labelUList& addressing,
            const bool isFlux
        );

        //- Disallow default bitwise copy construct
        fvFieldDecomposer(const fvFieldDecomposer&);

        //- Disallow default bitwise assignment
        void operator=(const fvFieldDecomposer&);


public:

    // Constructors

        //- Construct from components
        fvFieldDecomposer
        (
            const fvMesh& completeMesh,
            const fvMesh& procMesh,
            const labelList& faceAddressing,
            const labelList& cellAddressing,
            const surfaceLabelField::Boundary& faceAddressingBf
        );


    //- Destructor
    ~fvFieldDecomposer();


    // Member Functions

        //- Decompose volume field
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>>
        decomposeField
        (
            const GeometricField<Type, fvPatchField, volMesh>& field,
            const bool allowUnknownPatchFields = false
        ) const;

        //- Decompose surface field
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        decomposeField
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& field
        ) const;

        //- Decompose a list of fields
        template<class GeoField>
        void decomposeFields(const PtrList<GeoField>& fields) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvFieldDecomposerDecomposeFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
