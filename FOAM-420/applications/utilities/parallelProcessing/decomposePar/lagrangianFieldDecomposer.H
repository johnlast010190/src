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

Class
    Foam::lagrangianFieldDecomposer

Description
    Lagrangian field decomposer.

SourceFiles
    lagrangianFieldDecomposer.C
    lagrangianFieldDecomposerDecomposeFields.C

\*---------------------------------------------------------------------------*/

#ifndef lagrangianFieldDecomposer_H
#define lagrangianFieldDecomposer_H

#include "Cloud/Cloud.H"
#include "db/IOobjects/CompactIOField/CompactIOField.H"
#include "indexedParticle/indexedParticle.H"
#include "passiveParticle/passiveParticle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class IOobjectList;

/*---------------------------------------------------------------------------*\
                    Class lagrangianFieldDecomposer Declaration
\*---------------------------------------------------------------------------*/

class lagrangianFieldDecomposer
{
    // Private data

        //- Reference to processor mesh
        const polyMesh& procMesh_;

        //- Lagrangian positions for this processor
        Cloud<passiveParticle> positions_;

        //- The indices of the particles on this processor
        labelList particleIndices_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        lagrangianFieldDecomposer(const lagrangianFieldDecomposer&);

        //- Disallow default bitwise assignment
        void operator=(const lagrangianFieldDecomposer&);


public:

    // Constructors

        //- Construct from components
        lagrangianFieldDecomposer
        (
            const polyMesh& mesh,
            const polyMesh& procMesh,
            const labelList& faceProcAddressing,
            const labelList& cellProcAddressing,
            const word& cloudName,
            const Cloud<indexedParticle>& lagrangianPositions,
            const List<SLList<indexedParticle*>*>& cellParticles
        );


    // Member Functions

        // Read the fields and hold on the pointer list
        template<class Type>
        static void readFields
        (
            const label cloudI,
            const IOobjectList& lagrangianObjects,
            PtrList<PtrList<IOField<Type>>>& lagrangianFields
        );

        template<class Type>
        static void readFieldFields
        (
            const label cloudI,
            const IOobjectList& lagrangianObjects,
            PtrList
            <
                PtrList<CompactIOField<Field<Type>, Type>>
            >& lagrangianFields
        );


        //- Decompose volume field
        template<class Type>
        tmp<IOField<Type>> decomposeField
        (
            const word& cloudName,
            const IOField<Type>& field
        ) const;

        template<class Type>
        tmp<CompactIOField<Field<Type>, Type>> decomposeFieldField
        (
            const word& cloudName,
            const CompactIOField<Field<Type>, Type>& field
        ) const;


        template<class GeoField>
        void decomposeFields
        (
            const word& cloudName,
            const PtrList<GeoField>& fields
        ) const;

        template<class GeoField>
        void decomposeFieldFields
        (
            const word& cloudName,
            const PtrList<GeoField>& fields
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "lagrangianFieldDecomposerDecomposeFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
