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
    (c) 2015-2016 OpenFOAM Foundation

Class
    Foam::preservePatchesConstraint

Description
    Constraint to keep owner and neighbour of (cyclic) patch on same
    processor.

SourceFiles
    preservePatchesConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef preservePatchesConstraint_H
#define preservePatchesConstraint_H

#include "decompositionConstraints/decompositionConstraint/decompositionConstraint.H"
#include "primitives/strings/lists/wordReList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace decompositionConstraints
{

/*---------------------------------------------------------------------------*\
                Class preservePatchesConstraint Declaration
\*---------------------------------------------------------------------------*/

class preservePatchesConstraint
:
    public decompositionConstraint
{
    // Private data

        //- List of patches to keep together
        wordReList patches_;


public:

    //- Runtime type information
    TypeName("preservePatches");


    // Constructors

        //- Construct with generic dictionary with optional entry for type
        preservePatchesConstraint
        (
            const dictionary& constraintsDict,
            const word& type
        );

        //- Construct from components
        preservePatchesConstraint(const wordReList& patches);



    //- Destructor
    virtual ~preservePatchesConstraint()
    {}


    // Member Functions

        //- Add my constraints to list of constraints
        virtual void add
        (
            const polyMesh& mesh,
            boolList& blockedFace,
            PtrList<labelList>& specifiedProcessorFaces,
            labelList& specifiedProcessor,
            List<labelPair>& explicitConnections
        ) const;

        //- Apply any additional post-decomposition constraints
        virtual void apply
        (
            const polyMesh& mesh,
            const boolList& blockedFace,
            const PtrList<labelList>& specifiedProcessorFaces,
            const labelList& specifiedProcessor,
            const List<labelPair>& explicitConnections,
            labelList& decomposition
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace decompositionConstraints
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
