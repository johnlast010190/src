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
    (c) 2022 Esi Ltd.

Class
    Foam::fieldInits::volumeHeight

Description
    volumeHeight initialization method.

SourceFiles
    volumeHeight.C

\*--------------------------------------------------------------------*/

#ifndef volumeHeight_H
#define volumeHeight_H

#include "fieldInitialisation/fieldInit/fieldInit.H"
#include "sets/topoSetSource/topoSetSource.H"
#include "sets/topoSets/cellSet.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fieldInitializations
{

/*--------------------------------------------------------------------*\
                              Class volumeHeight
\*--------------------------------------------------------------------*/

class volumeHeight
:
    public fieldInit
{
    // Private data


    // Private Functions

        //- compute volume given height
        scalar compVolume
        (
            const labelList& cellIds,
            const List<scalar>& heights,
            const scalar maxHeight
        ) const;

        //- set field given height
        void setFieldValue
        (
            const labelList& cellIds,
            const List<scalar>& heights,
            const scalar minHeight = -GREAT,
            const scalar maxHeight = GREAT,
            const scalar setValue = 1.0
        ) const;

        //- set field values for given cellSets
        bool setCellValue
        (
            const PtrList<entry>& setDicts,
            const PtrList<cellSet>& cellSets
        ) const;

public:

    //- Runtime type information
    TypeName("volumeHeight");


    // Constructors

        //- construct from initMethod
        volumeHeight
        (
            const fvMesh& mesh,
            const fvSolutionRegistry& localDb,
            const dictionary& fieldDict,
            const word& fN
        );


    //- Destructor
    virtual ~volumeHeight()
    {}


    // Member Functions

        //- Correct the fieldInit model
        virtual void correct();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fieldInits
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
