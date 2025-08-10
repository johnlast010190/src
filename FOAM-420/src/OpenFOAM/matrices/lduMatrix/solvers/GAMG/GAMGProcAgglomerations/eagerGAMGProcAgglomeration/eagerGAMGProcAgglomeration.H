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
    (c) 2013 OpenFOAM Foundation

Class
    Foam::eagerGAMGProcAgglomeration

Description
    'Eager' processor agglomeration of GAMGAgglomerations: at every
    level agglomerates 'mergeLevels' number of processors onto the
    minimum processor number.

SourceFiles
    eagerGAMGProcAgglomeration.C

\*---------------------------------------------------------------------------*/

#ifndef eagerGAMGProcAgglomeration_H
#define eagerGAMGProcAgglomeration_H

#include "matrices/lduMatrix/solvers/GAMG/GAMGProcAgglomerations/GAMGProcAgglomeration/GAMGProcAgglomeration.H"
#include "containers/Lists/DynamicList/DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class GAMGAgglomeration;

/*---------------------------------------------------------------------------*\
              Class eagerGAMGProcAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class eagerGAMGProcAgglomeration
:
    public GAMGProcAgglomeration
{
    // Private data

        //- Agglpmeration level
        const label mergeLevels_;

        DynamicList<label> comms_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        eagerGAMGProcAgglomeration
        (
            const eagerGAMGProcAgglomeration&
        );

        //- Disallow default bitwise assignment
        void operator=(const eagerGAMGProcAgglomeration&);


public:

    //- Runtime type information
    TypeName("eager");


    // Constructors

        //- Construct given agglomerator and controls
        eagerGAMGProcAgglomeration
        (
            GAMGAgglomeration& agglom,
            const dictionary& controlDict
        );


    //- Destructor
    virtual ~eagerGAMGProcAgglomeration();


    // Member Functions

       //- Modify agglomeration. Return true if modified
        virtual bool agglomerate();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
