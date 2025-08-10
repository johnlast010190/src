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
    (c) 2011 OpenFOAM Foundation

Class
    Foam::patchObserver

Description

    Stores observed pressure fluctuations


SourceFiles
    patchObserver.C

\*---------------------------------------------------------------------------*/

#ifndef patchObserver_H
#define patchObserver_H

#include "acousticAnalogy/observer/observer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class patchObserver Declaration
\*---------------------------------------------------------------------------*/

class patchObserver
:
    public observer
{
    //- Name of patch
    word patchName_;

    //- Name of patch
    label patchID_;

    //- Write on file
    label writeFile_;

    //- Write as field
    label writeAsField_;

    //- Read/create and store pressure field
    void readOrCreatePressure();

    //- Gather to master and scatter to all the patch face centres
    vectorField gatherAndScatterPositions();

public:

    //- Runtime type information
    TypeName("patch");

    // Constructors

        // Construct from components
        patchObserver
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );

    // Destructors
        virtual ~patchObserver(){};

    // Member functions

        //- Return true to write on FO file
        virtual bool write() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
