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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2010-2012 Esi Ltd.

Class
    printDictionary

Group
    grpReportFunctionObjects

Description
    List of patches and constants to use for lift/drag calculations
    Call printDictionary.write() to calculate&write files.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    printDictionary.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_printDictionary_H
#define functionObjects_printDictionary_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class printDictionary Declaration
\*---------------------------------------------------------------------------*/

class printDictionary
:
    public fvMeshFunctionObject
{
protected:

    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        printDictionary(const printDictionary&);
        void operator=(const printDictionary&);


public:

    //- Runtime type information
    TypeName("printDictionary");


    // Constructors

        //- Construct for given objectRegistry and dictionary
        printDictionary
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    // Destructor

        virtual ~printDictionary();


    // Member Functions

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Read the parameters
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
