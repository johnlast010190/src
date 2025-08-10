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
    (c) 2011-2017 OpenFOAM Foundation

Class
    Foam::fileFormats::edgeMeshFormatsCore

Description
    A collection of helper functions for reading/writing edge formats.

SourceFiles
    edgeMeshFormatsCore.C

\*---------------------------------------------------------------------------*/

#ifndef edgeMeshFormatsCore_H
#define edgeMeshFormatsCore_H

#include "containers/HashTables/Map/Map.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "primitives/ints/lists/labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class IFstream;
class Time;

namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                       Class edgeMeshFormatsCore Declaration
\*---------------------------------------------------------------------------*/

class edgeMeshFormatsCore
{
protected:

    // Protected Member Functions

        //- Read non-comment line
        static string getLineNoComment(IFstream&);

public:

    // Static Data

        //- The file extension corresponding to 'native' edge format
        //  Normally "eMesh" (edge-mesh)
        static word nativeExt;

    // Static Member Functions

        static bool checkSupport
        (
            const wordHashSet& available,
            const word& ext,
            const bool verbose,
            const word& functionName
        );

       // //- Return the local file name (within time directory)
       // //  NEEDS FIXING
       // static fileName localMeshFileName(const word& edgeName="");

       // //- Find instance with edgeName
       // //  NEEDS FIXING
       // static fileName findMeshInstance
       // (
       //     const Time&,
       //     const word& edgeName=""
       // );

       // //- Find mesh file with edgeName
       // //  NEEDS FIXING
       // static fileName findMeshFile(const Time&, const word& edgeName="");


    // Constructors

        //- Construct null
        edgeMeshFormatsCore();


    //- Destructor
    virtual ~edgeMeshFormatsCore();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
