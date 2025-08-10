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
    (c) 2010-2012, 2020 Esi Ltd

Class
    Foam::functionObjects::fieldWrite

Group
    grpFieldFunctionObjects

Description
    List of patches and constants to use for lift/drag calculations
    Call fieldWrite.write() to calculate&write files.

SourceFiles
    fieldWrite.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fieldWrite_H
#define functionObjects_fieldWrite_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class fieldWrite Declaration
\*---------------------------------------------------------------------------*/

class fieldWrite
:
    public fvMeshFunctionObject
{
    // Protected data


public:

        enum spanControls
        {
            scStartAndEnd,
            scStartAndNwrites,
            scStart,
            scEnd,
            scNone
        };


private:

            static const NamedEnum<spanControls, 5> spanControlNames_;

        // Read from dictonary

            //- input
            spanControls spanControl_;
            label maxWrites_;

            scalar startWrite_;
            scalar endWrite_;

            wordList fieldNames_;

        //- Result values
            label nWrites_;
            bool outputTime_;
            scalar currentValue_;
            label outputIndex_;

    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        fieldWrite(const fieldWrite&);
        void operator=(const fieldWrite&);


public:

    //- Runtime type information
    TypeName("fieldWrite");


    // Constructors

        //- Construct for given objectRegistry and dictionary
        fieldWrite
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    // Destructor

        virtual ~fieldWrite();


    // Member Functions

        //- Update solver
        void calculate();

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
