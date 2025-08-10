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
    (c) 2015-2017 OpenCFD Ltd.

Class
    Foam::functionObjects::readFields

Group
    grpFieldFunctionObjects

Description
    Reads fields from the time directories and adds them to the mesh database
    for further post-processing.

Usage
    Example of function object specification:
    \verbatim
    readFields1
    {
        type        readFields;
        libs        ("libfieldFunctionObjects.so");
        ...
        fields      (U p);
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | type name: readFields   | yes         |
        fields       | list of fields to read  |  no         |
        readOnStart  | flag to start reading on start-up | no  | yes
        log          | Log to standard output  | no          | yes
    \endtable

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    readFields.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_readFields_H
#define functionObjects_readFields_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class readFields Declaration
\*---------------------------------------------------------------------------*/

class readFields
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        //- Fields to load
        wordList fieldSet_;

        //- Flag to read on construction
        bool readOnStart_;


    // Protected Member Functions

        template<class Type>
        bool loadField(const word&);


private:

    // Private member functions

        //- Disallow default bitwise copy construct
        readFields(const readFields&);

        //- Disallow default bitwise assignment
        void operator=(const readFields&);


public:

    //- Runtime type information
    TypeName("readFields");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        readFields
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~readFields();


    // Member Functions

        //- Read the set of fields from dictionary
        virtual bool read(const dictionary&);

        //- Read the fields
        virtual bool execute();

        //- Do nothing
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "readFields/readFieldsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
