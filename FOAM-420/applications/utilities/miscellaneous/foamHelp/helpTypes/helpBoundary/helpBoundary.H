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
    (c) 2012-2016 OpenFOAM Foundation

Class
    Foam::helpBoundary

Description
    This class provides help for boundary conditions (patch fields).  When no
    additional arguments are given, the utility outputs all known boundary
    conditions.

    Usage:
    \plaintable
        \c -browse \<patch type\> | open documentation for patch in browser
        \c -field \<name\>        | boundary conditions for a given field
        \c -constraint            | list constraint boundary conditions
        \c -fixedValue            | list \c fixedValue type boundary conditions
    \endplaintable

Note
    To use the \c -fixedValue option, the \c -field option should also be used

See also
    Foam::helpType
    Foam::fvPatchField

SourceFiles
    helpBoundary.C
    helpBoundaryTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef helpBoundary_H
#define helpBoundary_H

#include "helpTypes/helpType/helpType.H"
#include "db/IOobject/IOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace helpTypes
{

/*---------------------------------------------------------------------------*\
                        Class helpBoundary Declaration
\*---------------------------------------------------------------------------*/

class helpBoundary
:
    public helpType
{

protected:

    // Protected Member Functions

        //- Return/output the available boundary conditions for fields of Type
        template<class Type>
        wordList fieldConditions(const IOobject& io, const bool write) const;

        //- Output the available fixed boundary conditions for fields of Type
        template<class Type>
        void fixedValueFieldConditions(const IOobject& io) const;


public:

    //- Runtime type information
    TypeName("helpBoundary");

    //- Constructor
    helpBoundary();

    //- Destructor
    virtual ~helpBoundary();


    // Member Functions

        //- Initialise - typically setting static variables,
        //  e.g. command line arguments
        virtual void init();

        //- Execute the help
        virtual void execute(const argList& args, const fvMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace helpTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "helpBoundaryTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
