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
    (c) 2017 OpenCFD Ltd.

Namespace
    Foam::predicates

Description
    Various constant predicate types.

SourceFiles
    predicates.H

\*---------------------------------------------------------------------------*/

#ifndef predicates_H
#define predicates_H

#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace predicates
{

/*---------------------------------------------------------------------------*\
                           Class always Declaration
\*---------------------------------------------------------------------------*/

//- Unary and binary predicates returning true, useful for templating.
class always
{
public:
    typedef always value_type;

    //- Construct null
    inline always()
    {}

    //- Evaluated as a bool - return true
    inline operator bool() const
    {
        return true;
    }

    //- Unary predicate returning true
    template<class T>
    inline bool operator()(const T&) const
    {
        return true;
    }

    //- Binary predicate returning false
    template<class T1, class T2>
    inline bool operator()(const T1&, const T2&) const
    {
        return true;
    }

    //- String match returning true
    inline bool match(const std::string&, bool literal=false) const
    {
        return true;
    }
};


/*---------------------------------------------------------------------------*\
                            Class never Declaration
\*---------------------------------------------------------------------------*/

//- Unary and binary predicates returning false, useful for templating.
class never
{
public:
    typedef never value_type;

    //- Construct null
    inline never()
    {}

    //- Evaluated as a bool - return false
    inline operator bool() const
    {
        return false;
    }

    //- Unary predicate returning false
    template<class T>
    inline bool operator()(const T&) const
    {
        return false;
    }

    //- Binary predicate returning false
    template<class T1, class T2>
    inline bool operator()(const T1&, const T2&) const
    {
        return false;
    }

    //- String match returning false
    inline bool match(const std::string&, bool literal=false) const
    {
        return false;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace predicates

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
