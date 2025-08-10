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
    Foam::RegisterSwitch

Description
    Class and registration macros for InfoSwitches and OptimisationSwitches
    to support reading from system/controlDict and dynamic update.

\*---------------------------------------------------------------------------*/

#ifndef registerSwitch_H
#define registerSwitch_H

#include "global/debug/simpleRegIOobject.H"
#include "include/macros.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class RegisterSwitch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class RegisterSwitch
:
    public simpleRegIOobject
{
    Type& optSwitch_;

public:

    RegisterSwitch
    (
        void (*registryFn)(const char* name, simpleRegIOobject*),
        const char* name,
        Type& optSwitch
    )
    :
        simpleRegIOobject(registryFn, name),
        optSwitch_(optSwitch)
    {}

    virtual ~RegisterSwitch()
    {}

    virtual void readData(Foam::Istream& is)
    {
        is >> optSwitch_;
    }

    virtual void writeData(Foam::Ostream& os) const
    {
        os << optSwitch_;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define registerOptSwitch(Name, Type, Switch)                                  \
    static Foam::RegisterSwitch<Type> FILE_UNIQUE(_addToOpt_)                  \
        (Foam::debug::addOptimisationObject, Name, Switch)


#define registerInfoSwitch(Name, Type, Switch)                                 \
    static Foam::RegisterSwitch<Type> FILE_UNIQUE(_addToOpt_)                  \
        (Foam::debug::addInfoObject, Name, Switch)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
