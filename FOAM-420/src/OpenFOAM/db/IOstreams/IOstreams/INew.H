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
    Foam::INew

Description
    A helper class when constructing from an Istream or dictionary

\*---------------------------------------------------------------------------*/

#ifndef INew_H
#define INew_H

#include "memory/autoPtr/autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;

/*---------------------------------------------------------------------------*\
                           Class INew Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class INew
{

public:

    //- Construct null
    INew()
    {}

    //- Construct from Istream
    autoPtr<T> operator()(Istream& is) const
    {
        return T::New(is);
    }

    //- Construct from word and Istream
    autoPtr<T> operator()(const word&, Istream& is) const
    {
        return T::New(is);
    }

    //- Construct from dictionary
    autoPtr<T> operator()(const dictionary& dict) const
    {
        return T::New(dict);
    }

    //- Construct from word and dictionary
    autoPtr<T> operator()(const word&, const dictionary& dict) const
    {
        return T::New(dict);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
