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
    (c) 2019 Esi Ltd.
    (c) 2014-2016 OpenFOAM Foundation

Class
    Foam::PairKey

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef PairKey_H
#define PairKey_H

#include "primitives/Pair/Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template <class Type>
class PairKey;

template <class Type>
bool operator==(const PairKey<Type>&, const PairKey<Type>&);

template <class Type>
bool operator!=(const PairKey<Type>&, const PairKey<Type>&);

template <class Type>
Istream& operator>>(Istream&, PairKey<Type>&);

template <class Type>
Ostream& operator<<(Ostream&, const PairKey<Type>&);


/*---------------------------------------------------------------------------*\
                         Class PairKey Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class PairKey
:
    public Pair<Type>
{
public:

        class hash
        :
            public Hash<PairKey<Type>>
        {
        public:

            // Constructors

                // Construct null
                hash();


            // Member operators

                // Generate a hash from a phase pair key
                label operator()(const PairKey<Type>& key) const;
        };


private:

    // Private data

        //- Flag to indicate whether ordering is important
        bool ordered_;


public:

    // Constructors

        //- Construct null
        PairKey();

        //- Construct from names and the ordering flag
        PairKey(const word& name1, const word& name2, const bool ordered);


    // Destructor
    virtual ~PairKey();


    // Friend Operators

        //- Test if keys are equal
        template<class T>
        friend bool operator==(const PairKey<T>& a, const PairKey<T>& b);

        //- Test if keys are unequal
        template<class T>
        friend bool operator!=(const PairKey<T>& a, const PairKey<T>& b);

        //- Read from stdin
        template<class T>
        friend Istream& operator>>(Istream& is, PairKey<T>& key);

        //- Write to stdout
        template<class T>
        friend Ostream& operator<<(Ostream& os, const PairKey<T>& key);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "primitives/Pair/PairKey/PairKeyTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
