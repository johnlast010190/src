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

Class
    Foam::bufferedAccumulator

Description

SourceFiles
    bufferedAccumulatorI.H
    bufferedAccumulator.C
    bufferedAccumulatorIO.C

\*---------------------------------------------------------------------------*/

#ifndef bufferedAccumulator_H
#define bufferedAccumulator_H

#include "fields/Fields/Field/Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class bufferedAccumulator;

template<class Type>
Ostream& operator<<
(
    Ostream&,
    const bufferedAccumulator<Type>&
);

/*---------------------------------------------------------------------------*\
                      Class bufferedAccumulator Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class bufferedAccumulator
:
    public List<Field<Type>>
{
    // Private data

        label averagesTaken_;

        List<label> bufferOffsets_;


    // Private Member Functions

        inline Field<Type>& accumulationBuffer();

        inline const Field<Type>& accumulationBuffer() const;

        void accumulateAndResetBuffer(const label b);


public:

    //- Component type
    typedef typename pTraits<Type>::cmptType cmptType;


    // Static data members

        static const char* const typeName;


    // Constructors

        //- Construct null
        bufferedAccumulator();

        //- Construct from components
        bufferedAccumulator
        (
            const label nBuffers,
            const label bufferLength,
            const label bufferingInterval
        );

        //- Construct as copy
        bufferedAccumulator(const bufferedAccumulator<Type>&);


    //- Destructor
    ~bufferedAccumulator();


    // Member Functions

        label addToBuffers(const List<Type>& valuesToAdd);

        Field<Type> averaged() const;

        void resetAveraging();


        // Access

            inline label averagesTaken() const;

            inline label nBuffers() const;

            inline label bufferLength() const;

            inline const List<label>& bufferOffsets() const;


        // Edit

            void setSizes
            (
                const label nBuffers,
                const label bufferLength,
                const label bufferingInterval
            );


    // Member Operators

        void operator=(const bufferedAccumulator<Type>&);


    // IOstream Operators

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const bufferedAccumulator<Type>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bufferedAccumulator/bufferedAccumulatorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "bufferedAccumulator/bufferedAccumulator.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
