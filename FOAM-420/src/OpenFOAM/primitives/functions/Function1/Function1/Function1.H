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
    (c) 2018-2023 OpenCFD Ltd.
    (c) 2019-2022 Esi Ltd.

Class
    Foam::Function1

Description
    Top level data entry class for use in dictionaries. Provides a mechanism
    to specify a variable as a certain type, e.g. constant or table, and
    provide functions to return the (interpolated) value, and integral between
    limits.

SourceFiles
    Function1.C
    Function1New.C

\*---------------------------------------------------------------------------*/

#ifndef Function1_H
#define Function1_H

#include "db/dictionary/dictionary.H"
#include "fields/Fields/Field/Field.H"
#include "memory/refPtr/refPtr.H"
#include "containers/HashTables/HashPtrTable/HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class Time;

// Forward declaration of friend functions and operators
template<class Type> class Function1;
template<class Type> Ostream& operator<<(Ostream&, const Function1<Type>&);

/*---------------------------------------------------------------------------*\
                         Class Function1 Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Function1
:
    public tmp<Function1<Type>>::refCount
{
    // Private Member Functions

        //- Selector, with alternative entry, fallback redirection, etc
        static autoPtr<Function1<Type>> New
        (
            const word& entryName,   // Entry name for function
            const entry* eptr,       // Eg, dict.findEntry(entryName)
            const dictionary& dict,
            const word& redirectType,
            const bool mandatory
        );

        //- Disallow default bitwise assignment
        void operator=(const Function1<Type>&);


protected:

    // Protected data

        //- Name of entry
        const word name_;


public:

    //- Runtime type information
    TypeName("Function1")

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        Function1,
        dictionary,
        (
            const word& entryName,
            const dictionary& dict
        ),
        (entryName, dict)
    );


    // Constructor

        //- Construct from entry name
        Function1(const word& entryName);

        //- Copy constructor
        Function1(const Function1<Type>& de);

        //- Construct and return a clone
        virtual tmp<Function1<Type>> clone() const
        {
            return tmp<Function1<Type>>(new Function1<Type>(*this));
        }


    // Selectors

        //- Selector
        static autoPtr<Function1<Type>> New
        (
            const word& entryName,
            const dictionary& dict
        );

        //- An optional selector, without fallback redirection
        static autoPtr<Function1<Type>> NewIfPresent
        (
            const word& entryName,
            const dictionary& dict
        );

    // Caching Selectors - accept wildcards in dictionary

        //- Selector with external storage of Function1.
        //- This also allows wildcard matches in a dictionary
        //static autoPtr<Function1<Type>> New
        static refPtr<Function1<Type>> New
        (
            HashPtrTable<Function1<Type>>& cache,
            const word& entryName,
            const dictionary& dict,
            enum keyType::option matchOpt = keyType::LITERAL,
            const bool mandatory = true
        );


    //- Destructor
    virtual ~Function1();


    // Member Functions

        // Access

            //- Return the name of the entry
            const word& name() const;


        // Manipulation

            //- Convert time
            virtual void convertTimeBase(const Time& t);

            //- Shift time
            virtual void shiftTimeBase(const Time& t);


        // Evaluation

            //- Return value as a function of (scalar) independent variable
            virtual Type value(const scalar x) const;

            //- Return value as a function of (scalar) independent variable
            virtual Type value(const vector& xyz) const;

            //- Return value as a function of (scalar) independent variable
            virtual tmp<Field<Type>> value(const scalarField& x) const;

            //- Return value as a function of (scalar) independent variable
            virtual tmp<Field<Type>> value(const Field<vector>& xyz) const;

            //- Integrate between two (scalar) values
            virtual Type integrate(const scalar x1, const scalar x2) const;

            //- Integrate between two (scalar) values
            virtual tmp<Field<Type>> integrate
            (
                const scalarField& x1,
                const scalarField& x2
            ) const;

            //- Integrate function value divided by abscissa between two
            // (scalar) values
            virtual Type integrateYoverX
            (
                const scalar x1,
                const scalar x2
            ) const;

            //- Integrate function value divided by sbscissa between two
            // (scalar) values
            virtual tmp<Field<Type>> integrateYoverX
            (
                const scalarField& x1,
                const scalarField& x2
            ) const;

            //- Derivative
            virtual Type derivative(const scalar x) const;

            //- Derivative
            virtual tmp<Field<Type>> derivative(const scalarField& x) const;


        // I/O

            //- Ostream Operator
            friend Ostream& operator<< <Type>
            (
                Ostream& os,
                const Function1<Type>& func
            );

            //- Write in dictionary format
            virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeFunction1(Type)                                                    \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(Function1<Type>, 0);                   \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        Function1<Type>,                                                       \
        dictionary                                                             \
    );


#define makeFunction1Type(SS, Type)                                            \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(Function1Types::SS<Type>, 0);          \
                                                                               \
    Function1<Type>::adddictionaryConstructorToTable<Function1Types::SS<Type>> \
        add##SS##Type##ConstructorToTable_;


#define makeScalarFunction1(SS)                                                \
                                                                               \
    defineTypeNameAndDebug(SS, 0);                                             \
                                                                               \
    Function1<scalar>::adddictionaryConstructorToTable<SS>                     \
        add##SS##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "primitives/functions/Function1/Function1/Function1.C"
    #include "primitives/functions/Function1/Constant/Constant.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
