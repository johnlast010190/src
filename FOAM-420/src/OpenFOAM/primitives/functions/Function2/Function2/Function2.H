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
    (c) 2020-2021 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::Function2

Description
    Run-time selectable function of two variables

    with many options provided from simple constant values to complex
    functions, interpolated tabulated data etc. etc.

SourceFiles
    Function2.C
    Function2New.C

\*---------------------------------------------------------------------------*/

#ifndef Function2_H
#define Function2_H

#include "db/dictionary/dictionary.H"
#include "memory/tmp/tmp.H"
#include "db/typeInfo/typeInfo.H"
#include "fields/Fields/Field/Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
template<class Type> class Function2;
template<class Type> Ostream& operator<<(Ostream&, const Function2<Type>&);

/*---------------------------------------------------------------------------*\
                          Class Function2 Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Function2
:
    public tmp<Function2<Type>>::refCount
{

protected:

    // Protected data

        //- Name of entry
        const word name_;


public:

    typedef Type returnType;

    //- Runtime type information
    TypeName("Function2");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        Function2,
        dictionary,
        (
            const word& name,
            const dictionary& dict
        ),
        (name, dict)
    );


    // Constructors

        //- Construct from name
        Function2(const word& name);

        //- Copy constructor
        Function2(const Function2<Type>& f1);

        //- Construct and return a clone
        virtual tmp<Function2<Type>> clone() const
        {
            return tmp<Function2<Type>>(new Function2<Type>(*this));
        }


    //- Selector
    static autoPtr<Function2<Type>> New
    (
        const word& name,
        const dictionary& dict
    );


    //- Destructor
    virtual ~Function2();


    // Member Functions

        //- Return the name of the entry
        const word& name() const;

        //- Return value as a function of two scalars
        virtual Type value(const scalar x, const scalar y) const
        {
            NotImplemented;
        }

        //- Return value as a function of two scalar fields
        virtual tmp<Field<Type>> value
        (
            const scalarField& x,
            const scalarField& y
        ) const
        {
            NotImplemented;
        }

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;


    // Member Operators

        //- Assignment
        void operator=(const Function2<Type>&);


    // IOstream Operators

        //- Ostream Operator
        friend Ostream& operator<< <Type>
        (
            Ostream& os,
            const Function2<Type>& func
        );
};


/*---------------------------------------------------------------------------*\
                       Class FieldFunction2 Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class Function2Type>
class FieldFunction2
:
    public Function2<Type>
{

public:

    // Constructors

        //- Construct from name
        FieldFunction2(const word& name);

        //- Construct and return a clone
        virtual tmp<Function2<Type>> clone() const;


    //- Destructor
    virtual ~FieldFunction2();


    // Member Functions

        //- Return value as a function of two scalars
        virtual Type value(const scalar x, const scalar y) const = 0;

        //- Return value as a function of two scalar fields
        virtual tmp<Field<Type>> value
        (
            const scalarField& x,
            const scalarField& y
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeFunction2(Type)                                                    \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(Function2<Type>, 0);                   \
                                                                               \
    defineTemplateRunTimeSelectionTable(Function2<Type>, dictionary);


#define makeFunction2Type(SS, Type)                                            \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(Function2s::SS<Type>, 0);              \
                                                                               \
    Function2<Type>::adddictionaryConstructorToTable<Function2s::SS<Type>>     \
        addFunction2##SS##Type##ConstructorToTable_;


#define makeScalarFunction2(SS)                                                \
                                                                               \
    defineTypeNameAndDebug(SS, 0);                                             \
                                                                               \
    Function2<scalar>::adddictionaryConstructorToTable<SS>                     \
        addFunction2##SS##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "primitives/functions/Function2/Function2/Function2.C"
    #include "primitives/functions/Function2/Constant/Constant2.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
