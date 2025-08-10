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
    (c) 2019-2023 Esi Ltd.

Class
    Foam::Function1Types::NSRDSfunc4

Description
    Function1 version of NSRDSfunc4 given by
    a_ + b_/x + c_/pow(x, 3) + d_/pow(x, 8) + e_/pow(x, 9) .

    Evaluates the function with the supplied stream of coefficients
    and return the result, the value of the definite integral, or the value of
    the definite integral of the function divided by x, or the derivative.

    NSRDSfunc4 is NSRDS function number 104

    Source:
    \verbatim
                      NSRDS - AICHE
                 Data Compilation Tables
                    of Properties of
                     Pure Compounds

        Design Institute for Physical Property Data
          American Institute of Chemical Engineers
                  345 East 47th Street
                New York, New York 10017

         National Standard Reference Data System
         American Institute of Chemical Engineers

          T.E. Daubert       -       R.P. Danner

            Department of Chemical Engineering
            The Pennsylvania State University
                University Park, PA 16802
    \endverbatim

SourceFiles
    NSRDSfunc4.C

\*---------------------------------------------------------------------------*/

#ifndef NSRDSfunc4_H
#define NSRDSfunc4_H

#include "primitives/functions/Function1/Function1/Function1Fwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                        Class NSRDSfunc4 Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class NSRDSfunc4
:
    public Function1<Type>
{
    // Private data

        //- Equation coefficients
        FixedList<Type, 5> coeffs_;

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const NSRDSfunc4<Type>&);

        //- Return function value of the anti-derivative / indefinite integral
        Type indefiniteIntegral(const scalar x) const;

        //- Return function value of the anti-derivative / indefinite integral
        //  of the function divided by x
        Type indefiniteIntegralYoverX(const scalar x) const;



public:

    //- Runtime type information
    TypeName("NSRDSfunc4");


    // Constructors

        NSRDSfunc4(const word& entryName, const dictionary& dict);

        //- Construct from components
        NSRDSfunc4
        (
            const word& entryName,
            const FixedList<Type, 5>&
        );

        //- Copy constructor
        NSRDSfunc4(const NSRDSfunc4& poly);

        //- Construct and return a clone
        virtual tmp<Function1<Type>> clone() const
        {
            return tmp<Function1<Type>>(new NSRDSfunc4(*this));
        }


    //- Destructor
    virtual ~NSRDSfunc4();


    // Member Functions

        // Manipulation

            //- Convert time
            virtual void convertTimeBase(const Time& t);


        // Evaluation

            //- Return function value
            Type value(const scalar x) const;

            //- Integrate between two (scalar) values
            Type integrate(const scalar x1, const scalar x2) const;

            //- Integrate function value divided by abscissa between two
            // (scalar) values
            Type integrateYoverX(const scalar x1, const scalar x2) const;

            //- Derivative at the point
            Type derivative(const scalar x) const;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "primitives/functions/Function1/NSRDSfunc4/NSRDSfunc4.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
