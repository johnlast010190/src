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
    Foam::Function1Types::NSRDSfunc7

Description
    Function1 version of NSRDSfunc7 given by
    a_ + b_*sqr((c_/x)/sinh(c_/x)) + d_*sqr((e_/x)/cosh(e_/x)).

    Evaluates the function with the supplied stream of coefficients
    and return the result, or the value of the derivative.

    NSRDSfunc7 is NSRDS function number 107

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
    NSRDSfunc7.C

\*---------------------------------------------------------------------------*/

#ifndef NSRDSfunc7_H
#define NSRDSfunc7_H

#include "primitives/functions/Function1/Function1/Function1Fwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                        Class NSRDSfunc7 Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class NSRDSfunc7
:
    public Function1<Type>
{
    // Private data

        //- Equation coefficients
        FixedList<Type, 5> coeffs_;

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const NSRDSfunc7<Type>&);

        //- Return the value of the sinh() term
        Type sinhPart(const scalar x) const;

        //- Return the value of the cosh() term
        Type coshPart(const scalar x) const;


public:

    //- Runtime type information
    TypeName("NSRDSfunc7");


    // Constructors

        NSRDSfunc7(const word& entryName, const dictionary& dict);

        //- Construct from components
        NSRDSfunc7
        (
            const word& entryName,
            const FixedList<Type, 5>&
        );

        //- Copy constructor
        NSRDSfunc7(const NSRDSfunc7& poly);

        //- Construct and return a clone
        virtual tmp<Function1<Type>> clone() const
        {
            return tmp<Function1<Type>>(new NSRDSfunc7(*this));
        }


    //- Destructor
    virtual ~NSRDSfunc7();


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
    #include "primitives/functions/Function1/NSRDSfunc7/NSRDSfunc7.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
