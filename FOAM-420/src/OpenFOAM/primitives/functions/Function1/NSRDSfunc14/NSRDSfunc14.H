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
    Foam::Function1Types::NSRDSfunc14

Description
    Function1 version of NSRDSfunc14 given by
       a_*a_/(t + ROOTVSMALL) + b_ - t
      *(
          2.0*a_*c_
          + t*(a_*d_ + t*(c_*c_/3.0 + t*(0.5*c_*d_ + 0.2*d_*d_*t)))
      )

     where
     scalar t = 1.0 - Tdash/Tc_ , and
     Tdash = min(T, Tc_ - ROOTVSMALL) .

    Inputs are specified as
    NSRDSfunc6 Tc_ ( a_ b_ c_ d_ );

    Evaluates the function with the supplied stream of coefficients
    and return the result, or the value of the derivative.

    NSRDSfunc14 is NSRDS function number 114

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
    NSRDSfunc14.C

\*---------------------------------------------------------------------------*/

#ifndef NSRDSfunc14_H
#define NSRDSfunc14_H

#include "primitives/functions/Function1/Function1/Function1Fwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                        Class NSRDSfunc14 Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class NSRDSfunc14
:
    public Function1<Type>
{
    // Private data

        //- Reference temperature
        scalar Tc_;

        //- Equation coefficients
        FixedList<Type, 4> coeffs_;

    // Private Member Functions


public:

    //- Runtime type information
    TypeName("NSRDSfunc14");


    // Constructors

        NSRDSfunc14(const word& entryName, const dictionary& dict);

        //- Construct from components
        NSRDSfunc14
        (
            const word& entryName,
            const scalar&,
            const FixedList<Type, 4>&
        );

        //- Copy constructor
        NSRDSfunc14(const NSRDSfunc14& poly);

        //- Construct and return a clone
        virtual tmp<Function1<Type>> clone() const
        {
            return tmp<Function1<Type>>(new NSRDSfunc14(*this));
        }


    //- Destructor
    virtual ~NSRDSfunc14();


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
    #include "primitives/functions/Function1/NSRDSfunc14/NSRDSfunc14.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
