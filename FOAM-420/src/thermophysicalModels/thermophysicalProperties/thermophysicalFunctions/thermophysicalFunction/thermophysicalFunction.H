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
    Foam::thermophysicalFunction

Description
    Abstract base class for thermo-physical functions

\*---------------------------------------------------------------------------*/

#ifndef thermophysicalFunction_H
#define thermophysicalFunction_H

#include "primitives/Scalar/scalar/scalar.H"
#include "db/IOstreams/IOstreams.H"
#include "db/typeInfo/typeInfo.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "db/dictionary/dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class thermophysicalFunction Declaration
\*---------------------------------------------------------------------------*/

class thermophysicalFunction
{

public:

    //- Runtime type information
    TypeName("thermophysicalFunction");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            thermophysicalFunction,
            Istream,
            (Istream& is),
            (is)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            thermophysicalFunction,
            dictionary,
            (const dictionary& dict),
            (dict)
        );


    // Constructors

        //- Construct null
        thermophysicalFunction()
        {}

        //- Return pointer to new thermophysicalFunction created from input
        static autoPtr<thermophysicalFunction> New(Istream& is);

        //- Return pointer to new thermophysicalFunction created from dict
        static autoPtr<thermophysicalFunction> New(const dictionary& dict);


    //- Destructor
    virtual ~thermophysicalFunction()
    {}


    // Member Functions

        //- Evaluate the function and return the result
        virtual scalar f(scalar p, scalar T) const = 0;

        //- Write the function coefficients
        virtual void writeData(Ostream& os) const = 0;


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const thermophysicalFunction& f)
        {
            f.writeData(os);
            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
