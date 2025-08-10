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
    Foam::energyScalingFunction

Description

SourceFiles
    energyScalingFunction.C
    energyScalingFunctionNew.C

\*---------------------------------------------------------------------------*/

#ifndef energyScalingFunction_H
#define energyScalingFunction_H

#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "db/typeInfo/typeInfo.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "memory/autoPtr/autoPtr.H"
#include "pairPotential/basic/pairPotential.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class energyScalingFunction Declaration
\*---------------------------------------------------------------------------*/

class energyScalingFunction
{

protected:

    // Protected data

        word name_;

        dictionary energyScalingFunctionProperties_;

        const pairPotential& pairPot_;


    // Private Member Functions

        //- Disallow copy construct
        energyScalingFunction(const energyScalingFunction&);

        //- Disallow default bitwise assignment
        void operator=(const energyScalingFunction&);


public:

    //- Runtime type information
    TypeName("energyScalingFunction");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            energyScalingFunction,
            dictionary,
            (
                const word& name,
                const dictionary& energyScalingFunctionProperties,
                const pairPotential& pairPot
            ),
            (name, energyScalingFunctionProperties, pairPot)
        );


    // Selectors

        //- Return a reference to the selected viscosity model
        static autoPtr<energyScalingFunction> New
        (
            const word& name,
            const dictionary& energyScalingFunctionProperties,
            const pairPotential& pairPot
        );


    // Constructors

        //- Construct from components
        energyScalingFunction
        (
            const word& name,
            const dictionary& energyScalingFunctionProperties,
            const pairPotential& pairPot
        );


    //- Destructor
    virtual ~energyScalingFunction()
    {}


    // Member Functions

        virtual void scaleEnergy(scalar& e, const scalar r) const = 0;

        const dictionary& energyScalingFunctionProperties() const
        {
            return energyScalingFunctionProperties_;
        }

        //- Read energyScalingFunction dictionary
        virtual bool read
        (
            const dictionary& energyScalingFunctionProperties
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
