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
    Foam::NoComposition

Group
    grpLagrangianIntermediateCompositionSubModels

Description
    Dummy class for 'none' option - will raise an error if any functions are
    called that require return values.

SourceFiles
    NoComposition.C

\*---------------------------------------------------------------------------*/

#ifndef NoComposition_H
#define NoComposition_H

#include "submodels/Reacting/CompositionModel/CompositionModel/CompositionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class NoComposition Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoComposition
:
    public CompositionModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from dictionary
        NoComposition(const dictionary& dict, CloudType& owner);

        //- Construct copy
        NoComposition(const NoComposition<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<CompositionModel<CloudType>> clone() const
        {
            return autoPtr<CompositionModel<CloudType>>
            (
                new NoComposition<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoComposition();


    // Member Functions

        // Access

            // Mixture properties

                //- Return the list of mixture mass fractions
                //  If only 1 phase, return component fractions of that phase
                virtual const scalarField& YMixture0() const;

                // Indices of gas, liquid and solid phases in phase properties
                // list

                    //- Gas id
                    virtual label idGas() const;

                    //- Liquid id
                    virtual label idLiquid() const;

                    //- Solid id
                    virtual label idSolid() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "submodels/Reacting/CompositionModel/NoComposition/NoComposition.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
