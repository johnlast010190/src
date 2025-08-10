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

Class
    Foam::distributionModel

Description
    A library of runtime-selectable distribution models.

    Returns a sampled value given the expectation (nu) and variance (sigma^2)

    Current distribution models include:
    - exponential
    - fixedValue
    - general
    - multi-normal
    - normal
    - Rosin-Rammler
    - Mass-based Rosin-Rammler
    - uniform

    The distributionModel is tabulated in equidistant nPoints, in an interval.
    These values are integrated to obtain the cumulated distribution model,
    which is then used to change the distribution from unifrom to
    the actual distributionModel.

SourceFiles
    distributionModel.C
    distributionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef distributionModel_H
#define distributionModel_H

#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "memory/autoPtr/autoPtr.H"
#include "primitives/random/Random/Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class distributionModel Declaration
\*---------------------------------------------------------------------------*/

class distributionModel
{

protected:

    // Protected data

        //- Coefficients dictionary
        const dictionary distributionModelDict_;

        //- Reference to the random number generator
        Random& rndGen_;


    // Protected Member Functions

        //- Check that the distribution model is valid
        virtual void check() const;


public:

    //-Runtime type information
    TypeName("distributionModel");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        distributionModel,
        dictionary,
        (
            const dictionary& dict,
            Random& rndGen
        ),
        (dict, rndGen)
    );


    // Constructors

        //- Construct from dictionary
        distributionModel
        (
            const word& name,
            const dictionary& dict,
            Random& rndGen
        );

        //- Construct copy
        distributionModel(const distributionModel& p);

        //- Construct and return a clone
        virtual autoPtr<distributionModel> clone() const = 0;


    //- Selector
    static autoPtr<distributionModel> New
    (
        const dictionary& dict,
        Random& rndGen
    );


    //- Destructor
    virtual ~distributionModel();


    // Member Functions

        //- Sample the distributionModel
        virtual scalar sample() const = 0;

        //- Return the minimum value
        virtual scalar minValue() const = 0;

        //- Return the maximum value
        virtual scalar maxValue() const = 0;

        //- Return the maximum value
        virtual scalar meanValue() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
