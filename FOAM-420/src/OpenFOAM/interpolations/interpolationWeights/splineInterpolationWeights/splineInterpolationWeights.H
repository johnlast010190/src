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
    (c) 2012-2015 OpenFOAM Foundation

Class
    Foam::splineInterpolationWeights

Description
    Catmull-Rom spline interpolation.

SourceFiles
    splineInterpolationWeights.C

\*---------------------------------------------------------------------------*/

#ifndef splineInterpolationWeights_H
#define splineInterpolationWeights_H

#include "interpolations/interpolationWeights/interpolationWeights/interpolationWeights.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class splineInterpolationWeights Declaration
\*---------------------------------------------------------------------------*/

class splineInterpolationWeights
:
    public interpolationWeights
{

private:

    // Private data

        //- Cached index in samples from previous invocation
        mutable label index_;

public:

    //- Runtime type information
    TypeName("spline");

    // Constructors

        //- Construct from components. By default make sure samples are
        //  equidistant.
        splineInterpolationWeights
        (
            const scalarField& samples,
            const bool checkEqualDistance = true
        );


    //- Destructor
    virtual ~splineInterpolationWeights()
    {}


    // Member Functions

        //- Calculate weights and indices to calculate t from samples.
        //  Returns true if indices changed.
        virtual bool valueWeights
        (
            const scalar t,
            labelList& indices,
            scalarField& weights
        ) const;

        //- Calculate weights and indices to calculate integrand of t1..t2
        //  from samples. Returns true if indices changed.
        virtual bool integrationWeights
        (
            const scalar t1,
            const scalar t2,
            labelList& indices,
            scalarField& weights
        ) const
        {
            NotImplemented;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
