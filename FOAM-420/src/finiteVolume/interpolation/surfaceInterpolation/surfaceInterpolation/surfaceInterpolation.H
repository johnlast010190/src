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
    (c) 2010-2012 Esi Ltd.
    (c) 2011 OpenFOAM Foundation

Class
    Foam::surfaceInterpolation

Description
    Cell to surface interpolation scheme. Included in fvMesh.

SourceFiles
    surfaceInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceInterpolation_H
#define surfaceInterpolation_H

#include "memory/tmp/tmp.H"
#include "primitives/Scalar/scalar/scalar.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "db/typeInfo/className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class surfaceInterpolation Declaration
\*---------------------------------------------------------------------------*/

class surfaceInterpolation
{
    // Private data

        // Reference to fvMesh
        const fvMesh& mesh_;

        // Demand-driven data

            //- Linear difference weighting factors
            mutable surfaceScalarField* weights_;

            //- Cell-centre difference coefficients
            mutable surfaceScalarField* deltaCoeffs_;

            //- Non-orthogonal cell-centre difference coefficients
            mutable surfaceScalarField* nonOrthDeltaCoeffs_;

            //- Non-orthogonality correction vectors
            mutable surfaceVectorField* nonOrthCorrectionVectors_;


    // Private Member Functions

        //- Construct central-differencing weighting factors
        void makeWeights() const;

        //- Construct face-gradient difference factors
        void makeDeltaCoeffs() const;

        //- Construct face-gradient difference factors
        void makeNonOrthDeltaCoeffs() const;

        //- Construct non-orthogonality correction vectors
        void makeNonOrthCorrectionVectors() const;


protected:

    // Protected Member Functions

        // Storage management

            //- Clear all geometry and addressing
            void clearOut();


public:

    // Declare name of the class and its debug switch
    ClassName("surfaceInterpolation");


    // Constructors

        //- Construct given an fvMesh
        explicit surfaceInterpolation(const fvMesh&);


    //- Destructor
    virtual ~surfaceInterpolation();


    // Member functions

        //- Return reference to linear difference weighting factors
        virtual const surfaceScalarField& weights() const;

        //- Return reference to cell-centre difference coefficients
        virtual const surfaceScalarField& deltaCoeffs() const;

        //- Return reference to non-orthogonal cell-centre difference
        //  coefficients
        const surfaceScalarField& nonOrthDeltaCoeffs() const;

        //- Return reference to non-orthogonality correction vectors
        const surfaceVectorField& nonOrthCorrectionVectors() const;

        //- Do what is neccessary if the mesh has moved
        bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
