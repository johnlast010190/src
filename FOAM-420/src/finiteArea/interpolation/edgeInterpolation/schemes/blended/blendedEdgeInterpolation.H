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
    (c) 1991-2005 OpenCFD Ltd.

Class
    blendedEdgeInterpolation

Description
    linear/upwind blended differencing scheme.

SourceFiles
    blendedEdgeInterpolationMake.C

\*---------------------------------------------------------------------------*/

#ifndef blendedEdgeInterpolation_H
#define blendedEdgeInterpolation_H

#include "interpolation/edgeInterpolation/schemes/linear/linearEdgeInterpolation.H"
#include "interpolation/edgeInterpolation/schemes/upwind/upwindEdgeInterpolation.H"
#include "fields/areaFields/areaFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class blendedEdgeInterpolation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class blendedEdgeInterpolation
:
    public linearEdgeInterpolation<Type>,
    public upwindEdgeInterpolation<Type>
{
    // Private data

        const scalar blendingFactor_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        blendedEdgeInterpolation(const blendedEdgeInterpolation&);

        //- Disallow default bitwise assignment
        void operator=(const blendedEdgeInterpolation&);


public:

    //- Runtime type information
    TypeName("blended");


    // Constructors

        //- Construct from mesh, faceFlux and blendingFactor
        blendedEdgeInterpolation
        (
            const faMesh& mesh,
            const edgeScalarField& faceFlux,
            const scalar blendingFactor
        )
        :
            edgeInterpolationScheme<Type>(mesh),
            linearEdgeInterpolation<Type>(mesh),
            upwindEdgeInterpolation<Type>(mesh, faceFlux),
            blendingFactor_(blendingFactor)
        {}

        //- Construct from mesh and Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the database
        blendedEdgeInterpolation
        (
            const faMesh& mesh,
            Istream& is
        )
        :
            edgeInterpolationScheme<Type>(mesh),
            linearEdgeInterpolation<Type>(mesh),
            upwindEdgeInterpolation<Type>
            (
                mesh,
                mesh().objectRegistry::lookupObject<edgeScalarField>
                (
                    word(is)
                )
            ),
            blendingFactor_(readScalar(is))
        {}

        //- Construct from mesh, faceFlux and Istream
        blendedEdgeInterpolation
        (
            const faMesh& mesh,
            const edgeScalarField& faceFlux,
            Istream& is
        )
        :
            edgeInterpolationScheme<Type>(mesh),
            linearEdgeInterpolation<Type>(mesh),
            upwindEdgeInterpolation<Type>(mesh, faceFlux),
            blendingFactor_(readScalar(is))
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<edgeScalarField> weights
        (
            const GeometricField<Type, faPatchField, areaMesh>& vf
        ) const
        {
            return
                blendingFactor_*
                linearEdgeInterpolation<Type>::weights(vf)
              + (1 - blendingFactor_)*
                upwindEdgeInterpolation<Type>::weights(vf);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
