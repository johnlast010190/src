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
    (c) 2022 Esi Ltd.

Class
    upwindEdgeInterpolation

Description
    Upwind differencing scheme class.

SourceFiles
    upwindEdgeInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef upwindEdgeInterpolation_H
#define upwindEdgeInterpolation_H

#include "interpolation/edgeInterpolation/edgeInterpolationScheme/edgeInterpolationScheme.H"
#include "fields/areaFields/areaFields.H"
#include "fields/edgeFields/edgeFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class upwind Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class upwindEdgeInterpolation
:
    virtual public edgeInterpolationScheme<Type>
{
    // Private data

        const edgeScalarField& faceFlux_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const upwindEdgeInterpolation&);


public:

    //- Runtime type information
    TypeName("upwind");


    // Constructors

        //- Construct from faceFlux
        upwindEdgeInterpolation
        (
            const faMesh& mesh,
            const edgeScalarField& faceFlux
        )
        :
            edgeInterpolationScheme<Type>(mesh),
            faceFlux_(faceFlux)
        {}

        //- Construct from Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the database
        upwindEdgeInterpolation
        (
            const faMesh& mesh,
            Istream& is
        )
        :
            edgeInterpolationScheme<Type>(mesh),
            faceFlux_
            (
                mesh().objectRegistry::lookupObject<edgeScalarField>
                (
                    word(is)
                )
            )
        {}

        //- Construct from faceFlux and Istream
        upwindEdgeInterpolation
        (
            const faMesh& mesh,
            const edgeScalarField& faceFlux,
            Istream&
        )
        :
            edgeInterpolationScheme<Type>(mesh),
            faceFlux_(faceFlux)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<edgeScalarField> weights
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const
        {
            return pos0(faceFlux_);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
