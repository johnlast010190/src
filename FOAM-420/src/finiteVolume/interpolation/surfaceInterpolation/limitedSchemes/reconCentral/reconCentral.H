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
    (c) 2017 Wikki Ltd

Class
    Foam::reconCentral

Description
    reconCentral interpolation scheme uses extrapolated gradient-based
    correction from both sides onto the face, using 1/2 weighting.
    Class is derived from upwind and returns upwind weighting factors.
    Correction is calculated using the extrapolated gradients.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    reconCentral.C

\*---------------------------------------------------------------------------*/

#ifndef reconCentral_H
#define reconCentral_H

#include "interpolation/surfaceInterpolation/limitedSchemes/upwind/upwind.H"
#include "finiteVolume/gradSchemes/gaussGrad/gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class reconCentral Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class reconCentral
:
    public upwind<Type>
{
    // Private Data

        tmp<fv::gradScheme<Type>> gradScheme_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        reconCentral(const reconCentral&);

        //- Disallow default bitwise assignment
        void operator=(const reconCentral&);


public:

    //- Runtime type information
    TypeName("reconCentral");


    // Constructors

        //- Construct from faceFlux
        reconCentral
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux
        )
        :
            upwind<Type>(mesh, faceFlux),
            gradScheme_
            (
                new fv::gaussGrad<Type>(mesh)
            )
        {}

        //- Construct from Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        reconCentral
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& schemeData
        )
        :
            upwind<Type>(mesh, db, schemeData),
            gradScheme_
            (
                fv::gradScheme<Type>::New
                (
                    mesh,
                    db,
                    schemeData
                )
            )
        {}

        //- Construct from faceFlux and Istream
        reconCentral
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        )
        :
            upwind<Type>(mesh, faceFlux, schemeData),
            gradScheme_
            (
                fv::gradScheme<Type>::New
                (
                    mesh, faceFlux.db(), schemeData
                )
            )
        {}


    // Member Functions

        //- Return the face-interpolate of the given cell field
        //  with explicit correction
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        interpolate(const GeometricField<Type, fvPatchField, volMesh>&) const;

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
