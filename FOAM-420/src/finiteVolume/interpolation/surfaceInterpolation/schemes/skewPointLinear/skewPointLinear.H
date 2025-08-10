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
    (c) 2017 Esi Ltd.

Class
    skewPointLinear

Description
    Improved version of pointLinear to remove incorrect dependences on
    linear interpolated face centre values.

    Face-point interpolation scheme class derived from linear and
    returns linear weighting factors but also applies an explicit correction.

    Uses volPointInterpolation to obtain the field values at the face-points.


SourceFiles
    skewPointLinear.C

\*---------------------------------------------------------------------------*/

#ifndef skewPointLinear_H
#define skewPointLinear_H

#include "interpolation/surfaceInterpolation/schemes/linear/linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class skewPointLinear Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class skewPointLinear
:
    public linear<Type>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        skewPointLinear(const skewPointLinear&);

        //- Disallow default bitwise assignment
        void operator=(const skewPointLinear&);


public:

    //- Runtime type information
    TypeName("skewPointLinear");


    // Constructors

        //- Construct from mesh
        skewPointLinear(const fvMesh& mesh)
        :
            linear<Type>(mesh)
        {}


        //- Construct from mesh and Istream
        skewPointLinear
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream&
        )
        :
            linear<Type>(mesh)
        {}


        //- Construct from mesh, faceFlux and Istream
        skewPointLinear
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream&
        )
        :
            linear<Type>(mesh)
        {}


    // Member Functions

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
