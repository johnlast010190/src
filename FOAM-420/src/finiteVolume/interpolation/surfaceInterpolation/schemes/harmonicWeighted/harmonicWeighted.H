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
    (c) 2017 Esi Ltd.

Class
    Foam::harmonicWeighted

Description
    Two-scheme harmonicWeighted differencing scheme.

SourceFiles
    harmonicWeighted.C

\*---------------------------------------------------------------------------*/

#ifndef harmonicWeighted_H
#define harmonicWeighted_H

#include "interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.H"
#include "interpolation/surfaceInterpolation/schemes/linear/linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class harmonicWeighted Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class harmonicWeighted
:
    public surfaceInterpolationScheme<Type>
{
    // Private Member Functions

        //- Primary underlying scheme
        tmp<linear<Type>> tScheme1_;

        //- Weight interpolation
        tmp<linear<scalar>> tScheme2_;

        //- name of field used for weighting
        word weightName_;

        //- Disallow default bitwise copy construct
        harmonicWeighted(const harmonicWeighted&);

        //- Disallow default bitwise assignment
        void operator=(const harmonicWeighted&);


public:

    //- Runtime type information
    TypeName("harmonicWeighted");


    // Constructors

        //- Construct from mesh and Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        harmonicWeighted
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            tScheme1_
            (
                new linear<Type>(mesh)
            ),
            tScheme2_
            (
                new linear<scalar>(mesh)
            ),
            weightName_(is)
        {}

        //- Construct from mesh, faceFlux and Istream
        harmonicWeighted
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            tScheme1_
            (
                new linear<Type>(mesh)
            ),
            tScheme2_
            (
                new linear<scalar>(mesh)
            ),
            weightName_(is)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            NotImplemented;
        }

        //- Return the face-interpolate of the given cell field
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            const fvMesh& mesh = vf.mesh();

            const volScalarField& Cw
            (
                mesh.objectRegistry::template
                lookupObject<const volScalarField>(weightName_)
            );

            return
            (
                tScheme1_->interpolate(vf / Cw)
              / tScheme2_->interpolate(1.0/Cw)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
