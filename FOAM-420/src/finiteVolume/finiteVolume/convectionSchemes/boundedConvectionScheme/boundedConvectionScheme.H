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
    (c) 2012-2016 OpenFOAM Foundation
    (c) 2019-2022 Esi Ltd.

Class
    Foam::fv::boundedConvectionScheme

Group
    grpFvConvectionSchemes

Description
    Bounded form of the selected convection scheme.

    Boundedness is achieved by subtracting div(phi)*vf or Sp(div(phi), vf)
    which is non-conservative if div(phi) != 0 but conservative otherwise.

    Can be used for convection of bounded scalar properties in steady-state
    solvers to improve stability if insufficient convergence of the pressure
    equation causes temporary divergence of the flux field.

    In the case of a moving mesh:
    1. If phi has dimensions of volume/time, it is converted to absolute flux
    2. Otherwise, relative flux is used and 'bounded' should ordinarily be used
       on the corresponding ddt scheme in order to be consistent.

SourceFiles
    boundedConvectionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef boundedConvectionScheme_H
#define boundedConvectionScheme_H

#include "finiteVolume/convectionSchemes/convectionScheme/convectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class boundedConvectionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class boundedConvectionScheme
:
    public fv::convectionScheme<Type>
{
    // Private data

        tmp<fv::convectionScheme<Type>> scheme_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        boundedConvectionScheme(const boundedConvectionScheme&);

        //- Disallow default bitwise assignment
        void operator=(const boundedConvectionScheme&);

        //- Return absolute flux in moving mesh case
        tmp<surfaceScalarField> absoluteFlux
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;

public:

    //- Runtime type information
    TypeName("bounded");


    // Constructors
    //
        //- Construct from flux and Istream
        boundedConvectionScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            const tmp<fv::convectionScheme<Type>> scheme
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            scheme_
            (
                scheme
            )
        {}


        //- Construct from flux and Istream
        boundedConvectionScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            scheme_
            (
                fv::convectionScheme<Type>::New(mesh, faceFlux, is)
            )
        {}

        //- Construct multivariate from flux and Istream
        boundedConvectionScheme
        (
            const fvMesh& mesh,
            const typename multivariateSurfaceInterpolationScheme<Type>::
                fieldTable& fields,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            scheme_
            (
                fv::convectionScheme<Type>::New(mesh, fields, faceFlux, is)
            )
        {}


    // Member Functions

        const convectionScheme<Type>& scheme() const;

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> interpolate
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> flux
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<fvMatrix<Type>> fvmDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<fvMatrix<Type>> fvmDiv
        (
            const surfaceScalarField&,
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        virtual tmp<fvBlockMatrix<Type>> fvmBDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
