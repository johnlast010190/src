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
    (c) 2023 Esi Ltd

Class
    Foam::fv::correctedSnGrad

Group
    grpFvSnGradSchemes

Description
    Simple central-difference snGrad scheme with non-orthogonal correction.

SourceFiles
    correctedSnGrad.C

\*---------------------------------------------------------------------------*/

#ifndef correctedSnGrad_H
#define correctedSnGrad_H

#include "finiteVolume/snGradSchemes/snGradScheme/snGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class correctedSnGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class correctedSnGrad
:
    public snGradScheme<Type>
{

    // Private data

        const word& defaultGradSchemeName_;

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const correctedSnGrad&);


public:

    //- Runtime type information
    TypeName("corrected");


    // Constructors

        //- Construct from mesh
        correctedSnGrad
        (
            const fvMesh& mesh, const word& defaultGradSchemeName=word::null
        )
        :
            snGradScheme<Type>(mesh),
            defaultGradSchemeName_(defaultGradSchemeName)
        {}


        //- Construct from mesh and data stream
        correctedSnGrad
        (
            const fvMesh& mesh,
            Istream&,
            const word& defaultGradSchemeName
        )
        :
            snGradScheme<Type>(mesh),
            defaultGradSchemeName_(defaultGradSchemeName)
        {}


    //- Destructor
    virtual ~correctedSnGrad();


    // Member Functions

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> deltaCoeffs
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            return this->mesh().nonOrthDeltaCoeffs();
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the correctedSnGrad
        //  for the given field using the gradient of the field
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        fullGradCorrection
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Return the explicit correction to the correctedSnGrad
        //  for the given field using the gradients of the field components
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correction(const GeometricField<Type, fvPatchField, volMesh>&) const;
};


// * * * * * * * * Template Member Function Specialisations  * * * * * * * * //

template<>
tmp<surfaceScalarField> correctedSnGrad<scalar>::correction
(
    const volScalarField& vsf
) const;


template<>
tmp<surfaceVectorField> correctedSnGrad<vector>::correction
(
    const volVectorField& vvf
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/snGradSchemes/correctedSnGrad/correctedSnGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
