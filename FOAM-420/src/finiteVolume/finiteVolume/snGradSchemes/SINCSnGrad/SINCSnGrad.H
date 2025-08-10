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

Class
    Foam::fv::SINCSnGrad

SourceFiles
    SINCSnGrad.C

Authors
    Simon Hill         <simon.hill@linde.com>
    Daniel Deising     <deising@mma.tu-darmstadt.de>
    All rights reserved.

Description
    Simple central-difference scalar-type snGrad scheme with implicit
    bounded non-orthogonal correction.

\*---------------------------------------------------------------------------*/

#ifndef SINCSnGrad_H
#define SINCSnGrad_H

#include "finiteVolume/snGradSchemes/snGradScheme/snGradScheme.H"
#include "interpolation/surfaceInterpolation/schemes/linear/linear.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "finiteVolume/gradSchemes/gaussGrad/gaussGrad.H"
#include "finiteVolume/fvc/fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class SINCSnGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class SINCSnGrad
:
    public snGradScheme<Type>
{
    // Private data

        const word& defaultGradSchemeName_;

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const SINCSnGrad&);


public:

    //- Runtime type information
    TypeName("SINC");


    // Constructors

        //- Construct from mesh
        SINCSnGrad(const fvMesh& mesh)
        :
            snGradScheme<Type>(mesh)
        {}


        //- Construct from mesh and data stream
        SINCSnGrad
        (
            const fvMesh& mesh, Istream&, const word& defaultGradSchemeName
        )
        :
            snGradScheme<Type>(mesh),
            defaultGradSchemeName_(defaultGradSchemeName)
        {}


    //- Destructor
    virtual ~SINCSnGrad();


    // Member Functions

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> deltaCoeffs
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return false;
        }

        //- Return the explicit correction to the orthogonalSnGrad
        //  for the given field
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correction(const GeometricField<Type, fvPatchField, volMesh>&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "finiteVolume/snGradSchemes/SINCSnGrad/SINCSnGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
