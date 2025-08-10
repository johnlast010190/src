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
    (c) 2020 Esi Ltd.

Class
    Foam::BRICS

Group
    grpFvLimitedSurfaceInterpolationSchemes

Description
    BRICS differencing scheme class which returns the weights based on the
    normalised variable formulation.

SourceFiles
    BRICS.C

\*---------------------------------------------------------------------------*/

#ifndef BRICS_H
#define BRICS_H

#include "interpolation/surfaceInterpolation/limitedSchemes/limitedSurfaceInterpolationScheme/limitedSurfaceInterpolationScheme.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class BRICS Declaration
\*---------------------------------------------------------------------------*/

class BRICS
:
    public limitedSurfaceInterpolationScheme<scalar>
{
    protected:

        //- input data dictionary
        dictionary schemeData_;

        //- alpha field name
        mutable word blendingField_;

        //- model parameter: standard settings
        const scalar beta_GDS_ = 0.1;
        const scalar beta_IGDS_ = 0.5;
        const scalar p_GDS_ = 0.5;
        const scalar p_IGDS_ = 0.0;

        //- to avoid division by zero when calculating interface normal vector
        const dimensionedScalar deltaN_;

        //- blending of compressive and shape-preserving higher order scheme
        //  based on interface/cell-face orientation
        //- blending only in a small region around the interface to achieve
        //  consistent advection for volumetric phase fraction and species
        //  concentration fields
        mutable surfaceScalarField blendingFactor_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const BRICS&);


        //- Calculate individual limiter
        scalar limiter
        (
            const scalar cdWeight,
            const scalar faceFlux,
            const scalar& phiP,
            const scalar& phiN,
            const vector& gradcP,
            const vector& gradcN,
            const scalar CourantF,
            const vector d,
            const scalar blendingFactor
        ) const;

        //- Calculate individual limiter
        virtual scalar weight
        (
            const scalar faceFlux,
            const scalar& phiP,
            const scalar& phiN,
            const vector& gradcP,
            const vector& gradcN,
            const scalar CourantF,
            const vector d,
            const scalar blendingFactor
        ) const;


public:

    //- Runtime type information
    TypeName("BRICS");


    // Constructors

        //- Construct from faceFlux
        BRICS
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux
        )
        :
            limitedSurfaceInterpolationScheme<scalar>(mesh, faceFlux),
            schemeData_(),
            blendingField_("alpha.water"),
            deltaN_
            (
                "deltaN",
                1e-8/pow(average(this->mesh().V()), 1.0/3.0)
            ),
            blendingFactor_
            (
                IOobject
                (
                    "blendingFactorICS",
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh(),
                dimensionedScalar("blendingFactorICS", dimless, scalar(0.0))
            )
        {}

        //- Construct from Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        BRICS
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& is
        )
        :
            limitedSurfaceInterpolationScheme<scalar>(mesh, db, is),
            schemeData_(),
            blendingField_("alpha.water"),
            deltaN_
            (
                "deltaN",
                1e-8/pow(average(this->mesh().V()), 1.0/3.0)
            ),
            blendingFactor_
            (
                IOobject
                (
                    "blendingFactorICS",
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh(),
                dimensionedScalar("blendingFactorICS", dimless, scalar(0.0))
            )
        {
            if (schemeData_.found("blendingField"))
            {
                const word toRead(schemeData_.lookup("blendingField"));
                blendingField_ = toRead;
            }
        }

        //- Construct from faceFlux and Istream
        BRICS
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream&
        )
        :
            limitedSurfaceInterpolationScheme<scalar>(mesh, faceFlux),
            schemeData_(),
            blendingField_("alpha.water"),
            deltaN_
            (
                "deltaN",
                1e-8/pow(average(this->mesh().V()), 1.0/3.0)
            ),
            blendingFactor_
            (
                IOobject
                (
                    "blendingFactorICS",
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh(),
                dimensionedScalar("blendingFactorICS", dimless, scalar(0.0))
            )
        {
            if (schemeData_.found("blendingField"))
            {
                const word toRead(schemeData_.lookup("blendingField"));
                blendingField_ = toRead;
            }
        }


    // Member Functions

        //- Return the interpolation limiter
        virtual tmp<surfaceScalarField> limiter
        (
            const volScalarField&
        ) const;

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> calcWeights
        (
            const volScalarField&
        ) const;

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> weights
        (
            const volScalarField& vf
        ) const
        {
            return calcWeights(vf);
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return false;
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<surfaceScalarField> correction
        (
            const volScalarField&
        ) const;

        //- Calculate cell-face/interface orientation
        void calculateTheta() const;

        //- Calculate the scheme blending factor
        virtual void calculateBlendingFactor() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
