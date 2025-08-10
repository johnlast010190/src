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
    (c) 2017 Esi Ltd.

Class
    Foam::skewCorrected

Group
    grpFvSurfaceInterpolationSchemes

Description
    Skewness-corrected interpolation scheme that applies an explicit
    correction to given scheme.

SourceFiles
    skewCorrected.C

\*---------------------------------------------------------------------------*/

#ifndef skewCorrected_H
#define skewCorrected_H

#include "interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.H"
#include "interpolation/surfaceInterpolation/schemes/skewCorrected/skewCorrectionVectors.H"
#include "interpolation/surfaceInterpolation/schemes/linear/linear.H"
#include "finiteVolume/gradSchemes/gaussGrad/gaussGrad.H"
#include "finiteVolume/gradSchemes/gradScheme/gradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class skewCorrected Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class skewCorrected
:
    public surfaceInterpolationScheme<Type>
{
    // Private member data

        tmp<surfaceInterpolationScheme<Type>> tScheme_;

        tmp<fv::gradScheme<typename pTraits<Type>::cmptType>> tgradScheme_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        skewCorrected(const skewCorrected&);

        //- Disallow default bitwise assignment
        void operator=(const skewCorrected&);


public:

    //- Runtime type information
    TypeName("skewCorrected");


    // Constructors

        //- Construct from mesh and Istream
        skewCorrected
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            tScheme_
            (
                surfaceInterpolationScheme<Type>::New(mesh, db, is)
            )
        {
            if (is.eof())
            {
                tgradScheme_ =
                    tmp<fv::gradScheme<typename pTraits<Type>::cmptType>>
                    (
                        new fv::gaussGrad<typename pTraits<Type>::cmptType>(mesh)
                    );
            }
            else
            {
                tgradScheme_ =
                    tmp<fv::gradScheme<typename pTraits<Type>::cmptType>>
                    (
                        fv::gradScheme<typename pTraits<Type>::cmptType>::New
                        (
                            mesh, db, is
                        )
                    );
            }
        }


        //- Construct from mesh, faceFlux and Istream
        skewCorrected
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            tScheme_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            )
        {
            if (is.eof())
            {
                tgradScheme_ =
                    tmp<fv::gradScheme<typename pTraits<Type>::cmptType>>
                    (
                        new fv::gaussGrad<typename pTraits<Type>::cmptType>
                        (
                            mesh
                        )
                    );
            }
            else
            {
                tgradScheme_ =
                    tmp<fv::gradScheme<typename pTraits<Type>::cmptType>>
                    (
                        fv::gradScheme<typename pTraits<Type>::cmptType>::New
                        (
                            mesh,
                            faceFlux.db(),
                            is
                        )
                    );
            }
        }


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            return tScheme_().weights(vf);
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return
                tScheme_().corrected()
             || skewCorrectionVectors::New(this->mesh()).skew();
        }

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        skewCorrection
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            const fvMesh& mesh = this->mesh();

            const skewCorrectionVectors& scv = skewCorrectionVectors::New(mesh);

            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsfCorr
            (
                new GeometricField<Type, fvsPatchField, surfaceMesh>
                (
                    IOobject
                    (
                        "skewCorrected::skewCorrection(" + vf.name() + ')',
                        mesh.time().timeName(),
                        mesh
                    ),
                    mesh,
                    dimensioned<Type>
                    (
                        vf.name(),
                        vf.dimensions(),
                        Zero
                    )
                )
            );

            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
            {
                tmp<GeometricField<scalar, fvPatchField, volMesh>> vfcomponent
                (
                    vf.component(cmpt)
                );

                word vfcName(vfcomponent().name());

                tsfCorr.ref().replace
                (
                    cmpt,
                    scv() & linear
                    <
                        typename outerProduct
                        <
                            vector,
                            typename pTraits<Type>::cmptType
                        >::type
                    > (mesh).interpolate
                    (
                        tgradScheme_().grad
                        (
                            vfcomponent,
                            "skewCorrectionGrad("+ vfcName +")"
                        )
                    )
                );
            }

            return tsfCorr;
        }


        //- Return the explicit correction to the face-interpolate
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            if
            (
                tScheme_().corrected()
             && skewCorrectionVectors::New(this->mesh()).skew()
            )
            {
                return tScheme_().correction(vf) + skewCorrection(vf);
            }
            else if (tScheme_().corrected())
            {
                return tScheme_().correction(vf);
            }
            else if (skewCorrectionVectors::New(this->mesh()).skew())
            {
                return skewCorrection(vf);
            }
            else
            {
                return
                    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
                    (
                        nullptr
                    );
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
