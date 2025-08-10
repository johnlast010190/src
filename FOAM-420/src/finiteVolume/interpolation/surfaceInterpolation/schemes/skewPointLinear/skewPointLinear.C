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

\*---------------------------------------------------------------------------*/

#include "interpolation/surfaceInterpolation/schemes/skewPointLinear/skewPointLinear.H"
#include "fvMesh/fvMesh.H"
#include "interpolation/volPointInterpolation/volPointInterpolation.H"
#include "meshes/primitiveShapes/triangle/triangle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::skewPointLinear<Type>::
correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    GeometricField<Type, pointPatchField, pointMesh> pvf
    (
        volPointInterpolation::New(mesh).interpolate(vf)
    );

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsfCorr =
        linearInterpolate(vf);

    Field<Type>& sfCorr = tsfCorr.ref().primitiveFieldRef();

    typename GeometricField<Type, fvsPatchField, surfaceMesh>
        ::Boundary& sfCorrbf = tsfCorr.ref().boundaryFieldRef();

    //We dont want to use cell based interpolation sources to get the face
    //value at all, only volPoint interpolation
    //Thus we start with the face centroid and implicitly move the
    //value at the centroid to the left hand side of the integral
    //We assume the value at the centroid is the same as the
    //face average value

    //The scheme thus has integral skew correction and removes the
    //dependence on adjacent cell values that would have maintained
    //the inherent issues of the linear interpolation schemes


    const pointField& points = mesh.points();
    const surfaceVectorField& Cf = mesh.Cf();
    const faceList& faces = mesh.faces();

    forAll(sfCorr, facei)
    {

        scalar at = triangle<point, const point&>
        (
            Cf[facei],
            points[faces[facei][0]],
            points[faces[facei][faces[facei].size()-1]]
        ).mag();

        scalar sumAt = at;
        Type sumPsip = at*(1.0/3.0)*
        (
            pvf[faces[facei][0]]
          + pvf[faces[facei][faces[facei].size()-1]]
        );

        for (label pointi=1; pointi<faces[facei].size(); pointi++)
        {
            at = triangle<point, const point&>
            (
                Cf[facei],
                points[faces[facei][pointi]],
                points[faces[facei][pointi-1]]
            ).mag();

            sumAt += at;
            sumPsip += at*(1.0/3.0)*
            (
                pvf[faces[facei][pointi]]
              + pvf[faces[facei][pointi-1]]
            );

        }

        sfCorr[facei] = 1.5*sumPsip/sumAt - sfCorr[facei];
    }

    forAll(tsfCorr().boundaryField(),pI)
    {
        if (tsfCorr().boundaryField()[pI].coupled())
        {
            forAll(tsfCorr().boundaryField()[pI],pfI)
            {
                label psf = tsfCorr().boundaryField()[pI].patch().start();
                label facei = psf + pfI;

                scalar at = triangle<point, const point&>
                (
                    Cf.boundaryField()[pI][pfI],
                    points[faces[facei][0]],
                    points[faces[facei][faces[facei].size()-1]]
                ).mag();

                scalar sumAt = at;
                Type sumPsip = at*(1.0/3.0)*
                (
                    pvf[faces[facei][0]]
                  + pvf[faces[facei][faces[facei].size()-1]]
                );

                for (label pointi=1; pointi<faces[facei].size(); pointi++)
                {
                    at = triangle<point, const point&>
                    (
                        Cf.boundaryField()[pI][pfI],
                        points[faces[facei][pointi]],
                        points[faces[facei][pointi-1]]
                    ).mag();

                    sumAt += at;
                    sumPsip += at*(1.0/3.0)*
                    (
                        pvf[faces[facei][pointi]]
                      + pvf[faces[facei][pointi-1]]
                    );

                }
                sfCorrbf[pI][pfI] = 1.5*sumPsip/sumAt -
                    tsfCorr().boundaryField()[pI][pfI];
            }
        }
        else
        {
            sfCorrbf[pI] = pTraits<Type>::zero;
        }
    }

    return tsfCorr;
}


namespace Foam
{
    makeSurfaceInterpolationScheme(skewPointLinear);
}

// ************************************************************************* //
