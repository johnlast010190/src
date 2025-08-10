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

#include "finiteVolume/snGradSchemes/LENCSnGrad/LENCSnGrad.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::LENCSnGrad<Type>::~LENCSnGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fv::LENCSnGrad<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // I really dont know why the () is necessary but without
    // its a segmentation fault
    return tnonOrthCorr_();
}


template<class Type>
Foam::tmp<Foam::surfaceScalarField> Foam::fv::LENCSnGrad<Type>::deltaCoeffs
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    //Info<< "LENC-scheme used on: " << vf.name() << endl;
    // Implementation for implicitly limited 'Over-relaxed approach'
    const fvMesh& mesh = this->mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tdeltaCoeffs
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "deltaCoeffs",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh.nonOrthDeltaCoeffs()
        )
    );

    GeometricField<Type, fvsPatchField, surfaceMesh>& nonOrthCorr
        = tnonOrthCorr_.ref();
    nonOrthCorr *= dimensionedScalar("unitVF", vf.dimensions(), 1.);

    // Hard-coded gradient calculation of 'vf' with 'leastSquares'-scheme
    IStringStream gradSchemeLS("leastSquares");
    tmp<GeometricField<typename outerProduct<vector, Type>::type,
        fvPatchField, volMesh>> tgradc = fv::gradScheme<Type>::New
        (
            mesh,
            vf.db(),
            gradSchemeLS
        )().grad(vf, "grad(" + vf.name() + ')');
    GeometricField<typename outerProduct<vector, Type>::type,
        fvPatchField, volMesh>& gradc
        = tgradc.ref();

    nonOrthCorr =
        mesh.nonOrthCorrectionVectors()
      & linear<typename outerProduct<vector, Type>::type>(mesh).interpolate
        (gradc);

    // Skewness Correction
    if (true)
    {
        // gamma is needed because skewCorrectionVectors should have the unit [m]
        // but it is dimless.
        // The unit of mesh.Sf() as well as mesh.magSf() is [m^2]
        dimensioned<scalar> unit("unit", dimensionSet(0, 1, 0, 0, 0, 0, 0), 1);
        IStringStream gradSchemeLS("leastSquares");
        nonOrthCorr +=
            (mesh.Sf()/mesh.magSf()) &
            (
                 unit*skewCorrectionVectors::New(mesh)()
               & linear<typename outerProduct<vector, vector>::type>(mesh).interpolate
                 (
                     fv::gradScheme<vector>::New
                     (
                         mesh,
                         vf.db(),
                         gradSchemeLS
                     )().grad(gradc, "grad(grad(" + vf.name() + "))")
                 )
            );
    }

    return tdeltaCoeffs;
}

// ************************************************************************* //
