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

\*---------------------------------------------------------------------------*/

#include "sampledSurface/sampledPatch/sampledPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledPatch::sampleField
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
) const
{
    // One value per face
    tmp<Field<Type>> tvalues(new Field<Type>(patchFaceLabels_.size()));
    Field<Type>& values = tvalues.ref();
    forAll(patchFaceLabels_, i)
    {
        label patchi = patchIDs_[patchIndex_[i]];
        const Field<Type>& bField = vField.boundaryField()[patchi];
        values[i] = bField[patchFaceLabels_[i]];
    }

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledPatch::sampleField
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sField
) const
{
    // One value per face
    tmp<Field<Type>> tvalues(new Field<Type>(patchFaceLabels_.size()));
    Field<Type>& values = tvalues.ref();

    forAll(patchFaceLabels_, i)
    {
        label patchi = patchIDs_[patchIndex_[i]];
        values[i] = sField.boundaryField()[patchi][patchFaceLabels_[i]];
    }

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledPatch::interpolateField
(
    const interpolation<Type>& interpolator
) const
{
    // One value per vertex
    tmp<Field<Type>> tvalues(new Field<Type>(points().size()));
    Field<Type>& values = tvalues.ref();

    const labelList& own = mesh().faceOwner();

    boolList pointDone(points().size(), false);

    forAll(faces(), cutFacei)
    {
        label patchi = patchIDs_[patchIndex_[cutFacei]];
        const polyPatch& pp = mesh().boundaryMesh()[patchi];
        label patchFacei = patchFaceLabels()[cutFacei];
        const face& f = faces()[cutFacei];

        forAll(f, faceVertI)
        {
            label pointi = f[faceVertI];

            if (!pointDone[pointi])
            {
                label facei = patchFacei + pp.start();
                label celli = own[facei];

                values[pointi] = interpolator.interpolate
                (
                    points()[pointi],
                    celli,
                    facei
                );
                pointDone[pointi] = true;
            }
        }
    }

    return tvalues;
}


// ************************************************************************* //
