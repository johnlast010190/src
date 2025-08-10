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
    (c) 2022 Esi Ltd.

InClass
    Foam::surfaceFields

Description
    Foam::surfaceFields

SourceFiles
    surfaceFields.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceFields_H
#define surfaceFields_H

#include "fields/GeometricFields/GeometricField/GeometricFields.H"
#include "surfaceMesh/surfaceMesh.H"
#include "fvMesh/fvMesh.H"
#include "fields/fvsPatchFields/fvsPatchField/fvsPatchFields.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "fields/fvsPatchFields/basic/calculated/calculatedFvsPatchFields.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

namespace Foam
{

//- Check if surfaceField is a flux
//  i.e. an extensive field oriented with the face requiring
//  flipping after mapping.
//  Currently, it is assumed that only surfaceScalarFields can be fluxes
template<class Type>
inline bool isFlux(const DimensionedField<Type, surfaceMesh>& df)
{
    return false;
}

//- Check if surfaceScalarField is a flux
//  i.e. an extensive field oriented with the face requiring
//  flipping after mapping.
//  Currently, it is assumed that fluxes are either volumetric or mass
template<>
inline bool isFlux<scalar>(const DimensionedField<scalar, surfaceMesh>& df)
{
    return
        df.dimensions() == dimArea*dimVelocity
     || df.dimensions() == dimDensity*dimArea*dimVelocity;
}

//- Check if surfaceField is a face velocity or momentum
template<class Type>
inline bool isFaceVelocity(const DimensionedField<Type, surfaceMesh>& df)
{
    return false;
}

//- Check if surfaceVectorField is a face velocity or momentum
template<>
inline bool isFaceVelocity<vector>
(
    const DimensionedField<vector, surfaceMesh>& df
)
{
    return
        df.dimensions() == dimVelocity
     || df.dimensions() == dimDensity*dimVelocity;
}

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
