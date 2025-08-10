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

InClass
    Foam::transformGeometricField

Description
    Spatial transformation functions for FieldFields.

SourceFiles
    transformGeometricField.C

\*---------------------------------------------------------------------------*/

#ifndef transformGeometricField_H
#define transformGeometricField_H

#include "primitives/transform/transform.H"
#include "fields/GeometricFields/GeometricField/GeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void transform
(
    GeometricField<Type, PatchField, GeoMesh>&,
    const GeometricField<tensor, PatchField, GeoMesh>&,
    const GeometricField<Type, PatchField, GeoMesh>&
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> transform
(
    const GeometricField<tensor, PatchField, GeoMesh>&,
    const GeometricField<Type, PatchField, GeoMesh>&
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> transform
(
    const GeometricField<tensor, PatchField, GeoMesh>&,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>&
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> transform
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh>>&,
    const GeometricField<Type, PatchField, GeoMesh>&
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> transform
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh>>&,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>&
);


template<class Type, template<class> class PatchField, class GeoMesh>
void transform
(
    GeometricField<Type, PatchField, GeoMesh>&,
    const dimensionedTensor&,
    const GeometricField<Type, PatchField, GeoMesh>&
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> transform
(
    const dimensionedTensor&,
    const GeometricField<Type, PatchField, GeoMesh>&
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> transform
(
    const dimensionedTensor&,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>&
);


template<class Type, template<class> class PatchField, class GeoMesh>
void invTransform
(
    GeometricField<Type, PatchField, GeoMesh>&,
    const GeometricField<tensor, PatchField, GeoMesh>&,
    const GeometricField<Type, PatchField, GeoMesh>&
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> invTransform
(
    const GeometricField<tensor, PatchField, GeoMesh>&,
    const GeometricField<Type, PatchField, GeoMesh>&
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> invTransform
(
    const GeometricField<tensor, PatchField, GeoMesh>&,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>&
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> invTransform
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh>>&,
    const GeometricField<Type, PatchField, GeoMesh>&
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> invTransform
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh>>&,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>&
);


template<class Type, template<class> class PatchField, class GeoMesh>
void invTransform
(
    GeometricField<Type, PatchField, GeoMesh>&,
    const dimensionedTensor&,
    const GeometricField<Type, PatchField, GeoMesh>&
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> invTransform
(
    const dimensionedTensor&,
    const GeometricField<Type, PatchField, GeoMesh>&
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> invTransform
(
    const dimensionedTensor&,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>&
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fields/GeometricFields/transformGeometricField/transformGeometricField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
