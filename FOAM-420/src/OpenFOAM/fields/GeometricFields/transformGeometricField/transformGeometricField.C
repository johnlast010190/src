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
    (c) 2018 OpenCFD Ltd.
Description
    Spatial transformation functions for FieldFields.

\*---------------------------------------------------------------------------*/

#include "fields/GeometricFields/transformGeometricField/transformGeometricField.H"
#include "fields/Fields/transformField/transformField.H"
#include "fields/FieldFields/transformFieldField/transformFieldField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void transform
(
    GeometricField<Type, PatchField, GeoMesh>& rtf,
    const GeometricField<tensor, PatchField, GeoMesh>& trf,
    const GeometricField<Type, PatchField, GeoMesh>& tf
)
{
    transform
    (
        rtf.primitiveFieldRef(),
        trf.primitiveField(),
        tf.primitiveField()
    );
    transform
    (
        rtf.boundaryFieldRef(),
        trf.boundaryField(),
        tf.boundaryField()
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> transform
(
    const GeometricField<tensor, PatchField, GeoMesh>& trf,
    const GeometricField<Type, PatchField, GeoMesh>& tf
)
{
    tmp<GeometricField<Type, PatchField, GeoMesh>> tranf
    (
        new GeometricField<Type, PatchField, GeoMesh>
        (
            IOobject
            (
                "transform(" + trf.name() + ',' + tf.name() + ')',
                tf.instance(),
                tf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tf.mesh(),
            tf.dimensions()
        )
    );

    transform(tranf.ref(), trf, tf);

    return tranf;
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> transform
(
    const GeometricField<tensor, PatchField, GeoMesh>& trf,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& ttf
)
{
    tmp<GeometricField<Type, PatchField, GeoMesh>> tranf =
        transform(trf, ttf());
    ttf.clear();
    return tranf;
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> transform
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh>>& ttrf,
    const GeometricField<Type, PatchField, GeoMesh>& tf
)
{
    tmp<GeometricField<Type, PatchField, GeoMesh>> tranf =
        transform(ttrf(), tf);
    ttrf.clear();
    return tranf;
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> transform
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh>>& ttrf,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& ttf
)
{
    tmp<GeometricField<Type, PatchField, GeoMesh>> tranf =
        transform(ttrf(), ttf());
    ttf.clear();
    ttrf.clear();
    return tranf;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void transform
(
    GeometricField<Type, PatchField, GeoMesh>& rtf,
    const dimensionedTensor& t,
    const GeometricField<Type, PatchField, GeoMesh>& tf
)
{
    transform(rtf.primitiveFieldRef(), t.value(), tf.primitiveField());
    transform(rtf.boundaryFieldRef(), t.value(), tf.boundaryField());
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> transform
(
    const dimensionedTensor& t,
    const GeometricField<Type, PatchField, GeoMesh>& tf
)
{
    tmp<GeometricField<Type, PatchField, GeoMesh>> tranf
    (
        new GeometricField<Type, PatchField, GeoMesh>
        (
            IOobject
            (
                "transform(" + t.name() + ',' + tf.name() + ')',
                tf.instance(),
                tf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tf.mesh(),
            tf.dimensions()
        )
    );

    transform(tranf.ref(), t, tf);

    return tranf;
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> transform
(
    const dimensionedTensor& t,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& ttf
)
{
    tmp<GeometricField<Type, PatchField, GeoMesh>> tranf =
        transform(t, ttf());
    ttf.clear();
    return tranf;
}

// * * * * * * * * * * * invTransform Global Functions * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void invTransform
(
    GeometricField<Type, PatchField, GeoMesh>& result,
    const dimensionedTensor& rot,
    const GeometricField<Type, PatchField, GeoMesh>& fld
)
{
    invTransform
    (
        result.primitiveFieldRef(),
        rot.value(),
        fld.primitiveField()
    );
    invTransform
    (
        result.boundaryFieldRef(),
        rot.value(),
        fld.boundaryField()
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
void invTransform
(
    GeometricField<Type, PatchField, GeoMesh>& result,
    const GeometricField<tensor, PatchField, GeoMesh>& rot,
    const GeometricField<Type, PatchField, GeoMesh>& fld
)
{
    invTransform
    (
        result.primitiveFieldRef(),
        rot.primitiveField(),
        fld.primitiveField()
    );
    invTransform
    (
        result.boundaryFieldRef(),
        rot.boundaryField(),
        fld.boundaryField()
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>>
invTransform
(
    const GeometricField<tensor, PatchField, GeoMesh>& rot,
    const GeometricField<Type, PatchField, GeoMesh>& fld
)
{
    auto tresult = tmp<GeometricField<Type, PatchField, GeoMesh>>::New
    (
        IOobject
        (
            "invTransform(" + rot.name() + ',' + fld.name() + ')',
            fld.instance(),
            fld.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fld.mesh(),
        fld.dimensions()
    );

    invTransform(tresult.ref(), rot, fld);

    return tresult;
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>>
invTransform
(
    const GeometricField<tensor, PatchField, GeoMesh>& rot,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tfld
)
{
    auto tresult = invTransform(rot, tfld());
    tfld.clear();
    return tresult;
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>>
invTransform
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh>>& trot,
    const GeometricField<Type, PatchField, GeoMesh>& fld
)
{
    auto tresult = invTransform(trot(), fld);
    trot.clear();
    return tresult;
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>>
invTransform
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh>>& trot,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tfld
)
{
    auto tresult = invTransform(trot(), tfld());
    tfld.clear();
    trot.clear();
    return tresult;
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>>
invTransform
(
    const dimensionedTensor& rot,
    const GeometricField<Type, PatchField, GeoMesh>& fld
)
{
    auto tresult = tmp<GeometricField<Type, PatchField, GeoMesh>>::New
    (
        IOobject
        (
            "invTransform(" + rot.name() + ',' + fld.name() + ')',
            fld.instance(),
            fld.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fld.mesh(),
        fld.dimensions()
    );

    invTransform(tresult.ref(), rot, fld);

    return tresult;
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>>
invTransform
(
    const dimensionedTensor& rot,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tfld
)
{
    auto tresult = invTransform(rot, tfld());
    tfld.clear();
    return tresult;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
