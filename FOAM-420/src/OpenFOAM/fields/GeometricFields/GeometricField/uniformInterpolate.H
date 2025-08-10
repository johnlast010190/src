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
    (c) 2012-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fields/GeometricFields/GeometricField/GeometricField.H"
#include "containers/HashTables/HashPtrTable/HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * //

//- Interpolate selected fields (given by indices and corresponding weights)
//  (boundary type becomes calculated). Fields stored per index. Field gets name
//  "uniformInterpolate(" + fld.name() + ')'.
template<class GeoField>
tmp<GeoField> uniformInterpolate
(
    const HashPtrTable<GeoField, label, Hash<label>>& fields,
    const labelList& indices,
    const scalarField& weights
);

//- Interpolate fields. fieldsCache contains per timeName all loaded fields.
//  Resulting field gets properties according to fieldIO
template<class GeoField>
tmp<GeoField> uniformInterpolate
(
    const IOobject& fieldIO,
    const word& fieldName,
    const wordList& times,
    const scalarField& weights,
    const objectRegistry& fieldsCache
);

//- Interpolate fields. fieldsCache contains per timeName all loaded fields.
//  Resulting field gets properties according to fieldIO
template<class GeoField>
tmp<GeoField> uniformInterpolate
(
    const IOobject& fieldIO,
    const word& fieldName,
    const wordList& times,
    const scalarField& weights,
    const word& registryName = "fieldsCache"
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fields/GeometricFields/GeometricField/uniformInterpolate.C"
#endif

// ************************************************************************* //
