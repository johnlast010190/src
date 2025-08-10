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
    (c) 2017 Wikki Ltd.

Class
    volSurfaceMapping

Description
    Volume to surface and surface to volume mapping

SourceFiles
    volSurfaceMapping.C

\*---------------------------------------------------------------------------*/

#ifndef volSurfaceMapping_H
#define volSurfaceMapping_H

#include "faMesh/faMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type> class fvPatchField;

/*---------------------------------------------------------------------------*\
                        Class volSurfaceMapping Declaration
\*---------------------------------------------------------------------------*/

class volSurfaceMapping
{
    // Private data

        //- Reference to mesh
        const faMesh& mesh_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        volSurfaceMapping(const volSurfaceMapping&);

        //- Disallow default bitwise assignment
        void operator=(const volSurfaceMapping&);


public:

    // Constructors

        //- Construct null
        volSurfaceMapping();

        //- Construct from mesh
        volSurfaceMapping(const faMesh& mesh)
        :
            mesh_(mesh)
        {}


    // Destructor - default


    // Member Functions

        //- Map droplet cloud sources to surface
        template<template<class> class PatchField, class Type>
        tmp<Field<Type>> mapToSurface
        (
            const FieldField<PatchField, Type>& df
        ) const;

        //- Map surface field to volume boundary field
        template<class Type>
        void mapToVolume
        (
            const GeometricField<Type, faPatchField, areaMesh>& af,
            FieldField<fvPatchField, Type>& bf
        ) const;

        template<class Type>
        void mapToVolume
        (
            const tmp<GeometricField<Type, faPatchField, areaMesh>>& taf,
            FieldField<fvPatchField, Type>& bf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "interpolation/volSurfaceMapping/volSurfaceMapping.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
