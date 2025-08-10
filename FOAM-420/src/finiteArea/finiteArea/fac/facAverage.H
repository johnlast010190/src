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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016 Esi Ltd.

InNamespace
    Foam::fac

Description
    Edge-weighted average a edgeField creating a areaField

SourceFiles
    facAverage.C

\*---------------------------------------------------------------------------*/


#ifndef facAverage_H
#define facAverage_H

#include "fields/areaFields/areaFieldsFwd.H"
#include "fields/edgeFields/edgeFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fac functions Declaration
\*---------------------------------------------------------------------------*/

namespace fac
{
    //- Area-weighted average a edgeField creating a areaField
    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh>> average
    (
        const GeometricField<Type,  faePatchField, faEdgeMesh>&
    );

    //- Area-weighted average a edgeField creating a areaField
    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh>> average
    (
        const tmp<GeometricField<Type,  faePatchField, faEdgeMesh>>&
    );


    //- Area-weighted average a edgeField creating a areaField
    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh>> average
    (
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    //- Area-weighted average a edgeField creating a areaField
    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh>> average
    (
        const tmp<GeometricField<Type, faPatchField, areaMesh>>&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "finiteArea/fac/facAverage.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
