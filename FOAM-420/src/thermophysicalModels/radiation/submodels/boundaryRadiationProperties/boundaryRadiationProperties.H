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
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2010-2016 Esi Ltd.

Class
    Foam::radiation::boundaryRadiationProperties

Description
    Boundary radiation properties holder

SourceFiles
    boundaryRadiationProperties.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryRadiationProperties_H
#define boundaryRadiationProperties_H

#include "meshes/MeshObject/MeshObject.H"
#include "submodels/boundaryRadiationProperties/boundaryRadiationPropertiesPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

namespace radiation
{

/*---------------------------------------------------------------------------*\
                 Class boundaryRadiationProperties Declaration
\*---------------------------------------------------------------------------*/

class boundaryRadiationProperties
:
    public MeshObject
    <
        fvMesh,
        Foam::GeometricMeshObject,
        boundaryRadiationProperties
    >
{
    // Private data

        //- Ptr list of boundaryRadiationProperties
        List<autoPtr<boundaryRadiationPropertiesPatch>>
            radBoundaryPropertiesPtrList_;


public:

    // Declare name of the class and its debug switch
    TypeName("boundaryRadiationProperties");


    // Constructors

        //- Construct given fvMesh
        explicit boundaryRadiationProperties(const fvMesh&);


    // Member Functions

        //- Access boundary emissivity on patch
        tmp<scalarField> emissivity
        (
            const label patchI,
            const label bandI = 0
        ) const;


        //- Access boundary absorptivity on patch
        tmp<scalarField> absorptivity
        (
            const label patchI,
            const label bandI = 0
        ) const;


        //- Access boundary transmissivity on patch
        tmp<scalarField> transmissivity
        (
            const label patchI,
            const label bandI = 0
        ) const;

        //- Access boundary reflectivity on patch
        tmp<scalarField> reflectivity
        (
            const label patchI,
            const label bandI = 0
        ) const;

        //- Return emitted radiant power for given temperature
        tmp<scalarField> emittedRadiantFlux
        (
            const scalarField& T,
            const label patchI,
            const label bandI = 0
        ) const;

        //- Return non-interacting radiant source flux
        tmp<scalarField> radiantTransmissionSource
        (
            const label patchI,
            const label bandI = 0
        ) const;

        //- Access to radBoundaryProperties
        const List<autoPtr<boundaryRadiationPropertiesPatch>>&
        radBoundaryProperties() const;

    //- Destructor
    ~boundaryRadiationProperties();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
