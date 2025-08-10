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

\*---------------------------------------------------------------------------*/

#include "submodels/boundaryRadiationProperties/boundaryRadiationProperties.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(boundaryRadiationProperties, 0);
    }
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::radiation::boundaryRadiationProperties::boundaryRadiationProperties
(
    const fvMesh& mesh
)
:
    MeshObject
    <
        fvMesh,
        Foam::GeometricMeshObject,
        boundaryRadiationProperties
    >(mesh),
    radBoundaryPropertiesPtrList_(mesh.boundary().size())
{
    IOobject boundaryIO
    (
        boundaryRadiationProperties::typeName,
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (boundaryIO.typeHeaderOk<IOdictionary>(true))
    {
        const IOdictionary radiationDict(boundaryIO);

        forAll(mesh.boundary(), patchi)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchi];

            if (radiationDict.isDict(pp.name()))
            {
                const dictionary& dict = radiationDict.subDict(pp.name());

                radBoundaryPropertiesPtrList_[patchi].set
                (
                    new boundaryRadiationPropertiesPatch(pp, dict)
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member fucntions * * * * * * * * * * * * *  //


const Foam::List<Foam::autoPtr<Foam::radiation::boundaryRadiationPropertiesPatch>>&
Foam::radiation::boundaryRadiationProperties::radBoundaryProperties() const
{
    return radBoundaryPropertiesPtrList_;
}


Foam::tmp<Foam::scalarField>
Foam::radiation::boundaryRadiationProperties::emissivity
(
    const label patchi,
    const label bandi
) const
{
    if (!radBoundaryPropertiesPtrList_[patchi].empty())
    {
        return radBoundaryPropertiesPtrList_[patchi]->emissivity(bandi);
    }
    else
    {
        FatalErrorInFunction
            << "Patch : " << mesh().boundaryMesh()[patchi].name()
            << " is not found in the boundaryRadiationProperties. "
            << "Please add it"
            << exit(FatalError);

        return tmp<scalarField>(new scalarField());
    }
}


Foam::tmp<Foam::scalarField>
Foam::radiation::boundaryRadiationProperties::absorptivity
(
    const label patchi,
    const label bandi
) const
{
    if (!radBoundaryPropertiesPtrList_[patchi].empty())
    {
        return radBoundaryPropertiesPtrList_[patchi]->absorptivity(bandi);
    }
    else
    {
        FatalErrorInFunction
            << "Patch : " << mesh().boundaryMesh()[patchi].name()
            << " is not found in the boundaryRadiationProperties. "
            << "Please add it"
            << exit(FatalError);

        return tmp<scalarField>(new scalarField());
    }
}


Foam::tmp<Foam::scalarField>
Foam::radiation::boundaryRadiationProperties::transmissivity
(
    const label patchi,
    const label bandi
) const
{
    if (!radBoundaryPropertiesPtrList_[patchi].empty())
    {
        return radBoundaryPropertiesPtrList_[patchi]->transmissivity(bandi);
    }
    else
    {
        FatalErrorInFunction
            << "Patch : " << mesh().boundaryMesh()[patchi].name()
            << " is not found in the boundaryRadiationProperties. "
            << "Please add it"
            << exit(FatalError);

        return tmp<scalarField>(new scalarField());
    }
}


Foam::tmp<Foam::scalarField>
Foam::radiation::boundaryRadiationProperties::reflectivity
(
    const label patchi,
    const label bandi
) const
{
    if (!radBoundaryPropertiesPtrList_[patchi].empty())
    {
        return radBoundaryPropertiesPtrList_[patchi]->reflectivity(bandi);
    }
    else
    {
        FatalErrorInFunction
            << "Patch : " << mesh().boundaryMesh()[patchi].name()
            << " is not found in the boundaryRadiationProperties. "
            << "Please add it"
            << exit(FatalError);

        return tmp<scalarField>(new scalarField());
    }
}

Foam::tmp<Foam::scalarField> Foam::radiation::boundaryRadiationProperties::
emittedRadiantFlux
(
    const scalarField& T,
    const label patchI,
    const label bandI
) const
{
    if (!radBoundaryPropertiesPtrList_[patchI].empty())
    {
        return radBoundaryPropertiesPtrList_[patchI]->emittedRadiantFlux(T, bandI);
    }
    else
    {
        FatalErrorInFunction
            << "Patch : " << mesh().boundaryMesh()[patchI].name()
            << " is not found in the boundaryRadiationProperties. "
            << "Please add it"
            << exit(FatalError);

        return tmp<scalarField>(new scalarField());
    }
}

Foam::tmp<Foam::scalarField> Foam::radiation::boundaryRadiationProperties::
radiantTransmissionSource
(
    const label patchI,
    const label bandI
) const
{
    if (!radBoundaryPropertiesPtrList_[patchI].empty())
    {
        return radBoundaryPropertiesPtrList_[patchI]->radiantTransmissionSource(bandI);
    }
    else
    {
        FatalErrorInFunction
            << "Patch : " << mesh().boundaryMesh()[patchI].name()
            << " is not found in the boundaryRadiationProperties. "
            << "Please add it"
            << exit(FatalError);

        return tmp<scalarField>(new scalarField());
    }
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::radiation::boundaryRadiationProperties::~boundaryRadiationProperties()
{}


// ************************************************************************* //
