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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "interRegionOption.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(interRegionOption, 0);
    }
}


// * * * * * * * * * * * *  Protected member functions * * * * * * * * * * * //

void Foam::fv::interRegionOption::setMapper()
{
    if (master_)
    {
        Info<< indent << "- selecting inter region mapping" << endl;

        const fvMesh& nbrMesh =
            mesh_.time().lookupObject<fvMesh>(nbrRegionName_);

        if (mesh_.name() == nbrMesh.name())
        {
            FatalErrorInFunction
                << "Inter-region model selected, but local and "
                << "neighbour regions are the same: " << nl
                << "    local region: " << mesh_.name() << nl
                << "    secondary region: " << nbrMesh.name() << nl
                << exit(FatalError);
        }

        if (mesh_.bounds().overlaps(nbrMesh.bounds()))
        {
            meshInterpPtr_.reset
            (
                new meshToMesh
                (
                    mesh_,
                    nbrMesh,
                    meshToMesh::interpolationMethodNames_.read
                    (
                        coeffs_.lookup("interpolationMethod")
                    ),
                    false // not interpolating patches
                )
            );
        }
        else
        {
            FatalErrorInFunction
                << "regions " << mesh_.name() << " and "
                << nbrMesh.name() <<  " do not intersect"
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::interRegionOption::interRegionOption
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    option
    (
        name,
        modelType,
        dict,
        obr
    ),
    master_(coeffs_.lookupOrDefault<bool>("master", true)),
    nbrRegionName_(coeffs_.lookup("nbrRegion")),
    meshInterpPtr_()
{
    if (obr_.time().found("regionProperties"))
    {
        const IOdictionary& regionProperties =
            obr_.time().lookupObject<IOdictionary>("regionProperties");
        const wordList inactiveRegions
        (
            regionProperties.lookupOrDefault
            (
                "inactiveRegions",
                wordList::null()
            )
        );
        if
        (
            inactiveRegions.found(mesh_.name())
         || inactiveRegions.found(nbrRegionName_)
        )
        {
            active() = false;
        }
    }
    if (active())
    {
        setMapper();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::interRegionOption::~interRegionOption()
{}


// ************************************************************************* //
