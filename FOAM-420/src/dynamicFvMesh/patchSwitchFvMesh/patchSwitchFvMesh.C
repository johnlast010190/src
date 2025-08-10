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

\*---------------------------------------------------------------------------*/

#include "patchSwitchFvMesh/patchSwitchFvMesh.H"
#include "meshes/polyMesh/polyPatches/indirectPolyPatch/indirectPolyPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchSwitchFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, patchSwitchFvMesh, IOobject);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

patchSwitchFvMesh::patchSwitchFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    dynamicMeshCoeffs_(dynamicMeshDict().optionalSubDict(typeName + "Coeffs")),
    pSwitchPtr_()
{
    pSwitchPtr_.setSize(dynamicMeshCoeffs_.size());

    label zoneI = 0;

    forAllConstIter(dictionary, dynamicMeshCoeffs_, iter)
    {
        if (iter().isDict())
        {
            const dictionary& subDict = iter().dict();

            pSwitchPtr_.set
            (
                zoneI,
                patchSwitch::New
                (
                    *this,
                    subDict
                )
            );

            zoneI++;
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool patchSwitchFvMesh::update()
{
    bool meshChanged = false;
    forAll(pSwitchPtr_, psI)
    {
        bool changedI =  pSwitchPtr_[psI].update();
        if (changedI == true)
        {
            meshChanged = true;
        }
    }
    this->topoChanging(meshChanged);

    return meshChanged;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
