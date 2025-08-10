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
    (c) 2015 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/sensor/maxFaceWeight/maxFaceWeight.H"
#include "cfdTools/general/include/fvCFD.H"
#include "global/constants/mathematical/mathematicalConstants.H"
#include "meshes/polyMesh/polyMeshCheck/polyMeshTools.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::sensorTypes::maxFaceWeight<Type>::maxFaceWeight
(
    const fvMesh& mesh,
    const dictionary& sensorDict
)
:
    sensor<Type>(mesh, sensorDict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::sensorTypes::maxFaceWeight<Type>::~maxFaceWeight()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Foam::scalar>>
Foam::sensorTypes::maxFaceWeight<Type>::valueField() const
{
    const fvMesh& mesh = this->mesh_;

    tmp<scalarField> tfaceWeight =
        polyMeshTools::faceWeights(mesh, mesh.Cf(), mesh.Sf(), mesh.C());
    const scalarField& faceWeight = tfaceWeight();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    scalarField volMaxFaceWeight(mesh.nCells(), - SMALL);

    forAll(owner, fI)
    {
        label oI = owner[fI];
        label nI = neighbour[fI];
        scalar fO = faceWeight[fI];

        if (volMaxFaceWeight[oI] < fO)
        {
            volMaxFaceWeight[oI] = fO;
        }
        if (volMaxFaceWeight[nI] < fO)
        {
            volMaxFaceWeight[nI] = fO;
        }
    }

    return tmp<scalarField>
    (
        new scalarField(2.0*volMaxFaceWeight)
    );
}


// ************************************************************************* //
