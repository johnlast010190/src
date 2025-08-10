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

#include "cfdTools/general/sensor/set/set.H"
#include "sets/topoSetSource/topoSetSource.H"
#include "sets/topoSets/cellSet.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::sensorTypes::set<Type>::set
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    sensor<Type>(mesh, dict),
    setType_(dict.lookup("setType")),
    sensorOp_(dict.lookup("sensorOp"))
{}


template<class Type>
Foam::sensorTypes::set<Type>::set(const set<Type>& cv)
:
    sensor<Type>(cv),
    setType_(cv.setType_),
    sensorOp_(cv.sensorOp_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::sensorTypes::set<Type>::~set()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::sensorTypes::set<Type>::value
(
    const GeometricField<Type, fvPatchField, volMesh>& volField
) const
{
    // initialize cell set
    labelList cellList(0);

    // initialize sensor value
    Type sensorValue = pTraits<Type>::zero;

    // fill cell set from set or zone
    if (setType_ == "cellZone")
    {
        // read cell zone from mesh given the zone name
        word cellZoneName = this->dict_.lookup("cellZone");
        label cellZoneID = volField.mesh().cellZones().findZoneID(cellZoneName);
        cellList = volField.mesh().cellZones()[cellZoneID];
    }
    else if (setType_ == "cellSet")
    {
        // read cell set from constant/polyMesh/sets/<cellSetName>
        cellSet setCellSet
        (
            volField.mesh(),
            this->dict_.lookup("cellSet"),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );


        label idx = 0;
        cellList.resize(setCellSet.size());
        forAllConstIter(labelHashSet,  setCellSet, iter)
        {
            cellList[idx] = iter.key();
            idx++;
        }
    }
    else
    {
        FatalErrorInFunction
            << "setType " << setType_ << " not supported! " << nl
            << " Valid set types are:" << nl
            << "cellSet" << nl
            << "cellZone" << nl
            << exit(FatalError);
    }

    // compute sensor values
    if (sensorOp_ == "mean") // volumetric mean
    {
        scalarField volume(volField.mesh().V());
        scalar sumVolume = 0.0;
        forAll(cellList, idx)
        {
            sensorValue += volField[cellList[idx]] * volume[cellList[idx]];
            sumVolume += volume[cellList[idx]];
        }

        reduce(sumVolume, sumOp<scalar>());
        reduce(sensorValue, sumOp<Type>());

        sensorValue /= sumVolume;
    }
    else if (sensorOp_ == "min")
    {
        forAll(cellList, idx)
        {
            if (volField[cellList[idx]] < sensorValue)
            {
                sensorValue = volField[cellList[idx]];
            }
        }
        reduce(sensorValue, minOp<Type>());
    }
    else if (sensorOp_ == "max")
    {
        forAll(cellList, idx)
        {
            if (volField[cellList[idx]] > sensorValue)
            {
                sensorValue = volField[cellList[idx]];
            }
        }
        reduce(sensorValue, maxOp<Type>());
    }
    else
    {
        FatalErrorInFunction
            << "Invalid type " << sensorOp_ << " ! " << nl
            << " Valid sensor types are:" << nl
            << "min" << nl
            << "max" << nl
            << "mean" << nl
            << exit(FatalError);
    }

    return sensorValue;
}


// ************************************************************************* //
