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

#include "cfdTools/general/sensor/cellValue/cellValue.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::sensorTypes::cellValue<Type>::cellValue
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    sensor<Type>(mesh, dict),
    cellID_()
{

    Istream& is(dict.lookup("cellID"));
    is  >> cellID_;
}


template<class Type>
Foam::sensorTypes::cellValue<Type>::cellValue(const cellValue<Type>& cv)
:
    sensor<Type>(cv),
    cellID_(cv.cellID_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::sensorTypes::cellValue<Type>::~cellValue()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::sensorTypes::cellValue<Type>::value
(
    const GeometricField<Type, fvPatchField, volMesh>& volField
) const
{
    // cellID is global index
    // find local idx matching global ID in processor addressing
    label globalCellID(cellID_);
    label localCellID(-1);

    // initialize field value
    Type fieldValue = pTraits<Type>::zero;

    if (Pstream::parRun())
    {
        labelIOList cellProcAddressing
        (
            IOobject
            (
                "cellProcAddressing",
                volField.mesh().facesInstance(),
                volField.mesh().meshSubDir,
                volField.mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        forAll(cellProcAddressing, cI)
        {
            if (cellProcAddressing[cI] == globalCellID)
            {
                localCellID = cI;
            }
        }

        if (localCellID > -1)
        {
            fieldValue = volField[localCellID];
        }

        reduce(fieldValue, sumOp<Type>());
    }
    else
    {
        fieldValue = volField[globalCellID];
    }

    return fieldValue;
}


// ************************************************************************* //
