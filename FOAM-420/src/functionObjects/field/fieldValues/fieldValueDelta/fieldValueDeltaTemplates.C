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
    (c) 2012-2016 OpenFOAM Foundation
    (c) 2023      Esi Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::fieldValues::fieldValueDelta::applyOperation
(
    const word& resultType,
    const word& name1,
    const word& name2,
    const word& entryName1,
    const word& entryName2,
    bool& found
)
{
    if (pTraits<Type>::typeName != resultType)
    {
        return;
    }

    Type result = Zero;

    Type value1 = Zero;
    Type value2 = Zero;

    this->getObjectResult<Type>(name1, entryName1, value1);
    this->getObjectResult<Type>(name2, entryName2, value2);

    const word& opName = operationTypeNames_[operation_];

    switch (operation_)
    {
        case opAdd:
        {
            result = value1 + value2;
            break;
        }
        case opSubtract:
        {
            result = value1 - value2;
            break;
        }
        case opMin:
        {
            result = min(value1, value2);
            break;
        }
        case opMax:
        {
            result = max(value1, value2);
            break;
        }
        case opAverage:
        {
            result = 0.5*(value1 + value2);
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unable to process operation "
                << operationTypeNames_[operation_]
                << abort(FatalError);
        }
    }

    const word resultName
    (
        opName
        + '(' + name1 + '-' + entryName1 + ',' + name2 + '-' + entryName2 + ')'
    );
    Log << "    " << resultName << " = " << result << endl;

    this->file()<< tab << result;

    // Write state/results information
    this->setResult(resultName, result);

    found = true;
}


// ************************************************************************* //
