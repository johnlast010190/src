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

#include "fieldValues/fieldValueDelta/fieldValueDelta.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{
    defineTypeNameAndDebug(fieldValueDelta, 0);
    addToRunTimeSelectionTable(functionObject, fieldValueDelta, dictionary);
}
}
}

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::fieldValueDelta::operationType,
    5
>::names[] =
{
    "add",
    "subtract",
    "min",
    "max",
    "average"
};
const Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::fieldValueDelta::operationType,
    5
> Foam::functionObjects::fieldValues::fieldValueDelta::operationTypeNames_;


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::fieldValues::fieldValueDelta::writeFileHeader
(
    Ostream& os
) const
{
    writeHeaderValue(os, "Source1", operand1_);
    writeHeaderValue(os, "Source2", operand2_);
    writeHeaderValue(os, "Operation", operationTypeNames_[operation_]);
    writeCommented(os, "Time");


    if (columnSpecified_)
    {
        const wordList entries1 = {inputColumn1_};
        const wordList entries2 = {inputColumn2_};

        forAll(entries1, i)
        {
            os  << tab << entries1[i] << " " << operationTypeNames_[operation_]
                << " " << entries2[i];
        }
    }
    else
    {
        const wordList entries1 = objectResultEntries(operand1_);

        forAll(entries1, i)
        {
            os  << tab << entries1[i];
        }
    }

    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::fieldValueDelta::fieldValueDelta
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    operation_(opSubtract),
    region1Ptr_(nullptr),
    region2Ptr_(nullptr),
    columnSpecified_(false)
{
    read(dict);
    mustWriteHeader_ = true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::fieldValueDelta::~fieldValueDelta()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::fieldValues::fieldValueDelta::operate
(
    const wordList entries1, const wordList entries2
)
{
    if (entries1.size() == 0 || entries2.size() == 0)
    {
        FatalErrorInFunction
                << "Could not apply field value delta operation: "
                << "One or both of the function objects used with the "
                << "field value delta function object does not contain any "
                << "post-processing results."
                << exit(FatalError);
    }
    
    if ((entries1.size() != entries2.size()) && !columnSpecified_)
    {
        FatalErrorInFunction
            << name() << ": objects must generate the same number of results"
            << nl
            << "    " << operand1_ << " objects: " << entries1 << nl
            << "    " << operand2_ << " objects: " << entries2 << nl
            << exit(FatalError);
    }
    
    forAll(entries1, i)
    {
        const word& entry1(entries1[i]);
        const word& entry2(entries2[i]);
        const word type1 = objectResultType(operand1_, entry1);
        const word type2 = objectResultType(operand2_, entry2);

        if (type1 != type2)
        {
            FatalErrorInFunction
                << name()
                << ": input values for operation must be of the same type"
                << nl
                << "    " << entry1 << ": " << type1 << nl
                << "    " << entry2 << ": " << type2 << nl
                << exit(FatalError);
        }

        bool found = false;

        applyOperation<scalar>(type1, operand1_, operand2_, entry1, entry2, found);
        applyOperation<vector>(type1, operand1_, operand2_, entry1, entry2, found);
        applyOperation<sphericalTensor>
            (type1, operand1_, operand2_, entry1, entry2, found);
        applyOperation<symmTensor>(type1, operand1_, operand2_, entry1, entry2, found);
        applyOperation<tensor>(type1, operand1_, operand2_, entry1, entry2, found);

        if (!found)
        {
            Log << "Operation between "
                << operand1_ << " with result " << entry1 << " and "
                << operand2_ << " with result " << entry2 << " not applied"
                << endl;
        }
    }

    Log << (entries1.empty() ? "    none" : "") << endl;

    file()<< endl;
}

bool Foam::functionObjects::fieldValues::fieldValueDelta::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    if (dict.found("inputFO1") && dict.found("inputFO2"))
    {
        operand1_ = dict.lookup<word>("inputFO1");
        operand2_ = dict.lookup<word>("inputFO2");

        log = dict.lookupOrDefault<Switch>("log", false);
    }
    else
    {
        region1Ptr_.reset
        (
            fieldValue::New
            (
                name() + ".region1",
                obr_,
                dict.subDict("region1"),
                false
            ).ptr()
        );
        region2Ptr_.reset
        (
            fieldValue::New
            (
                name() + ".region2",
                obr_,
                dict.subDict("region2"),
                false
            ).ptr()
        );

        operand1_ = region1Ptr_->name();
        operand2_ = region2Ptr_->name();
    }

    if (dict.found("inputColumn1") && dict.found("inputColumn2"))
    {
        inputColumn1_ = dict.lookup<word>("inputColumn1");
        inputColumn2_ = dict.lookup<word>("inputColumn2");
        columnSpecified_ = true;
    }
    else if (dict.found("inputColumn1") && !dict.found("inputColumn2"))
    {
        FatalErrorInFunction
                    << "Keyword inputColumn2 is undefined."
                    << exit(FatalError);
    }
    else if (!dict.found("inputColumn1") && dict.found("inputColumn2"))
    {
        FatalErrorInFunction
                    << "Keyword inputColumn1 is undefined."
                    << exit(FatalError);
    }

    operation_ = operationTypeNames_.read(dict.lookup("operation"));

    return true;
}


bool Foam::functionObjects::fieldValues::fieldValueDelta::write()
{
    if (mustWriteHeader_)
    {
        writeFileHeader(file());
        mustWriteHeader_ = false;
    }

    if (region1Ptr_.valid())
    {
        region1Ptr_->write();
        region2Ptr_->write();
    }

    writeTime(file());

    Log << type() << " " << name() << " write:" << endl;
    
    if (columnSpecified_)
    {
        const wordList entries1 = {inputColumn1_};
        const wordList entries2 = {inputColumn2_};

        if 
        (
            objectResultEntries(operand1_).size() > 0 && objectResultType
            (
                operand1_, entries1[0]
            ) == word::null 
        )
        {   
            FatalErrorInFunction
                    << "Function object "
                    << operand1_
                    << " does not have a column called "
                    << entries1[0]
                    << " Valid options are: "
                    << objectResultEntries(operand1_)
                    << exit(FatalError);
        }
        else if 
        (
            objectResultEntries(operand2_).size() > 0 && objectResultType
            (
                operand2_, entries2[0]
            ) == word::null 
        )
        {
            FatalErrorInFunction
                    << "Function object "
                    << operand2_
                    << " does not have a column called "
                    << entries2[0]
                    << " Valid options are: "
                    << objectResultEntries(operand2_)
                    << exit(FatalError);
        }
        else if (objectResultEntries(operand1_).size() == 0)
        {
            FatalErrorInFunction
                    << "Function object "
                    << operand1_
                    << " cannot be found. Available function objects are: "
                    << stateDict().subDict("results").keys()
                    << exit(FatalError);
        }
        else if (objectResultEntries(operand2_).size() == 0)
        {
            FatalErrorInFunction
                    << "Function object "
                    << operand2_
                    << " cannot be found. Available function objects are: "
                    << stateDict().subDict("results").keys()
                    << exit(FatalError);
        }
        else
        {
            operate(entries1, entries2);
        }
    }
    else
    {
        const wordList entries1 = objectResultEntries(operand1_);
        const wordList entries2 = objectResultEntries(operand2_);

        operate(entries1, entries2);
    }

    return true;
}


bool Foam::functionObjects::fieldValues::fieldValueDelta::execute()
{
    return true;
}


// ************************************************************************* //
