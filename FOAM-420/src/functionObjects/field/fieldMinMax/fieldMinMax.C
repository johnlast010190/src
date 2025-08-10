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
    (c) 2015-2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "fieldMinMax/fieldMinMax.H"
#include "fields/Fields/fieldTypes.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldMinMax, 0);
    addToRunTimeSelectionTable(functionObject, fieldMinMax, dictionary);
}
}

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::fieldMinMax::modeType,
    2
>::names[] = {"magnitude", "component"};

const Foam::NamedEnum
<
    Foam::functionObjects::fieldMinMax::modeType,
    2
> Foam::functionObjects::fieldMinMax::modeTypeNames_;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fieldMinMax::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Field minima and maxima");
    writeCommented(os, "Time");

    if (location_)
    {
        writeDelimited(os, "field");
        writeDelimited(os, "min");
        writeDelimited(os, "location(min)");

        if (Pstream::parRun())
        {
            writeDelimited(os, "processor");
        }

        writeDelimited(os, "max");
        writeDelimited(os, "location(max)");

        if (Pstream::parRun())
        {
            writeDelimited(os, "processor");
        }
    }
    else
    {
        forAll(fieldSet_, fieldi)
        {
            writeDelimited(os, "min(" + fieldSet_[fieldi] + ')');
            writeDelimited(os, "max(" + fieldSet_[fieldi] + ')');
        }
    }

    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldMinMax::fieldMinMax
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    location_(true),
    mode_(mdMag),
    fieldSet_()
{
    read(dict);
    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldMinMax::~fieldMinMax()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldMinMax::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    location_ = dict.lookupOrDefault<Switch>("location", true);

    mode_ = modeTypeNames_[dict.lookupOrDefault<word>("mode", "magnitude")];
    dict.lookup("fields") >> fieldSet_;

    return true;
}


bool Foam::functionObjects::fieldMinMax::execute()
{
    if (!location_) writeTime(file());
    Log << type() << " " << name() <<  " execute:" << nl;

    forAll(fieldSet_, fieldi)
    {
        calcMinMaxFields<scalar>(fieldSet_[fieldi], mdCmpt);
        calcMinMaxFields<vector>(fieldSet_[fieldi], mode_);
        calcMinMaxFields<sphericalTensor>(fieldSet_[fieldi], mode_);
        calcMinMaxFields<symmTensor>(fieldSet_[fieldi], mode_);
        calcMinMaxFields<tensor>(fieldSet_[fieldi], mode_);
    }

    if (!location_) file()<< endl;
    Log << endl;

    return true;
}


bool Foam::functionObjects::fieldMinMax::write()
{
    return true;
}


// ************************************************************************* //
