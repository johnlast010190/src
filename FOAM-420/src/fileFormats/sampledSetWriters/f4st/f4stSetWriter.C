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
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "sampledSetWriters/f4st/f4stSetWriter.H"
#include "global/clock/clock.H"
#include "coordSet/coordSet.H"
#include "primitives/strings/fileName/fileName.H"
#include "db/IOstreams/Fstreams/OFstream.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
template<class Type>
void Foam::f4stSetWriter<Type>::createFile() const
{
    if (!isDir("f4st"))
    {
        mkDir("f4st");
    }

    word name_ = "surfaceData";
    char dataFile[80];
    sprintf(dataFile, "f4st/%s.f4", name_.c_str());

    if (outFilePtr_.valid())
    {
        outFilePtr_.clear();
    }

    if (!isFile(dataFile))
    {
        outFilePtr_ =
            new OFstream
            (
                dataFile,
                binary_ ? IOstream::BINARY : IOstream::ASCII,
                IOstream::currentVersion,
                compressed_ ? IOstream::COMPRESSED : IOstream::UNCOMPRESSED,
                false
            );
    }
    else
    {
        outFilePtr_ =
            new OFstream
            (
                dataFile,
                binary_ ? IOstream::BINARY : IOstream::ASCII,
                IOstream::currentVersion,
                compressed_ ? IOstream::COMPRESSED : IOstream::UNCOMPRESSED,
                true
            );
    }
}

template<class Type>
void Foam::f4stSetWriter<Type>::closeFile() const
{
    outFilePtr_().closed();
}

template<class Type>
void Foam::f4stSetWriter<Type>::writePoints
(
    const coordSet& points
) const
{
    OFstream& os = outFilePtr_();
//    OFstream os = getOutputFile();
    points.List<point>::writeEntry("points", os);
}

template<class Type>
void Foam::f4stSetWriter<Type>::writeData
(
    const Foam::wordList& valueSetNames,
    const Foam::List<const Foam::Field<Type>*>& valueSets,
    const Foam::scalar& time
) const
{
    OFstream& os = outFilePtr_();
//    OFstream os = getOutputFile();
    word timeName = std::to_string(time);
    forAll(valueSetNames, seti)
    {
        const Field<Type>& field = *valueSets[seti];
        word fieldName = valueSetNames[seti]+"_"+timeName;
        fieldName.erase(fieldName.find_last_not_of('0') + 1, std::string::npos);
        fieldName.erase(fieldName.find_last_not_of('.') + 1, std::string::npos);
        field.writeEntry(fieldName, os);
    }
}

template<class Type>
bool Foam::f4stSetWriter<Type>::meshIsWritten() const
{
    word name_ = "surfaceData";
    char dataFile[80];
    sprintf(dataFile, "f4st/%s.f4", name_.c_str());
    IFstream is(dataFile, binary_ ? IOstream::BINARY : IOstream::ASCII);

    dictionary inputDict(is);
    if (inputDict.found("points"))
    {
        return true;
    }

    return false;
}

template<class Type>
bool Foam::f4stSetWriter<Type>::timeIsWritten(const scalar& timeValue) const
{
    word name_ = "surfaceData";
    char dataFile[80];
    sprintf(dataFile, "f4st/%s.f4", name_.c_str());
    IFstream is(dataFile, binary_ ? IOstream::BINARY : IOstream::ASCII);
    dictionary inputDict(is);
    word currentTime = "Time_" + std::to_string(timeValue);
    if (inputDict.found(currentTime))
    {
        return true;
    }
    return false;
}

template<class Type>
void Foam::f4stSetWriter<Type>::writeTimeValue(const scalar& time) const
{
    OFstream& os = outFilePtr_();

    if (!timeIsWritten(time))
    {
//        word currentTime = "Time_" + std::to_string(time);
        word currentTime = "Time";
        os << currentTime << " " << time << token::END_STATEMENT<< endl;
    }
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::f4stSetWriter<Type>::f4stSetWriter()
:
    writer<Type>(),
    binary_(true),
    compressed_(false),
    meshIsWritten_(false),
    outFilePtr_(nullptr),
    timeList_(0)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::f4stSetWriter<Type>::~f4stSetWriter()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const char* Foam::f4stSetWriter<Type>::ext() const {
    return ".f4";
}

template<class Type>
void Foam::f4stSetWriter<Type>::write
    (
        const coordSet& points,
        const wordList& valueSetNames,
        const List<const Field<Type>*>& valueSets,
        Ostream& os
    ) const
{
    FatalErrorInFunction
        << "Not implemented function write on a given Ostream " << endl
        << exit(FatalError);
}

template<class Type>
void Foam::f4stSetWriter<Type>::write
(
    const coordSet& points,
    const wordList& valueSetNames,
    const List<const Field<Type>*>& valueSets,
    const scalar time
) const
{
    createFile();

    if (!meshIsWritten())
    {
        writePoints(points);
    }

    writeData(valueSetNames, valueSets, time);

    writeTimeValue(time);

    closeFile();
}
// ************************************************************************* //
