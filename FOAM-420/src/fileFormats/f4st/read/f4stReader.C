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

#include "f4st/read/f4stReader.H"

bool Foam::f4stReader::isUnique
(
    const Foam::scalar& time, const Foam::DynamicList<Foam::scalar>& timeList
)
{
    bool unique = true;
    forAll(timeList, tI)
    {
        if (timeList[tI] == time)
        {
            unique = false;
        }
    }
    return unique;
}

Foam::word Foam::f4stReader::getFileName()
{
    word name_ = "surfaceData";
    char dataFile[80];
    sprintf(dataFile, "f4st/%s.f4", name_.c_str());
    return dataFile;
}

Foam::pointField Foam::f4stReader::getPoints
(
    Foam::IFstream& is
)
{
    DynamicList<point> points;
    token keyToken;
    keyType keyword;
    Field<point> pointL;
    while (!is.eof())
    {
        is.read(keyToken);
        if (keyToken.isWord())
        {
            keyword = keyToken.wordToken();
            if (keyword == "points")
            {
                is >> pointL;
                break;
            }
        }
    }
    return  pointL;
}

Foam::List<Foam::instant> Foam::f4stReader::getSampleTimes
(
    const word& pathName
)
{
    token keyToken;
    keyType keyword;
    IFstream is(pathName, IOstream::BINARY);
    DynamicList<scalar> timeList;
    while (!is.eof())
    {
        is.read(keyToken);

        // If the token is a valid keyword set 'keyword' return true...
        if (keyToken.isWord())
        {
            keyword = keyToken.wordToken();
        }
        else if (keyToken.isString())
        {
            // Enable wildcards
            keyword = keyToken.stringToken();
        }

        if (keyword=="Time")
        {
            string timeLine;
            is.getLine(timeLine);
            scalar timeValue = std::stod(timeLine);
            if (isUnique(timeValue, timeList))
            {
                timeList.append(timeValue);
            }
        }
        token nextToken(is);
        is.putBack(nextToken);
    }
    List<instant> sampleTimes(timeList.size());
    forAll(sampleTimes, ti)
    {
        instant tInst(timeList[ti]);
        sampleTimes[ti]= tInst;
    }
    return sampleTimes;
}