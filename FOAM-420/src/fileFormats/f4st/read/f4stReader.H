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

Class
    Foam::f4stReader

Description
    Reader for f4st file format

SourceFiles
    f4stReader.C

\*---------------------------------------------------------------------------*/

#ifndef f4stReader_H
#define f4stReader_H

#include "db/IOstreams/Fstreams/OFstream.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/IOstreams/IOstreams/IOstream.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "db/Time/instant/instant.H"
namespace Foam
{

class f4stReader
{
protected:

    static bool isUnique(const scalar& time, const DynamicList<scalar>& timeList);

public:

    static word getFileName();

    static pointField getPoints(IFstream& is);

    template <typename T>
    static Field<T> getField
    (
        const word& keyword,
        IFstream& is,
        const label size
    )
    {
        if (size)
        {
            Field<T> field;
            token keyToken;
            keyType readingWord;
            while (!is.eof())
            {
                is.read(keyToken);
                if (keyToken.isWord())
                {
                    if (keyToken.wordToken() == keyword)
                    {
                        //Skip the next token
                        is.read(keyToken);

                        is >> field;
                        break;
                    }
                }
            }
            return field;
        }
        return Field<T>(0);
    }

    static List<instant> getSampleTimes(const word& pathName);

    //- Construct null
    f4stReader();

};

}

#endif