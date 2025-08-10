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
    (c) 2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "fieldValues/surfaceFieldValue/surfaceFieldValue.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

inline const Foam::labelList&
Foam::functionObjects::fieldValues::surfaceFieldValue::faceId() const
{
    return faceId_;
}


inline const Foam::labelList&
Foam::functionObjects::fieldValues::surfaceFieldValue::facePatch() const
{
    return facePatchId_;
}


inline const Foam::boolList&
Foam::functionObjects::fieldValues::surfaceFieldValue::faceFlip() const
{
    return faceFlip_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::functionObjects::fieldValues::surfaceFieldValue::regionTypes&
Foam::functionObjects::fieldValues::surfaceFieldValue::regionType() const
{
    return regionType_;
}


inline Foam::fileName
Foam::functionObjects::fieldValues::surfaceFieldValue::outputDir() const
{
    return baseFileDir()/name()/"surface"/time_.timeName();
}


// ************************************************************************* //
