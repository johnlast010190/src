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
    (c) 2022-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "coordinate/systems/frameData.H"
#include "coordinate/systems/coordinateSystem.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* NamedEnum
    <
        frameData::relationType,
        3
    >::names[] =
    {
        "global",
        "parent",
        "frame"
    };
}

const Foam::NamedEnum<Foam::frameData::relationType, 3>
Foam::frameData::relationTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::frameData::frameData()
:
    definedInFrame_(rtGlobal),
    relatedFrame_(word::null),
    isLocal_(false),
    isLocalOnRead_(false)
{}


Foam::frameData::frameData(const frameData& fd)
:
    definedInFrame_(fd.definedInFrame_),
    relatedFrame_(fd.relatedFrame_),
    isLocal_(fd.isLocal_),
    isLocalOnRead_(fd.isLocalOnRead_)
{}


Foam::frameData::frameData(const dictionary& dict)
:
    definedInFrame_(rtGlobal),
    relatedFrame_(word::null),
    isLocal_(false),
    isLocalOnRead_(false)
{
    frameData::read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::frameData::transformToLocal
(
    coordinateSystem& csys,
    const coordinateSystem& csysParent
)
{
    if (!csys.coordinateTransforms().isLocal())
    {
        vector origin = csys.origin();
        vector e1 = csys.e1();
        vector e2 = csys.e2();

        origin = csysParent.localPosition(origin);
        e1 = csysParent.localVector(e1);
        e2 = csysParent.localVector(e2);

        csys.update(origin, e1, e2);
    }
}


void Foam::frameData::transformToGlobal
(
    coordinateSystem& csys,
    const coordinateSystem& csysParent
)
{
    if (csys.coordinateTransforms().isLocal())
    {
        vector origin = csys.origin();
        vector e1 = csys.e1();
        vector e2 = csys.e2();

        origin = csysParent.globalPosition(origin);
        e1 = csysParent.globalVector(e1);
        e2 = csysParent.globalVector(e2);

        csys.update(origin, e1, e2);
    }
}


Foam::autoPtr<Foam::coordinateSystem>
Foam::frameData::frameToGlobal
(
    const dictionary& dict,
    const objectRegistry& obr
)
{
    const word dictName = "meshObjects";

    autoPtr<IOobject> meshObjectIO
    (
        new IOobject
        (
            dictName,
            obr.time().caseSystem(),
            obr,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    if (!meshObjectIO->typeHeaderOk<localIOdictionary>())
    {
        // Fallback: look for definition at the top level
        autoPtr<IOobject> defaultRegionIO
        (
            new IOobject
            (
                dictName,
                obr.time().caseSystem(),
                obr.time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        if (defaultRegionIO->typeHeaderOk<localIOdictionary>())
        {
            meshObjectIO = defaultRegionIO;
        }
    }

    IFstream is(meshObjectIO->objectPath());

    // Dict was not found skipping the localisation
    // needed for frames constructed directly from the dict with no reference
    // to a meshObjects file.
    if (!is.good())
    {
        return coordinateSystem::New(dict, coordinateSystem::typeName);
    }

    dictionary meshObjectsDict(is);
    const word frameName = dict.dictName();
    PtrList<coordinateSystem> coordinateSystems;
    wordList parentFrameNames;
    const dictionary* parentFrameDict = &dict;

    label maxIter = 0;
    while (maxIter < 100)
    {
        ++maxIter;
        coordinateSystems.append
        (
            coordinateSystem::New
            (
                (*parentFrameDict),
                coordinateSystem::typeName
            )
        );

        if
        (
            parentFrameDict->found("parentFrameName")
         && coordinateSystems.last().coordinateTransforms().isLocal()
        )
        {
            parentFrameNames.append(parentFrameDict->lookup("parentFrameName"));
            parentFrameDict =
                meshObjectsDict.subDictPtr(parentFrameNames.last());
        }
        else
        {
            break;
        }
    }

    for (label i = coordinateSystems.size() - 2; i >= 0; --i)
    {
        transformToGlobal
        (
            coordinateSystems[i],
            coordinateSystems[i + 1]
        );
    }
    coordinateSystems.first().coordinateTransforms().setLocal(false);
    if
    (
       !dict.found("parentFrameName")
     && coordinateSystems.first().coordinateTransforms().isLocalOnRead()
    )
    {
        WarningInFunction
            << "No parent frame found for frame " << frameName
            << " Local specification has no effect." << endl;
    }

    return autoPtr<coordinateSystem>(coordinateSystems.first().clone());
}


bool Foam::frameData::read(const dictionary& dict)
{
    if (dict.found("definedInFrame"))
    {
        definedInFrame_ =
            relationTypeNames_.read(dict.lookup("definedInFrame"));
    }
    if (dict.found("originallyLocal"))
    {
        isLocalOnRead_ = dict.lookup("originallyLocal");
    }
    else
    {
        switch (definedInFrame_)
        {
            case rtGlobal:
            {
                isLocalOnRead_ = false;
                break;
            }
            case rtParent:
            {
                isLocalOnRead_ = true;
                break;
            }
            case rtFrame:
            {
                FatalErrorInFunction
                    << "Frame data defined in frame not supported yet" << nl
                    << abort(FatalError);
                break;
            }
        }
    }

    switch (definedInFrame_)
    {
        case rtGlobal:
        {
            isLocal_ = false;
            break;
        }
        case rtParent:
        {
            isLocal_ = true;
            break;
        }
        case rtFrame:
        {
            FatalErrorInFunction
                << "Frame data defined in frame not supported yet" << nl
                << abort(FatalError);

            // isLocal_ = false;
            // isLocalOnRead_ = false;
            // relatedFrame_ = word(dict.lookup("referenceFrame"));
            break;
        }
    }

    return true;
}


void Foam::frameData::write(Ostream& os) const
{
    if (isLocal_)
    {
        os.writeEntry("definedInFrame", "parent");
    }
    if (isLocalOnRead_)
    {
        os.writeEntry("originallyLocal", true);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::frameData::operator=(const frameData& fd)
{
    definedInFrame_ = fd.definedInFrame_;
    relatedFrame_ = fd.relatedFrame_;
    isLocal_ = fd.isLocal_;
    isLocalOnRead_ = fd.isLocalOnRead_;
}

// ************************************************************************* //
