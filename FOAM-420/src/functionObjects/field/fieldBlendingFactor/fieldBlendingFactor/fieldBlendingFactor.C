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
    (c) 2011 Esi Ltd.

    - with the option of handling surfaceScalar field: Ntua 27.02.2011

\*---------------------------------------------------------------------------*/

#include "fieldBlendingFactor/fieldBlendingFactor/fieldBlendingFactor.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldBlendingFactor, 0);
    addToRunTimeSelectionTable(functionObject, fieldBlendingFactor, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldBlendingFactor::fieldBlendingFactor
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    lastUpdateIndex_(-1),
    staticBlend_(),
    fieldBlendingItems_(0)
{
    read(dict);

    //update field blending factors
    if (lastUpdateIndex_ != obr().time().timeIndex())
    {
        forAll(fieldBlendingItems_, fdi)
        {
            fieldBlendingItems_[fdi].update(&staticBlend_->field());
        }
        lastUpdateIndex_ = obr().time().timeIndex();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldBlendingFactor::~fieldBlendingFactor()
{
/*

    //store a list of blending field names so we can delete unused ones

    autoPtr<IOdictionary> oldBlendingFieldNames
    (
        new IOdictionary
        (
            IOobject
            (
                "oldBlendingFieldNames",
                obr_.time().system(),
                obr_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dictionary()
        )
    );

    wordList fieldNames(fieldBlendingItems_.size());

    forAll(fieldBlendingItems_, i)
    {
        fieldNames[i] = fieldBlendingItems_[i].name();
    }

    oldBlendingFieldNames->add("blendingFieldNames", fieldNames);

    Info<< "Old blending field names: " << oldBlendingFieldNames
         << endl;

    oldBlendingFieldNames->store(oldBlendingFieldNames);
*/
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldBlendingFactor::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    Log << type() << " " << name() <<  " read:" << nl;

    //read static blending item
    word staticBlendName(name()+"StaticBase");

    staticBlend_.reset
    (
        new fieldBlendingItem
        (
            obr(),
            mesh_,
            dict,
            staticBlendName
        )
    );
    staticBlend_->update();

    //read dynamic blending items

    if (dict.found("fieldBlending"))
    {
        PtrList<entry> fieldDicts(dict.lookup("fieldBlending"));

        fieldBlendingItems_.setSize
        (
            fieldDicts.size()
        );

        forAll(fieldDicts, fdi)
        {
            const entry& fieldDict = fieldDicts[fdi];
            word bfName
                = fieldDict.dict().lookupOrDefault<word>
                ("name", word(fieldDict.keyword() + "BlendingFactor"));

            fieldBlendingItems_.set
            (
                fdi,
                new fieldBlendingItem
                (
                    obr(),
                    mesh_,
                    fieldDict.dict(),
                    bfName
                )
            );
        }
    }
    return true;
}


bool Foam::functionObjects::fieldBlendingFactor::execute()
{
    //dont execute when read just occured
    if (lastUpdateIndex_ != obr().time().timeIndex())
    {
        Log << type() << " " << name() <<  " execute:" << nl;
        forAll(fieldBlendingItems_, fdi)
        {
            fieldBlendingItems_[fdi].update(&staticBlend_->field());
        }
        lastUpdateIndex_ = obr().time().timeIndex();
        Log << endl;
    }

    return true;
}


bool Foam::functionObjects::fieldBlendingFactor::write()
{
    Log << type() << " " << name() <<  " write:" << nl;
    forAll(fieldBlendingItems_,fbi)
    {
        fieldBlendingItems_[fbi].write();
    }
    Log << endl;

    return true;
}


bool Foam::functionObjects::fieldBlendingFactor::end()
{
    fieldBlendingItems_.clear();
    staticBlend_.clear();

    return true;
}

// ************************************************************************* //
