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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "profileModelList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::profileModelList::profileModelList
(
    const dictionary& dict,
    const bool readFields
)
:
    PtrList<profileModel>(),
    dict_(dict)
{
    if (readFields)
    {
        wordList modelNames(dict.toc());

        Info<< "    Constructing blade profiles:" << endl;

        if (modelNames.size() > 0)
        {
            this->setSize(modelNames.size());

            forAll(modelNames, i)
            {
                const word& modelName = modelNames[i];

                this->set
                (
                    i,
                    profileModel::New(dict.subDict(modelName))
                );
            }
        }
        else
        {
            Info<< "        none" << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::profileModelList::~profileModelList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::profileModelList::connectBlades
(
    const List<word>& names,
    List<label>& addr
) const
{
    // construct the addressing between blade sections and profiles
    forAll(names, bI)
    {
        label index = -1;
        const word& profileName = names[bI];

        forAll(*this, pI)
        {
            const profileModel& pm = this->operator[](pI);

            if (pm.name() == profileName)
            {
                index = pI;
                break;
            }
        }

        if (index == -1)
        {
            List<word> profileNames(size());
            forAll(*this, i)
            {
                const profileModel& pm = this->operator[](i);
                profileNames[i] = pm.name();
            }

            FatalErrorInFunction
                << "Profile " << profileName << " could not be found "
                << "in profile list.  Available profiles are"
                << profileNames << exit(FatalError);
        }
        else
        {
            addr[bI] = index;
        }
    }
}


// ************************************************************************* //
