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
    (c) 2020-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/fvOptions/fvOptionList.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "cfdTools/general/solutionControl/pisoControl/pisoControl.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"
#include "cfdTools/general/fvOptions/solveID.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(optionList, 0);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const Foam::dictionary& Foam::fv::optionList::optionsDict
(
    const dictionary& dict
) const
{
    if (dict.found("options"))
    {
        return dict.subDict("options");
    }
    else
    {
        return dict;
    }
}


bool Foam::fv::optionList::readOptions(const dictionary& dict)
{
    checkTimeIndex_ = mesh_.time().timeIndex() + 2;

    bool allOk = true;
    forAll(*this, i)
    {
        option& bs = this->operator[](i);
        bool ok = bs.read(dict.subDict(bs.name()));
        allOk = (allOk && ok);
    }
    return allOk;
}


void Foam::fv::optionList::checkApplied() const
{
    if (mesh_.time().timeIndex() == checkTimeIndex_)
    {
        forAll(*this, i)
        {
            const option& bs = this->operator[](i);
            bs.checkApplied();
        }
    }
}


void Foam::fv::optionList::checkBoundaryApplied() const
{
    if (mesh_.time().timeIndex() == checkTimeIndex_)
    {
        forAll(*this, i)
        {
            const option& bs = this->operator[](i);
            bs.checkBoundaryApplied();
        }
    }
}

const Foam::fvMesh& Foam::fv::optionList::getMesh(const objectRegistry& obr) const
{
    if (isA<fvSolutionRegistry>(obr))
    {
        return dynamic_cast<const fvSolutionRegistry&>(obr).mesh();
    }
    return dynamic_cast<const fvMesh&>(obr);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::optionList::optionList(const objectRegistry& obr, const dictionary& dict)
:
    PtrList<option>(),
    mesh_(getMesh(obr)),
    obr_(obr),
    checkTimeIndex_(mesh_.time().startTimeIndex() + 2)
{
    if (debug)
    {
        Info<< "Constructing optionList with dictionary:" << nl;
        Info<< dict << endl;
    }
    reset(optionsDict(dict), obr);
}


Foam::fv::optionList::optionList(const objectRegistry& obr)
:
    PtrList<option>(),
    mesh_(getMesh(obr)),
    obr_(obr),
    checkTimeIndex_(mesh_.time().startTimeIndex() + 2)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::optionList::reset(const dictionary& dict, const objectRegistry& obr)
{
    // Count number of active fvOptions
    label count = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            count++;
        }
    }

    this->setSize(count);
    label i = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            const word& name = iter().keyword();
            const dictionary& sourceDict = iter().dict();

            this->set
            (
                i++,
                option::New(name, sourceDict, obr)
            );
        }
    }
}


void Foam::fv::optionList::initialise()
{
    // Remove those options that return false
    labelList oldToNew(this->size());
    label newIdx = 0;
    label numRemoved = 0;
    forAll(*this, i)
    {
        if ((*this)[i].initialise())
        {
            oldToNew[i] = newIdx++;
        }
        else
        {
            Info<< "Removing fvOption " << (*this)[i].name() << endl;
            // Move to end for removal
            oldToNew[i] = this->size() - (++numRemoved);
        }
    }
    this->reorder(oldToNew);
    this->setSize(newIdx);

    // Initialise references to other regions
    const objectRegistry& timeDb = mesh_.time();
    HashTable<const solutionInstanceRegistry*> instanceRegistries =
        timeDb.lookupClass<solutionInstanceRegistry>();
    regionOptions_.resize(1);
    regionOptions_.set(0, this);
    regionNames_.resize(1);
    regionNames_[0] = (obr_.name() == word::null ? mesh_.name() : obr_.name());

    // If instancing active, find all regions from active instance
    for (const solutionInstanceRegistry* instance : instanceRegistries)
    {
        if (instance->isActive())
        {
            regionNames_ = instance->regionNames();
            regionOptions_.resize(regionNames_.size());
            forAll(regionNames_, regioni)
            {
                word meshName = instance->whichMesh(regionNames_[regioni]);
                objectRegistry& meshDb =
                    timeDb.lookupObjectRef<objectRegistry>(meshName);
                regionOptions_.set
                (
                    regioni,
                    &meshDb.lookupObjectRef<fvSolutionRegistry>
                    (
                        (regionNames_[regioni] == meshName)
                        ? word::null
                        : regionNames_[regioni]
                    ).registry().lookupObjectRef<optionList>("fvOptions")
                );
            }
            break;
        }
    }
}


bool Foam::fv::optionList::addOption(autoPtr<option>& opt)
{
    if (opt().initialise())
    {
        this->append(opt);
        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::fv::optionList::read(const dictionary& dict)
{
    return readOptions(optionsDict(dict));
}


bool Foam::fv::optionList::writeData(Ostream& os) const
{
    // Write list contents
    forAll(*this, i)
    {
        os  << nl;
        this->operator[](i).writeHeader(os);
        this->operator[](i).writeData(os);
        this->operator[](i).writeFooter(os);
    }

    // Check state of IOstream
    return os.good();
}


void Foam::fv::optionList::addSourceDependencies
(
    SolveTable<solveList>& depList
)
{
    forAll(*this, i)
    {
        option& optionI = this->operator[](i);

        if (optionI.isActive())
        {
            optionI.addSourceDependencies(depList);
        }
    }
}


void Foam::fv::optionList::correct()
{
    forAll(*this, i)
    {
        option& optionI = this->operator[](i);

        if (optionI.isActive())
        {
            if (optionI.execHook() == fv::option::ehtSolve)
            {
                bool finalOuterCorrector = true;
                if (mesh_.foundObject<pimpleControl>(solutionControl::typeName))
                {
                    const pimpleControl& pc =
                        mesh_.lookupObject<pimpleControl>
                        (
                            solutionControl::typeName
                        );
                    if (isA<pisoControl>(pc))
                    {
                        finalOuterCorrector = true;
                    }
                    else
                    {
                        finalOuterCorrector = pc.finalIter();
                    }
                }
                if (finalOuterCorrector)
                {
                    optionI.correct();
                }
            }
            else if (optionI.execHook() == fv::option::ehtOuterCorrect)
            {
                optionI.correct();
            }
        }
    }
}


void Foam::fv::optionList::checkCorrect
(
    const word& funcName,
    const word& fieldName
)
{
    forAll(*this, i)
    {
        option& optionI = this->operator[](i);

        label fieldI = optionI.applyToField(fieldName, thisRegionName());

        if (debug)
        {
            Info<< "calling checkCorrect within function: " << funcName << endl;
            Info<< "given arguments: " << funcName << " and " << fieldName << "\n"
                << "execHookName: " << optionI.execHookName() << endl;
        }
        if (optionI.execHookName() == funcName && fieldI != -1)
        {
            if (debug)
            {
                Info<< "solver hook found: calling solver ..." << endl;
            }

            optionI.setApplied(fieldI);

            if (optionI.isActive())
            {
                optionI.correct();
            }
        }
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const fv::optionList& options)
{
    options.writeData(os);
    return os;
}


// ************************************************************************* //
