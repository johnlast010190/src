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
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "modificationSwitches.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::modificationSwitches::modificationSwitches()
:
    resetInternalFields_(true),
    resetBoundaryFields_(true),
    resetSystemDicts_(true),
    resetConstDicts_(true),
    resetUniformDicts_(true),
    resetBoundaryMesh_(true),
    deleteUnusedFields_(true),
    strictPatchNameChecking_(true),
    reuseExistingDictionaries_(false),
    printBoundaries_(false),
    cleanInstanceFields_(false)
{}


Foam::modificationSwitches::modificationSwitches
(
    const dictionary& dict
)
:
    resetInternalFields_(true),
    resetBoundaryFields_(true),
    resetSystemDicts_(true),
    resetConstDicts_(true),
    resetUniformDicts_(true),
    resetBoundaryMesh_(true),
    deleteUnusedFields_(true),
    strictPatchNameChecking_(true),
    reuseExistingDictionaries_(false),
    printBoundaries_(false),
    cleanInstanceFields_(false)
{
    merge(dict);
}

Foam::modificationSwitches::modificationSwitches
(
    const dictionary& dict,
    const modificationSwitches& ms
)
:
    resetInternalFields_(ms.resetInternalFields_),
    resetBoundaryFields_(ms.resetBoundaryFields_),
    resetSystemDicts_(ms.resetSystemDicts_),
    resetConstDicts_(ms.resetConstDicts_),
    resetUniformDicts_(ms.resetUniformDicts_),
    resetBoundaryMesh_(ms.resetBoundaryMesh_),
    deleteUnusedFields_(ms.deleteUnusedFields_),
    strictPatchNameChecking_(ms.strictPatchNameChecking_),
    reuseExistingDictionaries_(ms.reuseExistingDictionaries_),
    printBoundaries_(ms.printBoundaries_),
    cleanInstanceFields_(ms.cleanInstanceFields_)
{
    merge(dict);
}


Foam::modificationSwitches::modificationSwitches
(
    const modificationSwitches& ms
)
:
    resetInternalFields_(ms.resetInternalFields_),
    resetBoundaryFields_(ms.resetBoundaryFields_),
    resetSystemDicts_(ms.resetSystemDicts_),
    resetConstDicts_(ms.resetConstDicts_),
    resetUniformDicts_(ms.resetUniformDicts_),
    resetBoundaryMesh_(ms.resetBoundaryMesh_),
    deleteUnusedFields_(ms.deleteUnusedFields_),
    strictPatchNameChecking_(ms.strictPatchNameChecking_),
    reuseExistingDictionaries_(ms.reuseExistingDictionaries_),
    printBoundaries_(ms.printBoundaries_),
    cleanInstanceFields_(ms.cleanInstanceFields_)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::modificationSwitches::~modificationSwitches()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::modificationSwitches::merge(const dictionary& dict)
{
    if (dict.found("modificationSwitches"))
    {
        const dictionary& msd(dict.subDict("modificationSwitches"));

        resetInternalFields_
            = msd.lookupOrDefault<Switch>
            ("resetInternalFields", resetInternalFields_);
        resetBoundaryFields_
            = msd.lookupOrDefault<Switch>
            ("resetBoundaryFields", resetBoundaryFields_);
        resetSystemDicts_
            = msd.lookupOrDefault<Switch>
            ("resetSystemDicts", resetSystemDicts_);
        resetConstDicts_
            = msd.lookupOrDefault<Switch>
            ("resetConstDicts", resetConstDicts_);
        resetUniformDicts_
            = msd.lookupOrDefault<Switch>
            ("resetUniformDicts", resetUniformDicts_);
        resetBoundaryMesh_
            = msd.lookupOrDefault<Switch>
            ("resetBoundaryMesh", resetBoundaryMesh_);
        deleteUnusedFields_
            = msd.lookupOrDefault<Switch>
            ("deleteUnusedFields", deleteUnusedFields_);
        strictPatchNameChecking_
            = msd.lookupOrDefault<Switch>
            ("strictPatchNameChecking", strictPatchNameChecking_);
        reuseExistingDictionaries_
            = msd.lookupOrDefault<Switch>
            ("reuseExistingDictionaries", reuseExistingDictionaries_);
        printBoundaries_
            = msd.lookupOrDefault<Switch>
            ("printBoundaries", printBoundaries_);
        cleanInstanceFields_
            = msd.lookupOrDefault<Switch>
            ("cleanInstanceFields", cleanInstanceFields_);
    }
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::modificationSwitches::operator=(const modificationSwitches& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorInFunction
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    resetInternalFields_ = rhs.resetInternalFields_;
    resetBoundaryFields_ = rhs.resetBoundaryFields_;
    resetSystemDicts_ = rhs.resetSystemDicts_;
    resetConstDicts_ = rhs.resetConstDicts_;
    resetUniformDicts_ = rhs.resetUniformDicts_;
    resetBoundaryMesh_ = rhs.resetBoundaryMesh_;
    deleteUnusedFields_ = rhs.deleteUnusedFields_;
    strictPatchNameChecking_ = rhs.strictPatchNameChecking_;
    reuseExistingDictionaries_ = rhs.reuseExistingDictionaries_;
    printBoundaries_ = rhs.printBoundaries_;
    cleanInstanceFields_ = rhs.cleanInstanceFields_;
}

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
