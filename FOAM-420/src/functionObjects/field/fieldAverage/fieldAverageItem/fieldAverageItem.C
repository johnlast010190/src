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
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fieldAverage/fieldAverageItem/fieldAverageItem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::functionObjects::fieldAverageItem::EXT_MEAN
(
    "Mean"
);

const Foam::word Foam::functionObjects::fieldAverageItem::EXT_PRIME2MEAN
(
    "Prime2Mean"
);

const Foam::word Foam::functionObjects::fieldAverageItem::EXT_MIN
(
    "Min"
);

const Foam::word Foam::functionObjects::fieldAverageItem::EXT_MAX
(
    "Max"
);

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::fieldAverageItem::baseType,
    2
>::names[] = { "iteration", "time"};

const Foam::NamedEnum
<
    Foam::functionObjects::fieldAverageItem::baseType,
    2
> Foam::functionObjects::fieldAverageItem::baseTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldAverageItem::fieldAverageItem()
:
    active_(false),
    fieldName_("unknown"),
    mean_(0),
    meanFieldName_("unknown"),
    prime2Mean_(0),
    prime2MeanFieldName_("unknown"),
    min_(0),
    minFieldName_("unknown"),
    max_(0),
    maxFieldName_("unknown"),
    base_(ITER),
    window_(-1.0),
    windowName_(""),
    timeIntegral_(0)
{}


Foam::functionObjects::fieldAverageItem::fieldAverageItem
(
    const fieldAverageItem& faItem
)
:
    active_(faItem.active_),
    fieldName_(faItem.fieldName_),
    mean_(faItem.mean_),
    meanFieldName_(faItem.meanFieldName_),
    prime2Mean_(faItem.prime2Mean_),
    prime2MeanFieldName_(faItem.prime2MeanFieldName_),
    min_(faItem.min_),
    minFieldName_(faItem.minFieldName_),
    max_(faItem.max_),
    maxFieldName_(faItem.maxFieldName_),
    base_(faItem.base_),
    window_(faItem.window_),
    windowName_(faItem.windowName_),
    timeIntegral_(faItem.timeIntegral_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldAverageItem::~fieldAverageItem()
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::functionObjects::fieldAverageItem::operator=
(
    const fieldAverageItem& rhs
)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorInFunction
            << "Attempted assignment to self" << nl
            << abort(FatalError);
    }

    // Set updated values
    active_ = rhs.active_;
    fieldName_ = rhs.fieldName_;
    mean_ = rhs.mean_;
    meanFieldName_ = rhs.meanFieldName_;
    prime2Mean_ = rhs.prime2Mean_;
    prime2MeanFieldName_ = rhs.prime2MeanFieldName_;
    min_ = rhs.min_;
    minFieldName_ = rhs.minFieldName_;
    max_ = rhs.max_;
    maxFieldName_ = rhs.maxFieldName_;
    base_ = rhs.base_;
    window_ = rhs.window_;
    windowName_ = rhs.windowName_;
    timeIntegral_ = rhs.timeIntegral_;
}


// ************************************************************************* //
