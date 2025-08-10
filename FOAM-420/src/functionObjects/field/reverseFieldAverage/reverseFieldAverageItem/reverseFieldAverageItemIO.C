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
    (c) 2016 OpenCFD Ltd.
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "reverseFieldAverage/reverseFieldAverageItem/reverseFieldAverageItem.H"
#include "db/IOstreams/IOstreams.H"
#include "db/dictionary/dictionaryEntry/dictionaryEntry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::reverseFieldAverageItem::reverseFieldAverageItem(Istream& is)
:
    active_(false),
    fieldName_("unknown"),
    mean_(0),
    meanFieldName_("unknown"),
    prime2Mean_(0),
    prime2MeanFieldName_("unknown"),
    base_(ITER),
    window_(-1.0)
{
    is.check(FUNCTION_NAME);

    const dictionaryEntry entry(dictionary::null, is);

    fieldName_ = entry.keyword();

    mean_ = entry.lookupOrDefault<Switch>("mean", true);
    prime2Mean_ = entry.lookupOrDefault<Switch>("prime2Mean", false);
    base_ = baseTypeNames_[entry.lookupOrDefault<word>("base", "time")];

    window_ = entry.lookupOrDefault<scalar>("window", -1.0);
    windowName_ = entry.lookupOrDefault<word>("windowName", "");

    meanFieldName_ = fieldName_ + EXT_MEAN;
    prime2MeanFieldName_ = fieldName_ + EXT_PRIME2MEAN;
    if ((window_ > 0) && (windowName_ != ""))
    {
        meanFieldName_ = meanFieldName_ + "_" + windowName_;
        prime2MeanFieldName_ = prime2MeanFieldName_ + "_" + windowName_;
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::functionObjects::operator>>
(
    Istream& is, reverseFieldAverageItem& faItem
)
{
    is.check(FUNCTION_NAME);

    const dictionaryEntry entry(dictionary::null, is);

    faItem.active_ = false;
    faItem.fieldName_ = entry.keyword();
    entry.lookup("mean") >> faItem.mean_;
    entry.lookup("prime2Mean") >> faItem.prime2Mean_;
    faItem.base_ = faItem.baseTypeNames_[entry.lookup("base")];
    faItem.window_ = entry.lookupOrDefault<scalar>("window", -1.0);
    faItem.windowName_ = entry.lookupOrDefault<word>("windowName", "");

    faItem.meanFieldName_ = faItem.fieldName_ + reverseFieldAverageItem::EXT_MEAN;
    faItem.prime2MeanFieldName_ =
        faItem.fieldName_ + reverseFieldAverageItem::EXT_PRIME2MEAN;

    if ((faItem.window_ > 0) && (faItem.windowName_ != ""))
    {
        faItem.meanFieldName_ =
            faItem.meanFieldName_ + "_" + faItem.windowName_;

        faItem.prime2MeanFieldName_ =
            faItem.prime2MeanFieldName_ + "_" + faItem.windowName_;
    }
    return is;
}


Foam::Ostream& Foam::functionObjects::operator<<
(
    Ostream& os, const reverseFieldAverageItem& faItem
)
{
    os.check(FUNCTION_NAME);

    os  << faItem.fieldName_ << nl << token::BEGIN_BLOCK << nl;
    os.writeEntry("mean", faItem.mean_);
    os.writeEntry("prime2Mean", faItem.prime2Mean_);
    os.writeEntry("base", faItem.baseTypeNames_[faItem.base_]);

    if (faItem.window_ > 0)
    {
        os.writeEntry("window", faItem.window_);

        if (faItem.windowName_ != "")
        {
            os.writeEntry("windowName", faItem.windowName_);
        }
    }

    os  << token::END_BLOCK << nl;

    os.check(FUNCTION_NAME);

    return os;
}


// ************************************************************************* //
