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
    (c) 2019-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function1/SelectEntry/SelectEntry.H"


// * * * * * * * * * * * * Private member functions  * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::Select<Type>::checkValid() const
{
    if (childFunctions_.size() != switchValues_.size() + 1)
    {
        FatalErrorInFunction
            << "Input for entry " << this->name_ << " is invalid:" << nl
            << "Number of functions must be one more than number of switch"
            << "values" << nl << exit(FatalError);
    }
    for (label i = 0; i < switchValues_.size()-1; i++)
    {
        if (switchValues_[i] > switchValues_[i+1])
        {
            FatalErrorInFunction
                << "Input for entry " << this->name_ << " is invalid:" << nl
                << "Switch values must be in ascending order" << nl
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::Select<Type>::Select
(
    const word& entryName,
    const dictionary& dict
)
:
    Function1<Type>(entryName),
    switchValues_(),
    childFunctions_()
{
    Istream& is(dict.lookup(entryName));
    word entryType(is);

    is  >> switchValues_;

    PtrList<dictionary> dictList(is);

    childFunctions_.resize(dictList.size());
    forAll(dictList, fi)
    {
        childFunctions_.set(fi, Function1<Type>::New(entryName, dictList[fi]));
    }

    checkValid();
}


template<class Type>
Foam::Function1Types::Select<Type>::Select
(
    const word& entryName,
    const List<scalar>& switchValues,
    const PtrList<Function1<Type>>& childFunctions
)
:
    Function1<Type>(entryName),
    switchValues_(switchValues),
    childFunctions_(childFunctions)
{
    checkValid();
}


template<class Type>
Foam::Function1Types::Select<Type>::Select(const Select& sel)
:
    Function1<Type>(sel),
    switchValues_(sel.switchValues_),
    childFunctions_(sel.childFunctions_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::Select<Type>::~Select()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::Select<Type>::convertTimeBase(const Time& t)
{
    forAll(childFunctions_, i)
    {
        childFunctions_[i].convertTimeBase(t);
    }
    forAll(switchValues_, i)
    {
        switchValues_[i] = t.userTimeToTime(switchValues_[i]);
    }
}


template<class Type>
Type Foam::Function1Types::Select<Type>::value(const scalar x) const
{
    forAll(switchValues_, i)
    {
        if (x < switchValues_[i])
        {
            return childFunctions_[i].value(x);
        }
    }
    return childFunctions_[switchValues_.size()].value(x);
}


template<class Type>
Type Foam::Function1Types::Select<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    Type intx(Zero);

    const scalar xmin = min(x1, x2);
    const scalar xmax = max(x1, x2);

    forAll(childFunctions_, i)
    {
        scalar startVal = xmin;
        scalar endVal = xmax;
        if (i > 0)
        {
            startVal = max(startVal, switchValues_[i-1]);
        }
        if (i < switchValues_.size())
        {
            endVal = min(endVal, switchValues_[i]);
        }
        if (startVal < endVal)
        {
            intx +=
                childFunctions_[i].integrate
                (
                    startVal,
                    endVal
                );
        }
    }

    if (xmin < xmax)
    {
        return intx;
    }
    else
    {
        return -intx;
    }
}


template<class Type>
Type Foam::Function1Types::Select<Type>::integrateYoverX
(
    const scalar x1,
    const scalar x2
) const
{
    Type intx(Zero);

    const scalar xmin = min(x1, x2);
    const scalar xmax = max(x1, x2);

    forAll(childFunctions_, i)
    {
        scalar startVal = xmin;
        scalar endVal = xmax;
        if (i > 0)
        {
            startVal = max(startVal, switchValues_[i-1]);
        }
        if (i < switchValues_.size())
        {
            endVal = min(endVal, switchValues_[i]);
        }
        if (startVal < endVal)
        {
            intx +=
                childFunctions_[i].integrateYoverX
                (
                    startVal,
                    endVal
                );
        }
    }

    if (xmin < xmax)
    {
        return intx;
    }
    else
    {
        return -intx;
    }
}


template<class Type>
Type Foam::Function1Types::Select<Type>::derivative(const scalar x) const
{
    forAll(switchValues_, i)
    {
        if (x < switchValues_[i])
        {
            return childFunctions_[i].derivative(x);
        }
    }
    return childFunctions_[switchValues_.size()].derivative(x);
}


template<class Type>
void Foam::Function1Types::Select<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);

    os  << nl << indent << switchValues_ << nl << indent;
    // TODO: This should be done by writing an actual list of dicts
    os  << token::BEGIN_LIST << nl << indent;
    forAll(childFunctions_, i)
    {
        os << token::BEGIN_BLOCK << nl << indent;
        childFunctions_[i].writeData(os);
        os << token::END_BLOCK << nl << indent;
    }
    os  << token::END_LIST << token::END_STATEMENT << nl;
}


// ************************************************************************* //
