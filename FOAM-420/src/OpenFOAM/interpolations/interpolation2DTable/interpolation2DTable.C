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
    (c) 2017-2023 Esi Ltd.
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "interpolations/interpolationTable/tableReaders/openFoam/openFoamTableReader.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Type>
void Foam::interpolation2DTable<Type>::readTable()
{
    fileName fName(fileName_);
    fName.expand();

    // Read data from file
    reader_()(fName, *this);

    if (this->empty())
    {
        FatalErrorInFunction
            << "table read from " << fName << " is empty" << nl
            << exit(FatalError);
    }

    // Check that the data are in ascending order
    checkOrder();

    checkTableShapeAndValues();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolation2DTable<Type>::interpolation2DTable()
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>(),
    boundsHandling_(interpolation2DTable::WARN),
    fileName_("fileNameIsUndefined"),
    reader_(nullptr)
{}


template<class Type>
Foam::interpolation2DTable<Type>::interpolation2DTable
(
    const List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>& values,
    const boundsHandling bounds,
    const fileName& fName
)
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>(values),
    boundsHandling_(bounds),
    fileName_(fName),
    reader_(nullptr)
{}


template<class Type>
Foam::interpolation2DTable<Type>::interpolation2DTable(const fileName& fName)
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>(),
    boundsHandling_(interpolation2DTable::WARN),
    fileName_(fName),
    reader_(new openFoamTableReader<Type>(dictionary()))
{
    readTable();
}


template<class Type>
Foam::interpolation2DTable<Type>::interpolation2DTable(const dictionary& dict)
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>(),
    boundsHandling_(wordToBoundsHandling(dict.lookup("outOfBounds"))),
    fileName_(dict.lookup("file")),
    reader_(tableReader<Type>::New(dict))
{
    readTable();
}


template<class Type>
Foam::interpolation2DTable<Type>::interpolation2DTable
(
     const interpolation2DTable& interpTable
)
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>(interpTable),
    boundsHandling_(interpTable.boundsHandling_),
    fileName_(interpTable.fileName_),
    reader_(interpTable.reader_)    // note: steals reader. Used in write().
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::interpolation2DTable<Type>::interpolateValue
(
    const List<Tuple2<scalar, Type>>& data,
    const scalar lookupValue
) const
{
    const scalar minLimit = data.first().first();
    const scalar maxLimit = data.last().first();

    if (lookupValue < minLimit)
    {
        switch (boundsHandling_)
        {
            case interpolation2DTable::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") less than lower "
                    << "bound (" << minLimit << ")" << nl
                    << exit(FatalError);
                break;
            }
            case interpolation2DTable::WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") less than lower "
                    << "bound (" << minLimit << ")" << nl
                    << "    Continuing with the first entry"
                    << endl;
                // Behaviour as per 'CLAMP'
                return data.first().second();
                break;
            }
            case interpolation2DTable::CLAMP:
            {
                return data.first().second();
                break;
            }
        }
    }
    else if (lookupValue >= maxLimit)
    {
        switch (boundsHandling_)
        {
            case interpolation2DTable::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") greater than upper "
                    << "bound (" << maxLimit << ")" << nl
                    << exit(FatalError);
                break;
            }
            case interpolation2DTable::WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") greater than upper "
                    << "bound (" << maxLimit << ")" << nl
                    << "    Continuing with the last entry"
                    << endl;
                // Behaviour as per 'CLAMP'
                return data.last().second();
                break;
            }
            case interpolation2DTable::CLAMP:
            {
                return data.last().second();
                break;
            }
        }
    }

    // look for the correct range in X
    label lo = 0;
    label hi = 0;

    forAll(data, i)
    {
        if (lookupValue >= data[i].first())
        {
            lo = hi = i;
        }
        else
        {
            hi = i;
            break;
        }
    }

    if (lo == hi)
    {
        return data[lo].second();
    }
    else
    {
        Type m =
            (data[hi].second() - data[lo].second())
           /(data[hi].first() - data[lo].first());

        // normal interpolation
        return data[lo].second() + m*(lookupValue - data[lo].first());
    }
}


template<class Type>
template<class BinaryOp>
Foam::label Foam::interpolation2DTable<Type>::Xi
(
    const BinaryOp& bop,
    const scalar x,
    const bool reverse
) const
{
    const table& Table = *this;

    label limitI = 0;
    if (reverse)
    {
        limitI = Table.size() - 1;
    }

    if (bop(x, Table[limitI].first()))
    {
        switch (boundsHandling_)
        {
            case interpolation2DTable::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << x << ") out of bounds"
                    << exit(FatalError);
                break;
            }
            case interpolation2DTable::WARN:
            {
                WarningInFunction
                    << "value (" << x << ") out of bounds"
                    << endl;
                // Behaviour as per 'CLAMP'
                return limitI;
                break;
            }
            case interpolation2DTable::CLAMP:
            {
                return limitI;
                break;
            }
            default:
            {
                FatalErrorInFunction
                    << "Unhandled enumeration " << boundsHandling_
                    << abort(FatalError);
            }
        }
    }

    label i = 0;
    if (reverse)
    {
        label sizeX = Table.size();
        i = 0;
        while ((i < sizeX) && (x > Table[i].first()))
        {
            i++;
        }
    }
    else
    {
        i = Table.size() - 1;
        while ((i > 0) && (x < Table[i].first()))
        {
            i--;
        }
    }

    return i;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
Type Foam::interpolation2DTable<Type>::operator()
(
    const scalar x,
    const scalar y
) const
{
    // Considers all of the list in Y being equal
    const label sizeX = this->size();

    const table& Table = *this;

    if (sizeX == 0)
    {
        WarningInFunction
            << "cannot interpolate a zero-sized table - returning zero"
            << endl;

        return Zero;
    }
    else if (sizeX == 1)
    {
        // only 1 column (in X) - interpolate to find Y value
        return interpolateValue(Table.first().second(), y);
    }
    else
    {
        // have 2-D data, interpolate

        // find low and high indices in the X range that bound x
        const label X0i = Xi(lessOp<scalar>(), x, false);
        const label X1i = Xi(greaterOp<scalar>(), x, true);

        if (X0i == X1i)
        {
            return interpolateValue(Table[X0i].second(), y);
        }
        else
        {
            Type y0(interpolateValue(Table[X0i].second(), y));
            Type y1(interpolateValue(Table[X1i].second(), y));

            // gradient in X
            const scalar x0 = Table[X0i].first();
            const scalar x1 = Table[X1i].first();
            Type mX = (y1 - y0)/(x1 - x0);

            // interpolate
            return y0 + mX*(x - x0);
        }
    }
}


template<class Type>
inline void Foam::interpolation2DTable<Type>::operator=
(
    const interpolation2DTable<Type>& interpolationTable
)
{
    table& Table = *this;
    const table& inputTable = interpolationTable;
    Table.setSize(inputTable.size());
    forAll(inputTable, i)
    {
        Table[i] = inputTable[i];
    }
}


template<class Type>
inline void Foam::interpolation2DTable<Type>::scale
(
    const Foam::scalar& Scale
)
{
    table& Table = *this;
    forAll(Table, i)
    {
        forAll(Table[i].second(), j)
        {
            Table[i].second()[j].second() *= Scale;
        }
    }
}


template<class Type>
inline void Foam::interpolation2DTable<Type>::idaxpy
(
    const scalar& a,
    const scalar& b,
    const interpolation2DTable<Type>& Table2
)
{
    // construct non-const reference to this data
    table& Table = *this;

    // perform y = a*y + b*x operation and interpolate from x
    forAll(Table, i)
    {
        forAll(Table[i].second(), j)
        {
            const scalar X = Table[i].first();
            const scalar Y = Table[i].second()[j].first();

            Table[i].second()[j].second() =
                a*Table[i].second()[j].second() + b*Table2(X, Y);
        }
    }
}


template<class Type>
inline Type Foam::interpolation2DTable<Type>::integrateOverY
(
    const scalar x,
    const scalar y0,
    const scalar y1
) const
{
    if (y1 < y0)
    {
        FatalErrorInFunction
            << " Integration has reversed integration limits:"
            << " y0 = " << y0 << " y1 = " << y1 << "."
            << nl << fileName_ << nl
            << exit(FatalError);
    }
    else if (y0 == y1)
    {
        return Zero;
    }

    const List<Tuple2<scalar, Type>>& tableY = (*this).first().second();
    scalarField Y(tableY.size(), 0.0);
    Field<Type> functionOfY(tableY.size(), Zero);
    forAll(tableY, i)
    {
        Y[i] = tableY[i].first();
        functionOfY[i] = (*this)(x, Y[i]);
    }

    Type integral = Zero;

    // Handle clamp values
    if (y0 < Y.first())
    {
        if
        (
            boundsHandling_ == interpolation2DTable::WARN
         || boundsHandling_ == interpolation2DTable::CLAMP
        )
        {
            const scalar dy = (y1 < Y.first()) ? (y1 - y0) : (Y.first() - y0);
            integral += functionOfY.first()*dy;
        }
        else
        {
            FatalErrorInFunction
                << "Table doesn't have clamp or warn definition."
                << "Integration outside with value x = " << x
                << " y0 = " << y0 << " y1 = " << y1 << " is outside of range."
                << nl << fileName_ << nl
                << exit(FatalError);
        }
    }

    if (y1 > Y.last())
    {
        if
        (
            boundsHandling_ == interpolation2DTable::WARN
         || boundsHandling_ == interpolation2DTable::CLAMP
        )
        {
            const scalar dy = (y0 > Y.last()) ? (y1 - y0) : (y1 - Y.last());
            integral += functionOfY.last()*dy;
        }
        else
        {
            FatalErrorInFunction
                << "Table doesn't have clamp or warn definition."
                << " Integration outside with value x = " << x
                << " y0 = " << y0 << " y1 = " << y1 << " is outside of range."
                << nl << fileName_ << nl
                << exit(FatalError);
        }
    }

    // Trapezoidal integration
    for (label i = 0; i < (Y.size() - 1); ++i)
    {
        if (y0 <= Y[i + 1] && y1 <= Y[i + 1] && y1 >= Y[i] && y0 >= Y[i])
        {
            integral +=
                (y1 - y0)
               *((*this)(x, y1) + (*this)(x, y0))*0.5;
        }
        else if (y0 <= Y[i] && y1 >= Y[i] && y1 <= Y[i + 1])
        {
            integral +=
                (y1 - Y[i])
               *((*this)(x, y1) + functionOfY[i])*0.5;
        }
        else if (y0 >= Y[i] && y0 <= Y[i + 1] && y1 >= Y[i + 1])
        {
            integral +=
                (Y[i + 1] - y0)
               *((*this)(x, y0) + functionOfY[i+1])*0.5;
        }
        else if (y0 <= Y[i] && y1 >= Y[i + 1])
        {
            integral +=
                (Y[i + 1] - Y[i])
               *(functionOfY[i + 1] + functionOfY[i])*0.5;
        }
    }

    return integral;
}


template<class Type>
Foam::word Foam::interpolation2DTable<Type>::boundsHandlingToWord
(
     const boundsHandling& bound
) const
{
    word enumName("warn");

    switch (bound)
    {
        case interpolation2DTable::ERROR:
        {
            enumName = "error";
            break;
        }
        case interpolation2DTable::WARN:
        {
            enumName = "warn";
            break;
        }
        case interpolation2DTable::CLAMP:
        {
            enumName = "clamp";
            break;
        }
    }

    return enumName;
}


template<class Type>
typename Foam::interpolation2DTable<Type>::boundsHandling
Foam::interpolation2DTable<Type>::wordToBoundsHandling
(
    const word& bound
)
{
    if (bound == "error")
    {
        return interpolation2DTable::ERROR;
    }
    else if (bound == "warn")
    {
        return interpolation2DTable::WARN;
    }
    else if (bound == "clamp")
    {
        return interpolation2DTable::CLAMP;
    }
    else
    {
        WarningInFunction
            << "bad outOfBounds specifier " << bound
            << " using 'warn'" << endl;

        return interpolation2DTable::WARN;
    }
}


template<class Type>
typename Foam::interpolation2DTable<Type>::boundsHandling
Foam::interpolation2DTable<Type>::outOfBounds
(
    const boundsHandling& bound
)
{
    boundsHandling prev = boundsHandling_;
    boundsHandling_ = bound;
    return prev;
}


template<class Type>
void Foam::interpolation2DTable<Type>::checkOrder() const
{
    const table& Table = *this;

    scalar prevValue = Table.first().first();

    for (label i=1; i < this->size(); ++i)
    {
        const scalar currValue = Table[i].first();

        // avoid duplicate values (divide-by-zero error)
        if (currValue <= prevValue)
        {
            FatalErrorInFunction
                << "out-of-order value: "
                << currValue << " at index " << i << nl
                << exit(FatalError);
        }
        prevValue = currValue;
    }
}


template<class Type>
void Foam::interpolation2DTable<Type>::checkTableShapeAndValues() const
{
    const table& Table = *this;
    const List<Tuple2<scalar, Type>>& subTable1 = Table.first().second();
    // Check that they have the same size (form rectangle/squeare)
    forAll(Table, i)
    {
        if (subTable1.size() != Table[i].second().size())
        {
            FatalErrorInFunction
                << "The table doesn't form square/rectangle." << nl
                << fileName_
                << exit(FatalError);
        }
    }

    forAll(subTable1, i)
    {
        const scalar value1 = subTable1[i].first();
        for (label j = 1; j < Table.size(); ++j)
        {
            const List<Tuple2<scalar, Type>>& subTablei = Table[j].second();
            const scalar value2 = subTablei[i].first();
            if (value1 != value2)
            {
                FatalErrorInFunction
                    << "The values " << value1 << " " << value2
                    << " in the 2D interpolation table are not the same." << nl
                    << fileName_
                    << exit(FatalError);
            }
        }
    }
}


template<class Type>
void Foam::interpolation2DTable<Type>::write(Ostream& os) const
{
    os.writeEntry("file", fileName_);
    os.writeEntry("outOfBounds", boundsHandlingToWord(boundsHandling_));

    os  << *this;
}


// ************************************************************************* //
