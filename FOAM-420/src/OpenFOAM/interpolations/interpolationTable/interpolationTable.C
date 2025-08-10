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
    (c) 2010-2022 Esi Ltd.
    (c) 2017-2020 OpenCFD Ltd.
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "interpolations/interpolationTable/interpolationTable.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "interpolations/interpolationTable/tableReaders/openFoam/openFoamTableReader.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolationTable<Type>::interpolationTable()
:
    List<Tuple2<scalar, Type>>(),
    boundsHandling_(interpolationTable::CLAMP),
    fileNamePtr_(new fileName("fileNameIsUndefined")),
    reader_(nullptr)
{}


template<class Type>
Foam::interpolationTable<Type>::interpolationTable
(
    const List<Tuple2<scalar, Type>>& values,
    const boundsHandling bounds,
    const fileName& fName
)
:
    List<Tuple2<scalar, Type>>(values),
    boundsHandling_(bounds),
    fileNamePtr_(),
    reader_(nullptr)
{
    if (fName.size())
    {
        fileNamePtr_.reset(new fileName(fName));
    }
}


template<class Type>
Foam::interpolationTable<Type>::interpolationTable(const fileName& fName)
:
    List<Tuple2<scalar, Type>>(),
    boundsHandling_(interpolationTable::CLAMP),
    fileNamePtr_(new fileName(fName)),
    reader_(new openFoamTableReader<Type>(dictionary()))
{
    readTable();
}


template<class Type>
Foam::interpolationTable<Type>::interpolationTable(const dictionary& dict)
:
    List<Tuple2<scalar, Type>>(),
    boundsHandling_
    (
        wordToBoundsHandling
        (
            dict.lookupOrDefault
            (
                "outOfBounds",
                word("warn")
            )
        )
    ),
    fileNamePtr_(),
    reader_(tableReader<Type>::New(dict))
{

    if (dict.found("file"))
    {
        fileNamePtr_.reset(new fileName(dict.lookup("file")));
        readTable();
    }
    else if (dict.found("fileName"))
    {
        fileNamePtr_.reset(new fileName(dict.lookup("fileName")));
        readTable();
    }
    else //data must be present
    {
        List<Tuple2<scalar, Type>>::operator=
            (List<Tuple2<scalar, Type>>(dict.lookup("data")));
        check();
    }

}


template<class Type>
Foam::interpolationTable<Type>::interpolationTable
(
     const interpolationTable& interpTable
)
:
    List<Tuple2<scalar, Type>>(interpTable),
    boundsHandling_(interpTable.boundsHandling_),
    fileNamePtr_(),
    reader_(interpTable.reader_)    // note: steals reader. Used in write()
{
    if (interpTable.fileNamePtr_.valid())
    {
        fileNamePtr_.reset(new fileName(interpTable.fileNamePtr_()));
    }
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::word Foam::interpolationTable<Type>::boundsHandlingToWord
(
     const boundsHandling& bound
) const
{
    word enumName("warn");

    switch (bound)
    {
        case interpolationTable::ERROR:
        {
            enumName = "error";
            break;
        }
        case interpolationTable::WARN:
        {
            enumName = "warn";
            break;
        }
        case interpolationTable::CLAMP:
        {
            enumName = "clamp";
            break;
        }
        case interpolationTable::REPEAT:
        {
            enumName = "repeat";
            break;
        }
    }

    return enumName;
}


template<class Type>
typename Foam::interpolationTable<Type>::boundsHandling
Foam::interpolationTable<Type>::wordToBoundsHandling
(
    const word& bound
) const
{
    if (bound == "error")
    {
        return interpolationTable::ERROR;
    }
    else if (bound == "warn")
    {
        return interpolationTable::WARN;
    }
    else if (bound == "clamp")
    {
        return interpolationTable::CLAMP;
    }
    else if (bound == "repeat")
    {
        return interpolationTable::REPEAT;
    }
    else
    {
        WarningInFunction
            << "bad outOfBounds specifier " << bound << " using 'warn'" << endl;

        return interpolationTable::WARN;
    }
}


template<class Type>
typename Foam::interpolationTable<Type>::boundsHandling
Foam::interpolationTable<Type>::outOfBounds
(
    const boundsHandling& bound
)
{
    boundsHandling prev = boundsHandling_;
    boundsHandling_ = bound;
    return prev;
}


template<class Type>
void Foam::interpolationTable<Type>::readTable()
{
    if (fileNamePtr_.valid())
    {
        // preserve the original (unexpanded) fileName to avoid absolute paths
        // appearing subsequently in the write() method
        fileName fName(fileNamePtr_());

        fName.expand();

        // Read data from file
        IFstream(fName).operator()() >> *this;

        // Check that the data are okay
        check();
    }
}


template<class Type>
void Foam::interpolationTable<Type>::check() const
{
    const label n = this->size();
    scalar prevValue = this->first().first();

    for (label i=1; i<n; ++i)
    {
        const scalar currValue =
            List<Tuple2<scalar, Type>>::operator[](i).first();

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

    if (this->empty())
    {
        FatalErrorInFunction
            << "table is empty" << nl
            << exit(FatalError);
    }
}


template<class Type>
void Foam::interpolationTable<Type>::write(Ostream& os) const
{
    if (fileNamePtr_.valid())
    {
        os.writeEntry("fileName", fileNamePtr_());
    }
    else
    {
        os.beginBlock("data");
        forAll(*this, ei)
        {
            os << indent << this->operator[](ei) << nl;
        }
        os.endBlock();
    }
    os.writeEntry("outOfBounds", boundsHandlingToWord(boundsHandling_));

    if (reader_.valid())
    {
        reader_->write(os);
    }
}


template<class Type>
Type Foam::interpolationTable<Type>::rateOfChange(const scalar value) const
{
    label n = this->size();

    if (n <= 1)
    {
        // There are not enough entries to provide a rate of change
        return 0;
    }

    const scalar minLimit = this->first().first();
    const scalar maxLimit = this->last().first();
    scalar lookupValue = value;

    if (lookupValue < minLimit)
    {
        switch (boundsHandling_)
        {
            case interpolationTable::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") underflow" << nl
                    << exit(FatalError);
                break;
            }
            case interpolationTable::WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") underflow" << nl
                    << "    Zero rate of change."
                    << endl;
                // behaviour as per 'CLAMP'
                return 0;
                break;
            }
            case interpolationTable::CLAMP:
            {
                return 0;
                break;
            }
            case interpolationTable::REPEAT:
            {
                // adjust lookupValue to >= minLimit
                scalar span = maxLimit-minLimit;
                lookupValue = fmod(lookupValue-minLimit, span) + minLimit;
                break;
            }
        }
    }
    else if (lookupValue >= maxLimit)
    {
        switch (boundsHandling_)
        {
            case interpolationTable::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") overflow" << nl
                    << exit(FatalError);
                break;
            }
            case interpolationTable::WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") overflow" << nl
                    << "    Zero rate of change."
                    << endl;
                // Behaviour as per 'CLAMP'
                return 0;
                break;
            }
            case interpolationTable::CLAMP:
            {
                return 0;
                break;
            }
            case interpolationTable::REPEAT:
            {
                // adjust lookupValue <= maxLimit
                scalar span = maxLimit-minLimit;
                lookupValue = fmod(lookupValue-minLimit, span) + minLimit;
                break;
            }
        }
    }

    label lo = 0;
    label hi = 0;

    // look for the correct range
    for (label i = 0; i < n; ++i)
    {
        if (lookupValue >= List<Tuple2<scalar, Type>>::operator[](i).first())
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
        // we are at the end of the table - or there is only a single entry
        return 0;
    }
    else if (hi == 0)
    {
        // this treatment should only occur under these conditions:
        //  -> the 'REPEAT' treatment
        //  -> (0 <= value <= minLimit)
        //  -> minLimit > 0
        // Use the value at maxLimit as the value for value=0
        lo = n - 1;

        return
        (
            (
                List<Tuple2<scalar, Type>>::operator[](hi).second()
              - List<Tuple2<scalar, Type>>::operator[](lo).second()
            )
           /(
               List<Tuple2<scalar, Type>>::operator[](hi).first()
             + minLimit
             - List<Tuple2<scalar, Type>>::operator[](lo).first()
            )
        );
    }
    else
    {
        // normal rate of change
        return
        (
            (
                List<Tuple2<scalar, Type>>::operator[](hi).second()
              - List<Tuple2<scalar, Type>>::operator[](lo).second()
            )
           /(
                List<Tuple2<scalar, Type>>::operator[](hi).first()
              - List<Tuple2<scalar, Type>>::operator[](lo).first()
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


template<class Type>
void Foam::interpolationTable<Type>::operator=
(
    const interpolationTable<Type>& rhs
)
{
    if (this == &rhs)
    {
        return;
    }

    typedef Tuple2<scalar, Type> value_type;
    static_cast<List<value_type>&>(*this) = rhs;
    boundsHandling_ = rhs.boundsHandling_;
    if (rhs.fileNamePtr_.valid())
    {
        fileNamePtr_.reset(new fileName(rhs.fileNamePtr_()));
    }
    reader_.reset(rhs.reader_.clone());
}



template<class Type>
const Foam::Tuple2<Foam::scalar, Type>&
Foam::interpolationTable<Type>::operator[](const label i) const
{
    label ii = i;
    label n  = this->size();

    if (n <= 1)
    {
        ii = 0;
    }
    else if (ii < 0)
    {
        switch (boundsHandling_)
        {
            case interpolationTable::ERROR:
            {
                FatalErrorInFunction
                    << "index (" << ii << ") underflow" << nl
                    << exit(FatalError);
                break;
            }
            case interpolationTable::WARN:
            {
                WarningInFunction
                    << "index (" << ii << ") underflow" << nl
                    << "    Continuing with the first entry"
                    << endl;
                // Behaviour as per 'CLAMP'
                ii = 0;
                break;
            }
            case interpolationTable::CLAMP:
            {
                ii = 0;
                break;
            }
            case interpolationTable::REPEAT:
            {
                while (ii < 0)
                {
                    ii += n;
                }
                break;
            }
        }
    }
    else if (ii >= n)
    {
        switch (boundsHandling_)
        {
            case interpolationTable::ERROR:
            {
                FatalErrorInFunction
                    << "index (" << ii << ") overflow" << nl
                    << exit(FatalError);
                break;
            }
            case interpolationTable::WARN:
            {
                WarningInFunction
                    << "index (" << ii << ") overflow" << nl
                    << "    Continuing with the last entry"
                    << endl;
                // Behaviour as per 'CLAMP'
                ii = n - 1;
                break;
            }
            case interpolationTable::CLAMP:
            {
                ii = n - 1;
                break;
            }
            case interpolationTable::REPEAT:
            {
                while (ii >= n)
                {
                    ii -= n;
                }
                break;
            }
        }
    }

    return List<Tuple2<scalar, Type>>::operator[](ii);
}


template<class Type>
Type Foam::interpolationTable<Type>::operator()(const scalar value) const
{
    label n = this->size();

    if (n <= 1)
    {
        return this->first().second();
    }

    const scalar minLimit = this->first().first();
    const scalar maxLimit = this->last().first();
    scalar lookupValue = value;

    if (lookupValue < minLimit)
    {
        switch (boundsHandling_)
        {
            case interpolationTable::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") underflow" << nl
                    << exit(FatalError);
                break;
            }
            case interpolationTable::WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") underflow" << nl
                    << "    Continuing with the first entry"
                    << endl;
                // Behaviour as per 'CLAMP'
                return this->first().second();
                break;
            }
            case interpolationTable::CLAMP:
            {
                return this->first().second();
                break;
            }
            case interpolationTable::REPEAT:
            {
                // adjust lookupValue to >= minLimit
                const scalar span = maxLimit-minLimit;
                lookupValue = fmod(lookupValue-minLimit, span) + minLimit;
                break;
            }
        }
    }
    else if (lookupValue >= maxLimit)
    {
        switch (boundsHandling_)
        {
            case interpolationTable::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") overflow" << nl
                    << exit(FatalError);
                break;
            }
            case interpolationTable::WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") overflow" << nl
                    << "    Continuing with the last entry"
                    << endl;
                // Behaviour as per 'CLAMP'
                return this->last().second();
                break;
            }
            case interpolationTable::CLAMP:
            {
                return this->last().second();
                break;
            }
            case interpolationTable::REPEAT:
            {
                // adjust lookupValue <= maxLimit
                const scalar span = maxLimit-minLimit;
                lookupValue = fmod(lookupValue-minLimit, span) + minLimit;
                break;
            }
        }
    }

    label lo = 0;
    label hi = 0;

    // look for the correct range
    for (label i = 0; i < n; ++i)
    {
        if (lookupValue >= List<Tuple2<scalar, Type>>::operator[](i).first())
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
        // we are at the end of the table - or there is only a single entry
        return List<Tuple2<scalar, Type>>::operator[](hi).second();
    }
    else if (hi == 0)
    {
        // this treatment should only occur under these conditions:
        //  -> the 'REPEAT' treatment
        //  -> (0 <= value <= minLimit)
        //  -> minLimit > 0
        // Use the value at maxLimit as the value for value=0
        lo = n - 1;

        return
        (
            List<Tuple2<scalar, Type>>::operator[](lo).second()
          + (
                List<Tuple2<scalar, Type>>::operator[](hi).second()
              - List<Tuple2<scalar, Type>>::operator[](lo).second()
            )
           *(lookupValue / minLimit)
        );
    }
    else
    {
        // normal interpolation
        return
        (
            List<Tuple2<scalar, Type>>::operator[](lo).second()
          + (
                List<Tuple2<scalar, Type>>::operator[](hi).second()
              - List<Tuple2<scalar, Type>>::operator[](lo).second()
            )
           *(
                lookupValue
              - List<Tuple2<scalar, Type>>::operator[](lo).first()
            )
           /(
                List<Tuple2<scalar, Type>>::operator[](hi).first()
              - List<Tuple2<scalar, Type>>::operator[](lo).first()
            )
        );
    }
}


// ************************************************************************* //
