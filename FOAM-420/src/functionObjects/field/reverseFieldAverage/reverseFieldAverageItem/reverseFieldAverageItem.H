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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

Class
    Foam::reverseFieldAverageItem

Description
    Helper class to describe what form of averaging to apply.  A set will be
    applied to each base field in Foam::reverseFieldAverage, of the form:

    \verbatim
    {
        mean            on;
        prime2Mean      on;
        base            time; // iteration
        window          200;  // optional averaging window
        windowName      w1;   // optional window name (default = "")
    }
    \endverbatim

    The averaging window corresponds to the averaging interval (iters or time)
    If not specified, the averaging is over 'all iters/time'

SourceFiles
    reverseFieldAverageItem.C
    reverseFieldAverageItemIO.C

\*---------------------------------------------------------------------------*/

#ifndef reverseFieldAverageItem_H
#define reverseFieldAverageItem_H

#include "primitives/enums/NamedEnum.H"
#include "primitives/bools/Switch/Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

namespace functionObjects
{

// Forward declaration of friend functions and operators
class reverseFieldAverageItem;
Istream& operator>>(Istream&, reverseFieldAverageItem&);
Ostream& operator<<(Ostream&, const reverseFieldAverageItem&);


/*---------------------------------------------------------------------------*\
                      Class reverseFieldAverageItem Declaration
\*---------------------------------------------------------------------------*/

class reverseFieldAverageItem
{
public:

    // Public data

        // File and field name extensions

            //- Mean average
            static const word EXT_MEAN;

            //- Prime-squared average
            static const word EXT_PRIME2MEAN;

        //- Enumeration defining the averaging base type
        enum baseType
        {
            ITER,
            TIME
        };


private:

    // Private data

        //- Active flag
        Switch active_;

        //- Field name
        word fieldName_;

        //- Compute mean flag
        Switch mean_;

        //- Name of mean field
        word meanFieldName_;

        //- Compute prime-squared mean flag
        Switch prime2Mean_;

        //- Name of prime-squared mean field
        word prime2MeanFieldName_;

        //- Averaging base type names
        static const NamedEnum<baseType, 2> baseTypeNames_;

        //- Averaging base type
        baseType base_;

        //- Averaging window - defaults to -1 for 'all iters/time'
        scalar window_;

        //- Averaging window name - defaults to 'window'
        word windowName_;


public:

    // Constructors

        //- Construct null
        reverseFieldAverageItem();

        //- Construct from Istream
        reverseFieldAverageItem(Istream&);

        //- Construct as copy
        reverseFieldAverageItem(const reverseFieldAverageItem&);


    //- Destructor
    ~reverseFieldAverageItem();


    // Member Functions

        // Access

            //- Return const access to the active flag
            const Switch& active() const
            {
                return active_;
            }

            //- Return non-const access to the active flag
            Switch& active()
            {
                return active_;
            }

            //- Return const access to the field name
            const word& fieldName() const
            {
                return fieldName_;
            }

            //- Return const access to the mean flag
            const Switch& mean() const
            {
                return mean_;
            }

            //- Return non-const access to the mean flag
            Switch& mean()
            {
                return mean_;
            }

            //- Return const access to the mean field name
            const word& meanFieldName() const
            {
                return meanFieldName_;
            }

            //- Return const access to the prime-squared mean flag
            const Switch& prime2Mean() const
            {
                return prime2Mean_;
            }

            //- Return non-const access to the prime-squared mean flag
            Switch& prime2Mean()
            {
                return prime2Mean_;
            }

            //- Return const access to the prime-squared mean field name
            const word& prime2MeanFieldName() const
            {
                return prime2MeanFieldName_;
            }

            //- Return averaging base type name
            const word base() const
            {
                return baseTypeNames_[base_];
            }

            //- Return true if base is ITER
            Switch iterBase() const
            {
                return base_ == ITER;
            }

            //- Return true if base is time
            Switch timeBase() const
            {
                return base_ == TIME;
            }

            scalar window() const
            {
                return window_;
            }

            const word& windowName() const
            {
                return windowName_;
            }


    // Modify

            //- Return const access to the mean field name
            void setMeanFieldName(const word& name)
            {
                meanFieldName_ = name;
            }

            //- Return const access to the prime-squared mean field name
            void setPrime2MeanFieldName(const word& name)
            {
                prime2MeanFieldName_ = name;
            }

    // Member Operators

        void operator=(const reverseFieldAverageItem&);


    // Friend Operators

        friend bool operator==
        (
            const reverseFieldAverageItem& a,
            const reverseFieldAverageItem& b
        )
        {
            return
                a.fieldName_ == b.fieldName_
             && a.mean_ == b.mean_
             && a.meanFieldName_ == b.meanFieldName_
             && a.prime2Mean_ == b.prime2Mean_
             && a.prime2MeanFieldName_ == b.prime2MeanFieldName_
             && a.base_ == b.base_
             && a.window_ == b.window_
             && a.windowName_ == b.windowName_;
        }

        friend bool operator!=
        (
            const reverseFieldAverageItem& a,
            const reverseFieldAverageItem& b
        )
        {
            return !(a == b);
        }


    // IOstream Operators

        friend Istream& operator>>(Istream&, reverseFieldAverageItem&);
        friend Ostream& operator<<(Ostream&, const reverseFieldAverageItem&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
