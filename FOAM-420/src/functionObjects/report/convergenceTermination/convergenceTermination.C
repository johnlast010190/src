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
    (c) 2013-2013 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "convergenceTermination/convergenceTermination.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::convergenceTermination::actionType,
        3
    >::names[] =
    {
        "noWriteNow",
        "writeNow",
        "nextWrite"
    };
}


const Foam::NamedEnum<Foam::convergenceTermination::actionType, 3>
    Foam::convergenceTermination::actionTypeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * /

void Foam::convergenceTermination::read(const dictionary& dict)
{
    if (dict.found("action"))
    {
        action_ = actionTypeNames_.read(dict.lookup("action"));
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::convergenceTermination::convergenceTermination
(
    const objectRegistry& obr,
    const dictionary& dict
)
:
    obr_(obr),
    fieldName_(dict.dictName()),
    action_(writeNow),
    stdDevTol_(dict.lookupOrDefault<scalar>("stdDevTol", -GREAT)),
    sampleSize_(dict.lookupOrDefault<scalar>("sampleSize", 100))
{
    read(dict);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::convergenceTermination::calculate(const scalar value)
{
    storedData_.push(value);

    while
    (
        storedData_.size() > sampleSize_
    )
    {
        storedData_.pop();
    }

    if (storedData_.size() > 0 && storedData_.size() == sampleSize_)
    {
        scalar meanValue = 0.;
        forAllConstIter(FIFOStack<scalar>, storedData_, iter)
        {
            meanValue += iter();
        }
        meanValue /= storedData_.size();

        scalar stdDev = 0.;
        forAllConstIter(FIFOStack<scalar>, storedData_, iter)
        {
            stdDev += sqr(iter()-meanValue);
        }
        stdDev = sqrt(stdDev/storedData_.size());

        Info<<"Convergence check for field: "<<fieldName_
            <<" stdDev : "<< stdDev << endl;

        if (stdDev <= stdDevTol_)
        {
            switch (action_)
            {
                case noWriteNow :
                {
                    if (obr_.time().stopAt(Time::saNoWriteNow))
                    {
                        Info<<"Stopping without writing data"
                            <<" Standard deviation : "<<stdDev
                            <<" for field : "<<fieldName_
                            <<endl;
                    }
                    break;
                }

                case writeNow :
                {
                    if (obr_.time().stopAt(Time::saWriteNow))
                    {
                        Info<<"Stopping and writing data."
                            <<" Standard deviation : "<<stdDev
                            <<" for field : "<<fieldName_
                            <<endl;
                    }
                    break;
                }

                case nextWrite :
                {
                    if (obr_.time().stopAt(Time::saNextWrite))
                    {
                        Info<<"Stopping after next data write."
                            <<" Standard deviation : "<<stdDev
                            <<" for field : "<<fieldName_
                            <<endl;
                    }
                    break;
                }
            }
            return true;
        }
    }
    return false;
}
