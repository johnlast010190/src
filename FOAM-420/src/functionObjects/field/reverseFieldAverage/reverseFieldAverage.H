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
    (c) 2010-2012 Esi Ltd.
    (c) 2015 OpenCFD Ltd.
    (c) 2017 Esi Ltd.

Class
    Foam::functionObjects::reverseFieldAverage

Group
    grpFieldFunctionObjects

Description
    This function object calculates average quantities for a user-specified
    selection of volumetric and surface fields.  Fields are entered as a list
    of sub-dictionaries, which indicate the type of averages to perform, and
    can be updated during the calculation.  The current options include:
    - \c mean: arithmetic mean:
        \f[
            \overline{x} = \frac{1}{N}\displaystyle\sum\limits_{i=0}^N x_i
        \f]
    - \c prime2Mean: prime-squared mean
        \f[
            \overline{x'}^2 = \frac{1}{N}\displaystyle\sum\limits_{i=0}^N
            (x_i - \overline{x})^2
        \f]
    - base: average over 'time', or 'iteration' (\f$N\f$ in the above)
    - window: optional averaging window, specified in 'base' units

    Average field names are constructed by concatenating the base field with
    the averaging type, e.g. when averaging field 'U', the resultant fields
    are:
    - arithmetic mean field, UMean
    - prime-squared field, UPrime2Mean

    Information regarding the number of averaging steps, and total averaging
    time are written on a per-field basis to the
    \c "<functionObject name>Properties" dictionary, located in \<time\>/uniform

    When restarting form a previous calculation, the averaging is continuous or
    may be restarted using the \c restartOnRestart option.

    The averaging process may be restarted after each calculation output time
    using the \c restartOnOutput option or restarted periodically using the \c
    periodicRestart option and setting \c restartPeriod to the required
    averaging period.

    Example of function object specification:
    \verbatim
    reverseFieldAverage1
    {
        type reverseFieldAverage;
        functionObjectLibs ("libfieldFunctionObjects.so");
        ...
        restartOnRestart    false;
        restartOnOutput     false;
        periodicRestart     false;
        restartPeriod       0.002;
        fields
        (
            U
            {
                mean            on;
                prime2Mean      on;
                base            time;
                window          10.0;
                windowName      w1;
            }
            p
            {
                mean            on;
                prime2Mean      on;
                base            time;
            }
        );
    }
    \endverbatim

    \heading Function object usage
    \table
        Property        | Description           | Required    | Default value
        type            | Type name: reverseFieldAverage | yes       |
        restartOnRestart| Restart the averaging on restart | no | no
        restartOnOutput | Restart the averaging on output | no | no
        periodicRestart | Periodically restart the averaging | no | no
        restartPeriod   | Periodic restart period | conditional |
        log             | Log to standard output  | no        | yes
        fields          | List of fields and averaging options | yes |
    \endtable


Note
    To employ the \c prime2Mean option, the \c mean option must be selecetd.

SeeAlso
    Foam::functionObjects::functionObject
    Foam::functionObjects::OutputFilterFunctionObject

SourceFiles
    reverseFieldAverage.C
    reverseFieldAverageTemplates.C
    reverseFieldAverageItem.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_reverseFieldAverage_H
#define functionObjects_reverseFieldAverage_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "fields/volFields/volFieldsFwd.H"
#include "primitives/bools/Switch/Switch.H"
#include "./reverseFieldAverageItem/reverseFieldAverageItem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class reverseFieldAverage Declaration
\*---------------------------------------------------------------------------*/

class reverseFieldAverage
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        //- Time at last call, prevents repeated averaging
        label prevTimeIndex_;

        //- Restart the averaging process on restart
        Switch restartOnRestart_;

        //- Restart the averaging process on output
        Switch restartOnOutput_;

        //- Periodically restart the averaging process
        Switch periodicRestart_;

        //- Restart period
        scalar restartPeriod_;

        //- Time to start averaging
        scalar averagingStartTime_;

        //- Initialised flag
        bool initialised_;

        //- List of field average items, describing what averages to be
        //  calculated and output
        List<reverseFieldAverageItem> faItems_;

        // Counters

            //- Iteration steps counter
            List<label> totalIter_;

            //- Total time counter
            List<scalar> totalTime_;

            //- Index for periodic restart
            label periodIndex_;


    // Private Member Functions

        // Initialisation routines

            //- Checkout fields (causes deletion) from the database
            //  and reset lists
            void resetFields();

            //- Reset lists (clear existing values) and initialize averaging.
            //  Check requested field averages are valid, populate field lists
            void initialize();

            //- Restart averaging for restartOnOutput
            void restart();

            //- Add mean average field to database
            template<class Type>
            void addMeanFieldType(const label fieldI);

            //- Add mean average field to database
            template<class Type>
            void addMeanField(const label fieldI);

            //- Add prime-squared average field to database
            template<class Type1, class Type2>
            void addPrime2MeanFieldType(const label fieldI);

            //- Add prime-squared average field to database
            template<class Type1, class Type2>
            void addPrime2MeanField(const label fieldI);


        // Calculation functions

            //- Main calculation routine
            virtual void calcAverages();

            //- Calculate mean average fields
            template<class Type>
            void calculateMeanFieldType(const label fieldI) const;

            //- Calculate mean average fields
            template<class Type>
            void calculateMeanFields() const;

            //- Calculate prime-squared average fields
            template<class Type1, class Type2>
            void calculatePrime2MeanFieldType(const label fieldI) const;

            //- Calculate prime-squared average fields
            template<class Type1, class Type2>
            void calculatePrime2MeanFields() const;

            //- Add mean-squared field value to prime-squared mean field
            template<class Type1, class Type2>
            void addMeanSqrToPrime2MeanType(const label fieldI) const;

            //- Add mean-squared field value to prime-squared mean field
            template<class Type1, class Type2>
            void addMeanSqrToPrime2Mean() const;


        // I-O

            //- Write averages
            virtual void writeAverages() const;

            //- Write fields
            template<class Type>
            void writeFieldType(const word& fieldName) const;

            //- Write fields
            template<class Type>
            void writeFields() const;

            //- Write averaging properties - steps and time
            void writeAveragingProperties();

            //- Read averaging properties - steps and time
            void readAveragingProperties();


        //- Disallow default bitwise copy construct
        reverseFieldAverage(const reverseFieldAverage&);

        //- Disallow default bitwise assignment
        void operator=(const reverseFieldAverage&);


public:

    //- Runtime type information
    TypeName("reverseFieldAverage");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        reverseFieldAverage
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~reverseFieldAverage();


    // Member Functions

        //- Read the field average data
        virtual bool read(const dictionary&);

        //- Execute the averaging
        virtual bool execute();

        //- Calculate the field average data and write
        virtual bool write();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "reverseFieldAverage/reverseFieldAverageTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
