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
    (c) 2022-2023 Esi Ltd.

Class
    Foam::functionObjects::fractionalEffectiveDose

Group
    grpFieldFunctionObjects

Description
    Calculation and reporting of fractional effective dose (FED)
    of incapacitation, for fire propagation simulations.

    The Fractional Effective Dose index (FED) was developed by
    D. A. Purser and is a commonly used measure of human incapacitation
    due to the combustion gases. The following equation is used
    to calculate the FED value:

    \f[
        FED_{IN} = (FED_{CO} \times HV_{CO_2}) + FED_{O_2}
    \f]

    in which the quantities on the right-hand side depend on
    species concentrations, the person's level of activity
    and some empirical coefficients of this model.

    Reference:
    \verbatim
        Purser, D. A. (2016).
        Handbook of Fire Protection Engineering.
        Springer, 5th edition.
    \endverbatim

    The FED calculation process may be reset periodically using the
    \c periodicReset option and setting \c resetPeriod to the required
    resetting period, or may be subjected to a user-defined one-time reset,
    through the setting \c resetTime.

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::volRegion
    Foam::functionObjects::writeFile

SourceFiles
    fractionalEffectiveDose.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fractionalEffectiveDose_H
#define functionObjects_fractionalEffectiveDose_H

#include "fieldExpression/fieldExpression.H"
#include "functionObjects/volRegion/volRegion.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class fractionalEffectiveDose Declaration
\*---------------------------------------------------------------------------*/

class fractionalEffectiveDose
:
    public fieldExpression,
    public volRegion,
    public writeFile
{
public:

    // Public Data

        //- List of input species names
        static wordList speciesFields;

        //- Activity type enumeration
        enum class activityLevel
        {
            activityRest,     //!< At rest
            activityLight,    //!< Light work
            activityHeavy     //!< Heavy work
        };

        //- Activity type names
        static const Enum<activityLevel> activityLevelNames;

private:

    // Private Data

        //- Average quantities

            //- Fractional effective dose of incapacitation [ppm x min]
            dimensionedScalar avgFedIN_;

            //- Fraction of an incapacitating dose of
            //  carbon monoxide [ppm x min]
            dimensionedScalar avgFedCO_;

            //- Fraction of an incapacitating dose due to
            //  oxygen depletion [ppm x min]
            dimensionedScalar avgFedO2_;

            //- Hyperventilation frequency factor [-]
            dimensionedScalar avgHVCO2_;

        //- Field of carbon monoxide FED
        volScalarField fedCO_;

        //- Field of oxygen depletion FED
        volScalarField fedO2_;

        //- Activity level for FED coefficients
        activityLevel activityLevel_;

        //- FED coefficient
        scalar FEDcoeff_;

        //- Time at last call
        scalar prevTime_;

        //- Time index at last call, prevents repeated integration on time
        label prevTimeIndex_;

        //- Periodically reset the FED cummulative value
        Switch periodicReset_;

        //- Index for periodic reset
        label periodIndex_;

        //- Reset period
        scalar resetPeriod_;

        //- Specific reset time
        scalar resetTime_;

        //- Initialised flag
        bool initialised_;

        //- Names of the relevant fields to be written
        wordList FEDFieldNames_;

    // Private Member Functions

        //- Calculate the fractional effective dose field and averages,
        //  returning true if successful
        virtual bool calc() override;

        //- Main driver function
        void calculateFED();

        //- Update field averages
        void updateAverages();

        //- Reset calculations
        void reset();

        //- Filter a field according to cellIDs
        template<typename Type>
        tmp<Field<Type>> filterField(const Field<Type>& field) const;

        //- Output file header information
        void writeFileHeader(Ostream& os) const;

        //- Disallow default bitwise copy construct
        fractionalEffectiveDose(const fractionalEffectiveDose&);

        //- Disallow default bitwise assignment
        void operator=(const fractionalEffectiveDose&);

public:

    //- Runtime type information
    TypeName("fractionalEffectiveDose");


    // Constructors

        //- Construct from Time and dictionary
        fractionalEffectiveDose
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~fractionalEffectiveDose();


    // Member Functions

        //- Read the data
        virtual bool read(const dictionary&) override;

        //- Write results
        virtual bool write() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
