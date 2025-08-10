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
    (c) 2016-2017 OpenFOAM Foundation

Class
    Foam::chemistryTabulationMethod

Description
    An abstract class for chemistry tabulation.

SourceFiles
    chemistryTabulationMethod.C

\*---------------------------------------------------------------------------*/

#ifndef chemistryTabulationMethod_H
#define chemistryTabulationMethod_H

#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "fields/Fields/scalarField/scalarField.H"
#include "primitives/bools/Switch/Switch.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
class TDACChemistryModel;

/*---------------------------------------------------------------------------*\
                           Class chemistryTabulationMethod Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class chemistryTabulationMethod
{

protected:

    const dictionary& dict_;

    const dictionary coeffsDict_;

    //- Is tabulation active?
    Switch active_;

    //- Switch to select performance logging
    Switch log_;

    TDACChemistryModel<CompType, ThermoType>& chemistry_;

    scalar tolerance_;


public:

    //- Runtime type information
    TypeName("chemistryTabulationMethod");


    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        chemistryTabulationMethod,
        dictionary,
        (
            const dictionary& dict,
            TDACChemistryModel<CompType, ThermoType>& chemistry
        ),
        (dict, chemistry)
    );


    // Constructors

        //- Construct from components
        chemistryTabulationMethod
        (
            const dictionary& dict,
            TDACChemistryModel<CompType, ThermoType>& chemistry
        );


    // Selectors

        static autoPtr<chemistryTabulationMethod> New
        (
            const IOdictionary& dict,
            TDACChemistryModel<CompType, ThermoType>& chemistry
        );


    //- Destructor
    virtual ~chemistryTabulationMethod();


    // Member Functions

        inline bool active()
        {
            return active_;
        }

        inline bool log()
        {
            return active_ && log_;
        }

        inline bool variableTimeStep()
        {
            return chemistry_.variableTimeStep();
        }

        inline scalar tolerance() const
        {
            return tolerance_;
        }

        virtual label size() = 0;

        virtual void writePerformance() = 0;

        // Retrieve function: (only virtual here)
        // Try to retrieve a stored point close enough (according to tolerance)
        // to a stored point. If successful, it returns true and store the
        // results in RphiQ, i.e. the result of the integration of phiQ
        virtual bool retrieve
        (
             const scalarField& phiQ,
             scalarField& RphiQ
        ) = 0;

        // Add function: (only virtual here)
        // Add information to the tabulation algorithm. Give the reference for
        // future retrieve (phiQ) and the corresponding result (RphiQ).
        virtual label add
        (
            const scalarField& phiQ,
            const scalarField& RphiQ,
            const scalar rho,
            const scalar deltaT
        ) = 0;

        // Update function: (only virtual here)
        // The underlying structure of the tabulation is updated/cleaned
        // to increase the performance of the retrieve
        virtual bool update() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "chemistryModel/TDACChemistryModel/tabulation/chemistryTabulationMethod/chemistryTabulationMethod.C"
    #include "chemistryModel/TDACChemistryModel/tabulation/chemistryTabulationMethod/chemistryTabulationMethodNew.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
