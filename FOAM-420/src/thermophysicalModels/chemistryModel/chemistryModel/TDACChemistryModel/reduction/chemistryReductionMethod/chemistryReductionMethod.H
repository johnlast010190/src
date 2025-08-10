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
    (c) 2016 OpenFOAM Foundation

Class
    Foam::chemistryReductionMethod

Description
    An abstract class for methods of chemical mechanism reduction

SourceFiles
    chemistryReductionMethod.C

\*---------------------------------------------------------------------------*/

#ifndef chemistryReductionMethod_H
#define chemistryReductionMethod_H

#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "primitives/bools/Switch/Switch.H"
#include "fields/Fields/scalarField/scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class CompType, class ThermoType>
class TDACChemistryModel;

/*---------------------------------------------------------------------------*\
                           Class chemistrySolver Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class chemistryReductionMethod
{

protected:

    const IOdictionary& dict_;

    //- Dictionary that store the algorithm data
    const dictionary coeffsDict_;

    //- Is mechanism reduction active?
    Switch active_;

    //- Switch to select performance logging
    Switch log_;

    TDACChemistryModel<CompType, ThermoType>& chemistry_;

    //- List of active species (active = true)
    List<bool> activeSpecies_;

    //- Number of active species
    label NsSimp_;

    //- Number of species
    const label nSpecie_;

    //- Tolerance for the mechanism reduction algorithm
    scalar tolerance_;


public:

    //- Runtime type information
    TypeName("chemistryReductionMethod");


    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        chemistryReductionMethod,
        dictionary,
        (
            const IOdictionary& dict,
            TDACChemistryModel<CompType, ThermoType>& chemistry
        ),
        (dict, chemistry)
    );


    // Constructors

        //- Construct from components
        chemistryReductionMethod
        (
            const IOdictionary& dict,
            TDACChemistryModel<CompType, ThermoType>& chemistry
        );


    // Selector

        static autoPtr<chemistryReductionMethod<CompType, ThermoType>> New
        (
            const IOdictionary& dict,
            TDACChemistryModel<CompType, ThermoType>& chemistry
        );


    //- Destructor
    virtual ~chemistryReductionMethod();


    // Member Functions

        //- Is mechanism reduction active?
        inline bool active() const;

        //- Is performance data logging enabled?
        inline bool log() const;

        //- Return the active species
        inline const List<bool>& activeSpecies() const;

        //- Return the number of active species
        inline label NsSimp();

        //- Return the initial number of species
        inline label nSpecie();

        //- Return the tolerance
        inline scalar tolerance() const;

        //- Reduce the mechanism
        virtual void reduceMechanism
        (
            const scalarField &c,
            const scalar T,
            const scalar p
        ) = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "chemistryModel/TDACChemistryModel/reduction/chemistryReductionMethod/chemistryReductionMethodI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "chemistryModel/TDACChemistryModel/reduction/chemistryReductionMethod/chemistryReductionMethod.C"
    #include "chemistryModel/TDACChemistryModel/reduction/chemistryReductionMethod/chemistryReductionMethodNew.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
