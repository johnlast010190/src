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
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::XiEqModel

Description
    Base-class for all XiEq models used by the b-XiEq combustion model.
    The available models are :
        \link basicXiSubXiEq.H \endlink
        \link Gulder.H \endlink
        \link instabilityXiEq.H \endlink
        \link SCOPEBlendXiEq.H \endlink
        \link SCOPEXiEq.H \endlink

SourceFiles
    XiEqModel.C

\*---------------------------------------------------------------------------*/

#ifndef XiEqModel_H
#define XiEqModel_H

#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "psiuReactionThermo/psiuReactionThermo.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class XiEqModel Declaration
\*---------------------------------------------------------------------------*/

class XiEqModel
{

protected:

    // Protected data

        //- Dictionary
        dictionary XiEqModelCoeffs_;

        //- Thermo
        const psiuReactionThermo& thermo_;

        //- Turbulence
        const compressible::RASModel& turbulence_;

        //- Laminar burning velocity
        const volScalarField& Su_;


private:

    // Private Member Functions

        //- Disallow copy construct
        XiEqModel(const XiEqModel&);

        //- Disallow default bitwise assignment
        void operator=(const XiEqModel&);


public:

    //- Runtime type information
    TypeName("XiEqModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            XiEqModel,
            dictionary,
            (
                const dictionary& XiEqProperties,
                const psiuReactionThermo& thermo,
                const compressible::RASModel& turbulence,
                const volScalarField& Su
            ),
            (
                XiEqProperties,
                thermo,
                turbulence,
                Su
            )
        );


    // Selectors

        //- Return a reference to the selected XiEq model
        static autoPtr<XiEqModel> New
        (
            const dictionary& XiEqProperties,
            const psiuReactionThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su
        );


    // Constructors

        //- Construct from components
        XiEqModel
        (
            const dictionary& XiEqProperties,
            const psiuReactionThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su
        );


    //- Destructor
    virtual ~XiEqModel();


    // Member Functions

        //- Return the flame-wrinking XiEq
        virtual tmp<volScalarField> XiEq() const
        {
            return turbulence_.muEff();
        }

        //- Return the sub-grid Schelkin effect
        tmp<volScalarField> calculateSchelkinEffect(const scalar) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& XiEqProperties) = 0;

        //- Write fields
        void writeFields() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
