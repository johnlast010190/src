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
    (c) 2010-2020 Esi Ltd.
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::phase

Description
    Single incompressible phase derived from the phase-fraction.
    Used as part of the multiPhaseMixture for interface-capturing multi-phase
    simulations.

SourceFiles
    phase.C

\*---------------------------------------------------------------------------*/

#ifndef phase_H
#define phase_H

#include "fields/volFields/volFields.H"
#include "db/dictionary/dictionaryEntry/dictionaryEntry.H"
#include "incompressible/physicalProperties/physicalProperties.H"
#include "relativeVelocityModel/relativeVelocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phase Declaration
\*---------------------------------------------------------------------------*/

class phase
:
    public volScalarField
{
    // Private data

        word name_;
        dictionary phaseDict_;
        scalar alphaMax_;
        bool useTotalSolids_;
        autoPtr<dimensionedScalar> rho_;
        autoPtr<dimensionedScalar> nu_;
        autoPtr<relativeVelocityModel> UdmModelPtr_;


public:

    // Constructors

        //- Construct from components
        phase
        (
            const word& name,
            const dictionary& phaseDict,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const scalar& alphaMax
        );

        //- Return clone
        autoPtr<phase> clone() const;

        //- Return a pointer to a new phase created on freestore
        //  from Istream
        class iNew
        {
            const volVectorField& U_;
            const surfaceScalarField& phi_;
            const scalar& alphaMax_;

        public:

            iNew
            (
                const volVectorField& U,
                const surfaceScalarField& phi,
                const scalar& alphaMax
            )
            :
                U_(U),
                phi_(phi),
                alphaMax_(alphaMax)
            {}

            autoPtr<phase> operator()(Istream& is) const
            {
                dictionaryEntry ent(dictionary::null, is);
                return autoPtr<phase>(new phase(ent.keyword(), ent, U_, phi_, alphaMax_));
            }
        };


    // Member Functions

        const word& name() const
        {
            return name_;
        }

        const word& keyword() const
        {
            return name();
        }

        const dictionary& dict() const
        {
            return phaseDict_;
        }

        //- Return const-access to phase density
        const dimensionedScalar rho() const
        {
            return rho_();
        }

        //- Return const-access to phase viscosity
        const dimensionedScalar nu() const
        {
            return nu_();
        }

        //- Return const-access to drift velocity
        const volVectorField& Udm() const
        {
            return UdmModelPtr_->Udm();
        }

        //- Return relative velocity
        tmp<volVectorField> Urel() const
        {
            return -Udm();
        }

        //- Return alpha23 from drift flux model
        scalar alpha23() const
        {
            return UdmModelPtr_->alpha23();
        }

        //- Correct the phase properties
        void correct()
        {}

        //- Correct the drift velocity
        void correctUdm();

        //- use total solids y/n
        bool useTotalSolids() const
        {
            return useTotalSolids_;
        }

        //- Inherit read from volScalarField
        using volScalarField::read;

        //- Read base transportProperties dictionary
        bool read(const dictionary& phaseDict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
