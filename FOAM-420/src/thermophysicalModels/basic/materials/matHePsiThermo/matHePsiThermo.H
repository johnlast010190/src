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
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2021 Esi Ltd.

Class
    Foam::matHePsiThermo

Description
    Energy for a mixture based on density

SourceFiles
    matHePsiThermo.C

\*---------------------------------------------------------------------------*/

#ifndef matHePsiThermo_H
#define matHePsiThermo_H

#include "materials/matHeThermo/matHeThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class matHePsiThermo Declaration
\*---------------------------------------------------------------------------*/
template<class BasicPsiThermo, class MixtureType>
class matHePsiThermo
:
    public matHeThermo<BasicPsiThermo, MixtureType>
{
    // Private Member Functions

        //- Calculate the thermo variables
        void calculate
        (
            const volScalarField& p,
            volScalarField& T,
            volScalarField& he,
            volScalarField& psi,
            volScalarField& mu,
            volScalarField& alpha,
            const bool doOldTimes
        );

        //- Construct as copy (not implemented)
        matHePsiThermo(const matHePsiThermo&);

        //- Object registry reference
        const objectRegistry& obr_;


public:

        TypeName("matHePsiThermo");
    // Constructors

        //- Construct from obr and phase name
        matHePsiThermo
        (
            const objectRegistry& obr,
            const word& phaseName = word::null
        );


    //- Destructor
    virtual ~matHePsiThermo();


    // Member functions

        //- Update properties
        virtual void correct();

        //- Density [kg/m^3] - uses current value of pressure
        //  This should not be called for psiThermo
        virtual tmp<volScalarField> rho() const
        {
            // This expression is only valid if rho is linear with respect to p
            // If the solver calls rho(), is is making that assumption, which is
            // dangerous
            FatalErrorInFunction
                << "This solver does not support the 'psiFluid' material type; "
                << "please use 'fluid' type instead."
                << exit(FatalError);
            return BasicPsiThermo::rho();
        }

        //- Density for patch [kg/m^3]
        //  This should not be called for psiThermo
        virtual tmp<scalarField> rho(const label patchi) const
        {
            // This expression is only valid if rho is linear with respect to p
            // If the solver calls rho(), is is making that assumption, which is
            // dangerous
            FatalErrorInFunction
                << "This solver does not support the 'psiFluid' material type; "
                << "please use 'fluid' type instead."
                << exit(FatalError);
            return BasicPsiThermo::rho(patchi);
        }


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "matHePsiThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
