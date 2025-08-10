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
    (c) 2016 OpenCFD.Ltd

Class
    Foam::functionObjects::thermoCoupleProbes

Group
    grpUtilitiesFunctionObjects

Description
    Sample probe for temperature using a thermocouple.

    Uses the correlation:

    \f[
        Nu = 2.0 + \left(0.4 Re^{0.5} + 0.06 Re^{2/3}\right)*Pr^{0.4}
    \f]

Usage
    Example of function object specification:
    \verbatim
    probes
    {
        type            thermoCoupleProbes;
        libs            ("libutilityFunctionObjects.so");
        writeControl    timeStep;
        writeInterval   1;

        solver          rodas23;
        absTol          1e-12;
        relTol          1e-8;

        interpolationScheme cellPoint;

        // thermocouple properties
        rho             8908;
        Cp              440;
        d               1e-3;
        epsilon         0.85;

        radiationField  G;

        probeLocations
        (
            (0.5 0.5 0.5)
        );
        fields
        (
            T
        );
    }
    \endverbatim


SourceFiles
    thermoCoupleProbes.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_thermoCoupleProbes_H
#define functionObjects_thermoCoupleProbes_H

#include "probes/probes.H"
#include "ODESystem/ODESystem.H"
#include "ODESolvers/ODESolver/ODESolver.H"
#include "basicThermo/basicThermo.H"
#include "fluidThermo/fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                     Class thermoCoupleProbes Declaration
\*---------------------------------------------------------------------------*/

class thermoCoupleProbes
:
    public probes,
    public ODESystem
{
protected:

    // Protected data

        //- Thermocouple density
        scalar rho_;

        //- Thermocouple heat capacity
        scalar Cp_;

        //- Thermocouple diameter
        scalar d_;

        //- Thermocouple emissivity
        scalar epsilon_;

        //- Name of the velocity field
        word UName_;

        //- Name of the incident radiation field
        word radiationFieldName_;

         //- Fluid thermo reference
        const fluidThermo& thermo_;

        //- ODESolver
        autoPtr<ODESolver> odeSolver_;

        //- Cached thermocouple temperature
        scalarField Ttc_;


    // Protected Member Functions

        //- Sample and write a particular volume field
        template<class Type>
        void sampleAndWrite
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );


private:

        //- Disallow default bitwise copy construct
        thermoCoupleProbes(const thermoCoupleProbes&);

        //- Disallow default bitwise assignment
        void operator=(const thermoCoupleProbes&);


public:

    //- Runtime type information
    TypeName("thermoCoupleProbes");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        thermoCoupleProbes
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool loadFromFiles = false,
            const bool readFields = true
        );


    //- Destructor
    virtual ~thermoCoupleProbes();


    // ODE functions (overriding abstract functions in ODE.H)

    //- Number of ODE's to solve
    virtual label nEqns() const;

    virtual void derivatives
    (
        const scalar x,
        const scalarField& y,
        scalarField& dydx
    ) const;

    virtual void jacobian
    (
        const scalar t,
        const scalarField& y,
        scalarField& dfdt,
        scalarSquareMatrix& dfdy
    ) const;


    // Public Member Functions

        //- Sample and write
        virtual bool write();

        //- Execute, currently does nothing
        virtual bool execute();

        //- Read
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "thermoCoupleProbes/thermoCoupleProbesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
