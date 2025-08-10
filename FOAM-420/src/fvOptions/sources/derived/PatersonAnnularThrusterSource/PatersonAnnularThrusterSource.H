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
    (c) 2016-2021 Esi Ltd

Class
    Foam::fv::PatersonAnnularThrusterSource

Description
    Actuation disk source including radial thrust

SourceFiles
    PatersonAnnularThrusterSource.C
    PatersonAnnularThrusterSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef PatersonAnnularThrusterSource_H
#define PatersonAnnularThrusterSource_H

#include "primitives/functions/Polynomial/Polynomial.H"
#include "forces/forces.H"
#include "cellSetOption/cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                  Class PatersonAnnularThrusterSource Declaration
\*---------------------------------------------------------------------------*/

class PatersonAnnularThrusterSource
:
    public cellSetOption
{
    // Private data


        //- Ship data

            //- drag force calculation
            functionObjects::forces shipForce_;

            //- Uship - ship reference velocity
            vector Uship_;

            //- reference density
            scalar rhoRef_;

        //Propeller data

            //- Disk Radius
            scalar R_;

            //- hub Radius
            scalar rh_;

            //- propeller root position
            vector propPosition_;

            //- disk direction
            vector diskDir_;

            //- propeller thickness from root
            scalar tProp_;

            //- torque direction around diskDir, can be positive or negative
            scalar qPropDir_;

        //- Propeller thrust coefficients

            //- diskFraction - undocumented
            scalar diskFraction_;

            //- J - undocumented
            scalar J_;

            //- "a" coefficients - undocumented
            Polynomial<6> a_;

            //- "b" coefficients - undocumented
            Polynomial<6> b_;

            //- normalise forces to produce equivalence with integral forces
            Switch normaliseForces_;

            //- compute and output wake fraction
            Switch wakeFraction_;


    // Private Member Functions

        //- calculate thrust
        vector calcDrag();

        //- calculate torque
        vector calcTorque(const vector&);

        //- analytical annular volume
        scalar annularVolume() const;

        //- Check data
        void checkData() const;

        //- Add resistance to the UEqn
        template<class RhoFieldType>
        void addBodyForce
        (
            vectorField& Usource,
            const labelList& cells,
            const scalarField& V,
            const RhoFieldType& rho,
            const vectorField& U
        );

        //- Disallow default bitwise copy construct
        PatersonAnnularThrusterSource(const PatersonAnnularThrusterSource&);

        //- Disallow default bitwise assignment
        void operator=(const PatersonAnnularThrusterSource&);

protected:

        //- Set the cell set for annular disk
        virtual void setCellSet();



public:

    //- Runtime type information
    TypeName("PatersonAnnularThrusterSource");


    // Constructors

        //- Construct from components
        PatersonAnnularThrusterSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~PatersonAnnularThrusterSource()
    {}


    // Member Functions

        //- Source term to momentum equation
        virtual void addSup
        (
            fvMatrix<vector>& eqn,
            const label fieldI
        );

        //- Source term to compressible momentum equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldI
        );

    // I-O

        //- Read dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PatersonAnnularThrusterSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
