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
    (c) 2018-2019 Esi Ltd.
    (c) 2017 D.Segersson SMHI
    (c) 2011-2017 OpenFOAM Foundation

Class
    Foam::fv::vegetationSource

Group
    grpFvOptionsSources

Original Author
    -momentum and turbulence source terms
    David Segersson,
    Swedish Meteorological and Hydrological Institute,
    david.segersson@smhi.se

Modifications (by ESI Ltd)
    - Class refactoring. Inherits from cellSetOption to be able to specify the
    source term on cellZones.
    - Implementation of the leaf heat balance equation for temperature and
    humidity source terms.


Description
    Base class for tree canopy source term. Includes implementation for
    momentum, temperature and humidity source terms.
    The leaf heat balance equation is solved for each cell of the vegetation
    cell zone. The model requires radiation, solar radiation and humidity
    equations active.

    type            vegetationSource;
    active          true;
    hookOp          solve; // add this only to solve leaf heat balance
                           // when T and w are active
    vegetationSourceCoeffs
    {
        selectionMode   cellZone;
        cellZone        tree;
        fields          ( U k epsilon );
        origin          ( 0 0 1 ); // at the tree base
        direction       ( 0 0 1 ); // along the tree height
        Cd              0.2;   // drag coefficient of trees
        LAI             1; // Leaf Area Index (used if LADmax is not specified)
        LADmax          4; // Maximum Leaf Area Density [m ^-1]
        LADProfile ( 0.05 0.1 0.15 0.35 1.1 0.9 0.5 0.2 0.15 0.05 0.01 );
        as              0.5; // Absorbed fraction of the solar radiation
        at              0.04;// Absorbed fraction of the thermal radiation
        leafSize        0.1; // leaf size [m]
        uniformLAD      false;
        }

    // Cd drag coefficient
    // LADProfile is the vertical profile of Leaf Area Density [m2/m3]
    // first value is closest to ground.
    // Each value represents an equal share of the tree height
    // e.g. for a 4 m tree and 4 values, each value will represent 1 m.
    // Values are scaled so that the highest will correspond to LADMax.

    // The user needs to set:
    // 1) LADProfile and LADmax, or
    // 2) LADProfile and LAI (LADmax calculated accordingly) or
    // 3) uniformLAD true and LADmax

    // The origin is a point at the base of the tree, with the
    // tree height measured along the direction


    References:
    D.Sergesson
    A tutorial to urban wind flow using OpenFOAM for course CFD. Proceedings
    of CFD with OpenSource Software, 2017.
    J.Carmeliet et all.
    Transpirative cooling potential of vegetation in urban environment using
    coupled CFD and leaf energy balance. Building Simulation 2017.


SourceFiles
    VegetationSource.C

\*---------------------------------------------------------------------------*/

#ifndef vegetationSource_H
#define vegetationSource_H

#include "cfdTools/general/fvOptions/fvOption.H"
#include "cellSetOption/cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class vegetationSource Declaration
\*---------------------------------------------------------------------------*/

class vegetationSource
:
    public cellSetOption
{

protected:

    // Protected data

        //- Cell Origin
        vector origin_;

        //- Axis direction
        vector zoneDir_;

        // U field name
        word UName_;

        // T field name
        word TName_;

        // p field name
        word PName_;

        // energy field name
        word eName_;

        // uniform LAD distribution
        bool uniformLAD_;

        // Tree height
        scalar height_;

        // Drag coefficient
        scalar Cd_;

        // Characteristic leaf size [m]
        scalar lsize_;

        // Leaf area index
        scalar LAI_;

        // Leaf area density max
        scalar LADmax_;

        // Canopy absorption coeff (solar)
        scalar as_;

        // Canopy absorption coeff (G)
        scalar at_;

        //- Molecular weight of vapour phase
        scalar Mvap_;

        //- Molecular weight of air
        scalar Mair_;

        //  Latent heat of vaporisation [J/kg]
        scalar lambdaVap_;

        // Leaf area density profile (normalised)
        scalarList LADProfile_;

        // canopy term Cd*LAD
        autoPtr<volScalarField> canopy_;

        // Leaf Area Density field [m2/m3]
        scalarField LAD;

        // sensible heat flux [W/m2]
        scalarField Qs;

        // latent heat flux [W/m2]
        scalarField Ql;

        // Solve leaf heat balance equation to calculate fluxes
        void leafHeatBalance();

        //- Return density
        virtual tmp<volScalarField> getRho() const;

        //- Return Cp
        virtual tmp<volScalarField> getCp() const;

private:

    // Private Member Functions

        // returns LAD for a given height
        scalar setLAD(scalar z);

        // calculates LADmax from LAI and LADProfile_
        void LADmaxFromLAI();

        // function used to integrate LAD
        scalar integrateLAD();

        // calculates the canopy term as Cd*LAD
        void calculateCanopy();

        // calculate distance w.r.t. cell centers
        void calculatePointDistance(scalarField &distance);

        // Init function
        void init();

        // member functions
        const scalarList& LADProfile()
        {
            return LADProfile_;
        }

        //- Disallow default bitwise copy construct
        vegetationSource(const vegetationSource&);

        //- Disallow default bitwise assignment
        void operator=(const vegetationSource&);


public:

    //- Runtime type information
    TypeName("vegetationSource");


    // Constructors

        //- Construct from components
        vegetationSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


        //- Destructor
        virtual ~vegetationSource()
        {}


    // Member Functions

        // Access

         // Add explicit and implicit contributions

            //- Source term to momentum equation
        virtual void addSup
        (
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Source term to compressible momentum equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Source term to temperature and humidity equations - ico
        virtual void addSup
        (
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Source term to energy and humidity equations - cmp
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        // Solve function - used for the solution of the leaf heat balance
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
