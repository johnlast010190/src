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
    (c) 2011 OpenFOAM Foundation
    (c) 2010-2022 Esi Ltd.

Class
    Foam::thrusterActuationDiskSource

Description

    Example usage:
    \verbatim
    thrusterActuationDiskSourceCoeffs
    {
        fieldNames ( U );

        // Polynomial coefficients for pressure jump as function of
        // normal velocity [m/s]
        // (if fanMode: averageVelocity or localVelocity)
        // or flow rate (ico [m3/s], cmp [kg/s]) (if fanMode: flowRate)
        f List<scalar> 7(40.0753 -7.61342 4.80957 -20.6546 14.5827 -3.77995 0.334983);

        // Fan mode (averageVelocity, localVelocity, flowRate)
        fanMode averageVelocity;

        // Clip the value to these bounds before applying to polynomial.
        // If fanMode = averageVelocity or localVelocity, the entry is
        // the velocity value to be clipped.
        // If fanMode = flowRate, the entry is the flowRate to be clipped.
        lowerBound 0;
        upperBound 3.92;

        // Load profile as a function of relative radius
        loadProfile ( 1  1 );

        // Blade angle in degrees as a function of relative radius
        alphaProfile ( 0 0 );

        // Reference frame name
        referenceFrame thrusterFrame;
    }
    \endverbatim


SourceFiles
    thrusterActuationDiskSource.C
    thrusterActuationDiskSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef thrusterActuationDiskSource_H
#define thrusterActuationDiskSource_H

#include "containers/Lists/FixedList/FixedList.H"
#include "cellSetOption/cellSetOption.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "primitives/functions/Polynomial/polynomialFunction.H"
#include "referenceFrames/coordinateFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                  Class thrusterActuationDiskSource Declaration
\*---------------------------------------------------------------------------*/

class thrusterActuationDiskSource
:
    public cellSetOption
{

public:

    // Data types

        //- Enumeration listing the possible fan input model
        enum fanModelType
        {
            fsAverageVelocity,
            fsLocalVelocity,
            fsFlowRate
        };

private:

    // Private data

        //- Fan source type names
        static const NamedEnum<fanModelType, 3> fanModelTypeNames_;

        //- Fan source type
        fanModelType fanModel_;

        //- Pressure gradient
        scalarField load_;

        //- velocity jump
        vectorField jumpV_;

        //- Fan pressure rise polynomial
        polynomialFunction f_;

        //- Lower bound value to be clipped (velocity or flowRate)
        const scalar lowerBound_;

        //- Upper bound value to be clipped (velocity or flowRate)
        const scalar upperBound_;

        //- Reference frame
        coordinateFrame* coorFramePtr_;

        //- Load profile as a function of relative radius
        autoPtr<Function1<scalar>> loadProfile_;

        //- Blade angle profile as a function of relative radius, in degrees
        autoPtr<Function1<scalar>> alphaProfile_;

        //- Density for incompressible case
        scalar rhoRef_;
        bool rhoRefFound_;

       //- Disk front area
       scalar diskArea_;

       //- Is compressible?
       bool compressible_;

       //- zone boundary faces for flow rate calc
       labelList zoneBoundaryFaces_;


    // Private Member Functions

        //- Add resistance to the UEqn
        template<class RhoFieldType>
        void addThrusterActuationDiskAxialInertialResistance
        (
            vectorField& Usource,
            const labelList& cells,
            const scalarField& V,
            const RhoFieldType& rho,
            const vectorField& U
        ) const;

        //- Disallow default bitwise copy construct
        thrusterActuationDiskSource(const thrusterActuationDiskSource&);

        //- Disallow default bitwise assignment
        void operator=(const thrusterActuationDiskSource&);

        //- cellZone boundary faces
        labelList zoneBoundaryFaces() const;

        //- return sum of mag of face flux: ico [m3/s], cmp [kg/s]
        scalar zoneFlux(const fvMesh& mesh, const labelList& faces) const;

public:

    //- Runtime type information
    TypeName("thrusterActuationDiskSource");


    // Constructors

        //- Construct from components
        thrusterActuationDiskSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~thrusterActuationDiskSource()
    {}


        // Public Functions

            //- Source term to fvMatrix<vector>
            virtual void addSup(fvMatrix<vector>& eqn, const label fieldI);

            //- Source term to fvMatrix<vector>
            virtual void addSup(const volScalarField& rho, fvMatrix<vector>& eqn, const label fieldI);

        // I-O

            //- Write data
            virtual void writeData(Ostream&) const;

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "thrusterActuationDiskSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
