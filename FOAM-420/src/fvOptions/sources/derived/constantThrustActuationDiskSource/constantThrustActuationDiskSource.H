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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2013 Kevin Maki
    (c) 2013-2019 Esi Ltd.

Class
    Foam::fv::constantThrustActuationDiskSource

Description
    Constant thrust actuation disk source

SourceFiles
    constantThrustActuationDiskSource.C
    constantThrustActuationDiskSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef constantThrustActuationDiskSource_H
#define constantThrustActuationDiskSource_H

#include "cellSetOption/cellSetOption.H"
#include "coordinate/systems/coordinateSystem.H"
#include "coordinate/systems/cartesianCS.H"
#include "referenceFrames/coordinateFrame.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class constantThrustActuationDiskSource Declaration
\*---------------------------------------------------------------------------*/

class constantThrustActuationDiskSource
:
    public cellSetOption
{
public:

    // Public data

        //- Enumeration for volume types
        enum volumeModeType
        {
            vmAbsolute,
            vmSpecific
        };

        //- Word list of volume mode type names
        static const NamedEnum<volumeModeType, 2> volumeModeTypeNames_;


protected:

    // Protected data

        //- Volume mode
        volumeModeType volumeMode_;

        //- Volume normalisation
        scalar VDash_;

        //- Disk area normal
        vector diskDir_;

        //- Thrust
        autoPtr<Function1<scalar>> TF1_;
        scalar T_;

        //- Local coordinate System
        autoPtr<coordinateSystem> coorSysPtr_;

        //- Reference frame ptr
        coordinateFrame* coorFramePtr_;

        //- Access to local coordinate system
        inline const coordinateSystem& csys() const;

private:

    // Private Member Functions

        //- Check data
        void checkData() const;

        //- Create coordinate system
        void createCS();

        //- Add resistance to the UEqn
        template<class RhoFieldType>
        void addActuationDiskAxialInertialResistance
        (
            vectorField& Usource,
            const labelList& cells,
            const scalarField& V,
            const RhoFieldType& rho,
            const vectorField& U
        ) const;

        //- Disallow default bitwise copy construct
        constantThrustActuationDiskSource
        (
            const constantThrustActuationDiskSource&
        );

        //- Disallow default bitwise assignment
        void operator=(const constantThrustActuationDiskSource&);


public:

    //- Runtime type information
    TypeName("constantThrustActuationDiskSource");


    // Constructors

        //- Construct from components
        constantThrustActuationDiskSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~constantThrustActuationDiskSource()
    {}


    // Member Functions

        // Access

            //- Return Thrust
            scalar T() const
            {
                return T_;
            }

            //- Normal disk direction
            const vector& diskDir() const
            {
                return diskDir_;
            }

        // Public Functions

            //- Source term to fvMatrix<vector>
            virtual void addSup(fvMatrix<vector>& eqn, const label fieldI);
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            );


        // I-O

            //- Write data
            virtual void writeData(Ostream&) const;

            //- Read dictionary
            virtual bool read(const dictionary& dict);
            bool readThrustData(const dictionary& dict);

            //- Update Thrust
            bool updateThrust(const scalar& thrust);

            //- Update Disk direction
            bool updateDiskDirection(const vector& diskDirection);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "constantThrustActuationDiskSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
