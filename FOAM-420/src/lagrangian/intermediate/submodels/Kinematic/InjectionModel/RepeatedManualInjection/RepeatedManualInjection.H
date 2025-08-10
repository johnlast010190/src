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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

Class
    Foam::RepeatedManualInjection

Description
    Manual injection
    - User specifies
      - Total mass to inject
      - Parcel positions in file \c positionsFile
      - Initial parcel velocity
    - Parcel diameters obtained by distribution model
    - All parcels introduced at SOI

SourceFiles
    RepeatedManualInjection.C

\*---------------------------------------------------------------------------*/

#ifndef RepeatedManualInjection_H
#define RepeatedManualInjection_H

#include "submodels/Kinematic/InjectionModel/InjectionModel/InjectionModel.H"
#include "primitives/bools/Switch/Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class distributionModel;

/*---------------------------------------------------------------------------*\
                       Class RepeatedManualInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class RepeatedManualInjection
:
    public InjectionModel<CloudType>
{
    // Private data
        //- Injection duration [s]
        scalar duration_;

        //- Name of file containing positions data
        const word positionsFile_;

        //- Field of parcel positions
        vectorIOField positions_;

        //- Field of parcel diameters
        scalarList diameters_;

        //- List of cell labels corresponding to injector positions
        labelList injectorCells_;

        //- List of tetFace labels corresponding to injector positions
        labelList injectorTetFaces_;

        //- List of tetPt labels corresponding to injector positions
        labelList injectorTetPts_;

        //- Initial parcel velocity
        const vector U0_;

        //- injection interval
        scalar injectionInterval_;

        //- injection interval
        label injectionTimeIndex_;

        scalar deltaT_;

        bool inject_;

        //- Parcel size distribution model
        const autoPtr<distributionModel> sizeDistribution_;

        //- Flag to suppress errors if particle injection site is out-of-bounds
        Switch ignoreOutOfBounds_;


public:

    //- Runtime type information
    TypeName("repeatedManualInjection");


    // Constructors

        //- Construct from dictionary
        RepeatedManualInjection
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        RepeatedManualInjection(const RepeatedManualInjection<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType>> clone() const
        {
            return autoPtr<InjectionModel<CloudType>>
            (
                new RepeatedManualInjection<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~RepeatedManualInjection();


    // Member Functions

        //- Set injector locations when mesh is updated
        virtual void updateMesh();

        //- Return the end-of-injection time
        scalar timeEnd() const;

        //- Number of parcels to introduce relative to SOI
        virtual label parcelsToInject(const scalar time0, const scalar time1);

        //- Volume of parcels to introduce relative to SOI
        virtual scalar volumeToInject(const scalar time0, const scalar time1);


        // Injection geometry

            //- Set the injection position and owner cell, tetFace and tetPt
            virtual void setPositionAndCell
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                vector& position,
                label& cellOwner,
                label& tetFaceI,
                label& tetPtI
            );

            //- Set the parcel properties
            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                typename CloudType::parcelType& parcel
            );

            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const;

            //- Return flag to identify whether or not injection of parcelI is
            //  permitted
            virtual bool validInjection(const label parcelI);

        // recursive modulus
        double modulus(double,double);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "submodels/Kinematic/InjectionModel/RepeatedManualInjection/RepeatedManualInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
