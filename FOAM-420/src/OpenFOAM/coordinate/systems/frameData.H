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
    (c) 2022-2023 Esi Ltd.

Class
    Foam::cartesian

Description
    A Cartesian coordinate system

    \heading Dictionary entries
    \table
        Property    | Description                           | Required | Default
        type        | Type name: cartesian                  | yes   |
    \endtable

SourceFiles
    cartesianCS.C

\*---------------------------------------------------------------------------*/

#ifndef frameData_H
#define frameData_H

#include "primitives/enums/NamedEnum.H"
#include "primitives/Vector/vector/vector.H"
#include "memory/autoPtr/autoPtr.H"
#include "memory/tmp/tmp.H"
#include "db/objectRegistry/objectRegistry.H"
#include "db/IOobjects/IOdictionary/localIOdictionary.H"
#include "db/IOstreams/Fstreams/IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class coordinateSystem;

class coordinateFrame;

/*---------------------------------------------------------------------------*\
                   Class cartesian Declaration
\*---------------------------------------------------------------------------*/

class frameData
{
public:

    // Public data

        //- Enumeration for relation types
        enum relationType
        {
            rtGlobal,
            rtParent,
            rtFrame
        };

        //- Word list of relation type names
        static const NamedEnum<relationType, 3> relationTypeNames_;

protected:

    // Protected Data

        //- Frame coordinate system is defined in
        relationType definedInFrame_;

        //- Name of related frame (rtFrame only)
        word relatedFrame_;

        //- Coordinate system definition is local y/n
        bool isLocal_;

        //- Coordinate system definition is local y/n
        bool isLocalOnRead_;


public:


    // Constructors

        //- Construct null. This is an identity coordinateSystem.
        frameData();

        //- Copy construct
        frameData(const frameData& fd);

        //- Construct from dictionary without a name
        frameData(const dictionary& dict);


    //- Destructor
    virtual ~frameData(){};


    // Member Functions

        //- Transform coordinate system into parent coordinates
        // Note: not necessarily into direct parent
        // needs to be called successively starting with global/frame
        // and recursively transform until direct parent
        static void transformToLocal
        (
            coordinateSystem& csys,
            const coordinateSystem& csysParent
        );

        //- Transform coordiante system into global coordinates
        static void transformToGlobal
        (
            coordinateSystem& csys,
            const coordinateSystem& csysParent
        );

        //- Transform coordiante system into global coordinates
        static autoPtr<coordinateSystem> frameToGlobal
        (
            const dictionary& dict,
            const objectRegistry& obr
        );

        //- Set isLocal flag once transformation is done
        void setLocal(const bool& local)
        {
            isLocal_ = local;
        }

        //- Is frame local?
        bool isLocal() const
        {
            return isLocal_;
        }

        //- Is frame local?
        bool isLocalOnRead() const
        {
            return isLocalOnRead_;
        }


    // IO

        //- Read dictionary
        bool read(const dictionary& dict);

        //- Write flags
        void write(Ostream& os) const;


    // Member Operators

        //- Copy assignment
        void operator=(const frameData& fd);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
