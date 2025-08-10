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
    (c) 2009-2009 OpenCFD Ltd.
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/mappedProfile/pointProfileData.H"
#include "algorithms/indexedOctree/indexedOctree.H"
#include "indexedOctree/treeDataPoint.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointProfileData::pointProfileData()
:
    regionName_(""),
    N_(0),
    fieldNames_(0),
    fields_(),
    points_(0)
{}


Foam::pointProfileData::pointProfileData(ISstream& is)
:
    regionName_(""),
    N_(0),
    fieldNames_(0),
    fields_(),
    points_(0)
{
    read(is);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointProfileData::~pointProfileData()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalarField* Foam::pointProfileData::mapData
(
    const pointField& targetPoints,
    word fieldName
)
{
    scalarField* mf(nullptr);

    label fieldIndex = -1;

    forAll(fieldNames_, fnI)
    {
        if (fieldName == fieldNames_[fnI])
        {
            fieldIndex = fnI;
            break;
        }
    }

    if (fieldIndex != -1)
    {
        mf = new scalarField(targetPoints.size(), -1);

        scalarField& mfR = *mf;

        //construct octree
        treeBoundBox overallBb(points_);
        Random rndGen(123456);
        overallBb = overallBb.extend(rndGen, 1E-4);
        overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

        indexedOctree<treeDataPoint> pointTree
        (
            treeDataPoint(points_),
            overallBb,  // overall search domain
            8,                              // maxLevel
            10,                             // leafsize
            3.0                             // duplicity
        );

        //loop over points, find nearest and assign mapped data
        scalar spanSqr = Foam::sqr(pointTree.bb().mag());

        forAll(targetPoints, pI)
        {
            const point& p = targetPoints[pI];

            pointIndexHit info = pointTree.findNearest
            (
                p,
                spanSqr
            );

            if (!info.hit())
            {
                info = pointTree.findNearest
                (
                    p,
                    GREAT*spanSqr
                );
            }

            mfR[pI] = fields_[fieldIndex][info.index()];
        }


    }
    else
    {
        WarningInFunction
            << "Specified target field " << fieldName << " not found in "
            << "profile data." << nl << "Available fields are: " << fieldNames_
            << endl;
    }

    return mf;
}

void Foam::pointProfileData::read(Foam::ISstream& is)
{

    //read leading bracket
    token lbToken(is);
    checkIO(lbToken, token::BEGIN_LIST);

    //read boundary region identifier
    regionName_ = word(is);

    //read "points" type id
    word pointID(is);

    //read size
    N_ = readLabel(is);

    //closing bracket
    token cbToken(is);
    checkIO(cbToken, token::END_LIST);

    // now loop until all data for region has been read
    bool endRegion = false;
    label nFields = 0;
    while (!endRegion)
    {
        token currToken(is);

        if (currToken == token::BEGIN_LIST)
        {
            nFields++;
            fieldNames_.setSize(nFields);
            fieldNames_[nFields-1] = word(is);
            fields_.setSize(nFields);
            fields_.set(nFields-1, new scalarField(N_, 0));
            scalarField& csf = fields_[nFields-1];

            forAll(csf, pI)
            {
                csf[pI] = readScalar(is);
            }

            //read last bracket
            token endToken(is);
        }
        else if (currToken == token::END_LIST)
        {
            endRegion = true;
        }
        else
        {
            FatalIOErrorIn("pointProfileData::read(Istream&)", is)
                << "Unexpected token " << currToken
                << " in input stream. Expected '(' or ')'"
                << exit(FatalIOError);

        }
    }

    // build points from xyz
    points_.setSize(N_);
    label foundXYZ = 0;

    forAll(fieldNames_, fI)
    {
        if (fieldNames_[fI] == word("x"))
        {
            foundXYZ++;

            forAll(points_, pI)
            {
                points_[pI].x() = fields_[fI][pI];
            }
        }
        if (fieldNames_[fI] == word("y"))
        {
            foundXYZ++;

            forAll(points_, pI)
            {
                points_[pI].y() = fields_[fI][pI];
            }
        }
        if (fieldNames_[fI] == word("z"))
        {
            foundXYZ++;

            forAll(points_, pI)
            {
                points_[pI].z() = fields_[fI][pI];
            }
        }
    }

    if (foundXYZ != 3)
    {
        FatalIOErrorIn("pointProfileData::read(Istream&)", is)
            << "Could not find one or more input point coordinates x, y, z."
            << nl << tab << "Available point fields are: " << fieldNames_
            << exit(FatalIOError);

    }

}

void Foam::pointProfileData::checkIO
(
    const token& t,
    const token::punctuationToken& pt
)
{
    if (t != pt)
    {
        FatalErrorInFunction
            << "Unexpected token " << t
            << " in input stream. Expected " << pt
            << exit(FatalIOError);
    }
}

void Foam::pointProfileData::nextValidToken(ISstream& is)
{
    token cToken(is);

    while
    (
        cToken == token::NULL_TOKEN
        || cToken == token::SPACE
        || cToken == token::TAB
        || cToken == token::NL
    )
    {
        //Info<< cToken << endl;
        cToken = token(is);

    }

    //Info<< cToken << endl;
    is.putBack(cToken);
}

// ************************************************************************* //
