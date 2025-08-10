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

#include "fields/fvPatchFields/derived/mappedFromFile/fileData.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * //

namespace Foam
{
    template<class Type>
    PtrList<Field<Type>> fileData<Type>::fields_(0);

    template<class Type>
    triSurface fileData<Type>::surface_ = triSurface();

    template<class Type>
    fileName fileData<Type>::filename_("");

    template<class Type>
    word fileData<Type>::fieldName_("");

    template<class Type>
    bool fileData<Type>::convertToMillimeters_(false);

    template<class Type>
    autoPtr<indexedOctree<treeDataPoint>> fileData<Type>::pointTreePtr_(nullptr);

    template<class Type>
    bool fileData<Type>::geometryLoaded_(false);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fileData<Type>::fileData()
:
    mapBack_()
{}


template<class Type>
Foam::fileData<Type>::fileData(const dictionary& dict)
:
    mapBack_(dict.lookupOrDefault<bool>("mapBack", false))
{
    if (filename_ == "")
    {
        filename_ = fileName(dict.lookup("file"));
        /*
        surface_ = triSurface(filename_);
        surface_.write("exported.stl",true,false);

        Info<< "Constructing Octree" << endl;
        const vectorField& pnts = surface_.Cf();

        //construct octree
        treeBoundBox overallBb(pnts);
        Random rndGen(123456);
        overallBb = overallBb.extend(rndGen, 1E-4);
        overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

        indexedOctree<treeDataPoint> pointTree
        (
            treeDataPoint(pnts),
            overallBb,  // overall search domain
            8,          // maxLevel
            10,         // leafsize
            3.0         // duplicity
        );

        pointTreePtr_ = pointTree.clone();
        */
    }
    if (fieldName_ == "") fieldName_ = word(dict.lookup("field"));
    convertToMillimeters_ = dict.lookupOrDefault<bool>
        (
            "convertToMillimeters",
            false
        );

    /*
    if (fields_.size() == 0)
    readFields();
    */
}

template<class Type>
Foam::fileData<Type>::fileData(const fileData& fd)
:
    mapBack_(fd.mapBack_)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::fileData<Type>::~fileData()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::fileData<Type>::readFields()
{
    if (filename_.substr(filename_.size()-3,filename_.size()) == "nas")
    {
        readFieldsNAS();
    }
    else if (filename_.substr(filename_.size()-3,filename_.size()) == "ntl")
    {
        readFieldsNTL();
    }
}

template<typename Type>
void Foam::fileData<Type>::readFieldsNTL()
{
    // Nonsense, but preserves original insane semantics without breaking the linker.
    if constexpr (!std::is_same_v<Type, scalar>) {
        return;
    }

    //const Field<point>& pnts = surface_.localPoints();
    const Field<point>& pnts = surface_.points();
    PtrList<scalarField> pntFields(2);
    pntFields.set(0, new Field<scalar> (pnts.size(), pTraits<scalar>::zero));
    pntFields.set(1, new Field<scalar> (pnts.size(), pTraits<scalar>::zero));

    Field<scalar> pntField= Field<scalar>
    (
        pnts.size(),
        pTraits<scalar>::zero
    );

    Info<< "Reading fields values from file " << filename_ << endl;
    IFstream is(filename_);

    label ctr(0);

    while (is.good())
    {
        size_t linei = 2;
        string line;
        is.getLine(line);

        // Reading Node Temperatures
        if (line.substr(0, 2) == "10")
        {
            if (ctr == pntField.size())
            {
                Info<< "exceeded list!!" << endl;
            }

            is.getLine(line);
            linei = 0;

            // Get front and back temperature on nodes
            pntFields[0][ctr] = parseScalar(readToken(line, 16, linei));
            pntFields[1][ctr] = parseScalar(readToken(line, 16, linei));

            ctr++;
        }
        else if (line.substr(0, 2) == " 1")
        {
            is.getLine(line);
            is.getLine(line);
        }
        // Reading Material Properties
        else if (line.substr(0, 2) == " 3")
        {
            is.getLine(line);
            is.getLine(line);
        }
        // Reading Element Properties
        else if (line.substr(0, 2) == " 4")
        {
            is.getLine(line);
        }
        // Reading Coordinate Frames
        else if (line.substr(0, 2) == " 5")
        {
            is.getLine(line);
            is.getLine(line);
            is.getLine(line);
            is.getLine(line);
        }
        // Reading Distributed Loads (Data Card number varies)
        else if (line.substr(0, 2) == " 6")
        {
            is.getLine(line);
        }
        // Reading Node Forces (Data Card number varies)
        else if (line.substr(0, 2) == " 7")
        {
            is.getLine(line);
        }
        // Reading Node Displacements
        else if (line.substr(0, 2) == " 8")
        {
            is.getLine(line);
            is.getLine(line);
        }
        // Reading Element Temperatures
        else if (line.substr(0, 2) == "11")
        {
            is.getLine(line);
        }
        // Reading MPC Data
        else if (line.substr(0, 2) == "14")
        {
            is.getLine(line);
            is.getLine(line);
            is.getLine(line);
        }
        // Reading Nodal Heat Source
        else if (line.substr(0, 2) == "15")
        {
            is.getLine(line);
        }
        // Reading Distributed Heat Source (Data Card number varies)
        else if (line.substr(0, 2) == "16")
        {
            is.getLine(line);
        }
        // Reading Convection Coefficients
        else if (line.substr(0, 2) == "17")
        {
            is.getLine(line);
            is.getLine(line);
        }
        // Reading Emissivity Values
        else if (line.substr(0, 2) == "18")
        {
            is.getLine(line);
            is.getLine(line);
        }
        if (line.substr(0, 2) == "25")
        {
            is.getLine(line);
        }
        // Reading/Skipping Summary Data
        else if (line.substr(0, 2) == "26")
        {
            is.getLine(line);
        }
        else if (line.substr(0, 2) == "98")
        {
            is.getLine(line);
        }
        else if (line.substr(0, 2) == "99")
        {
            break;
        }
        // Skip
        else
        {
            continue;
        }
    }

    Info<< "Interpolating fields from face points to centres" << endl;

    fields_.setSize(2);
    fields_.set
    (
        0,
        new Field<scalar> (surface_.Cf().size(), pTraits<scalar>::zero)
    );
    fields_.set
    (
        1,
        new Field<scalar> (surface_.Cf().size(), pTraits<scalar>::zero)
    );

    forAll(surface_, fI)
    {
        forAll(surface_[fI],pI)
        {
            label pntI = surface_[fI][pI];
            fields_[0][fI] += pntFields[0][pntI]/3.;
            fields_[1][fI] += pntFields[1][pntI]/3.;
        }
    }
    Info<< "Done" << endl;
}


template<typename Type>
void Foam::fileData<Type>::readFieldsNAS()
{
    // Nonsense, but preserves original insane semantics without breaking the linker.
    if constexpr (!std::is_same_v<Type, scalar>) {
        return;
    }

    fields_.setSize(2);
    forAll(fields_,fldI)
    {
        const Field<point>& pnts = surface_.points();
        Field<scalar> pntField= Field<scalar>
        (
            pnts.size(),
            pTraits<scalar>::zero
        );

        fileName fn=filename_;
        if (fldI==1) fn = fn.substr(0,fn.size()-4)+"_back.nas";

        Info<< "Reading fields values from file " << fn << endl;
        IFstream is(filename_);

        label ctr(0);

        while (is.good())
        {
            size_t linei = 0;
            string line;
            is.getLine(line);

            if (line.empty() || line[0] == '$')
            {
                // Skip empty or comment
                continue;
            }

            // Read first word
            word cmd(IStringStream(readToken(line, 8, linei))());

            if (cmd == fieldName_)
            {
                if (ctr == pntField.size())
                {
                    Info<< "exceeded list!!" << endl;
                }

                readToken(line, 8, linei);
                readToken(line, 8, linei);

                scalar fieldValue = parseScalar(readToken(line, 8, linei));

                pntField[ctr] = fieldValue;

                ctr++;
            }
        }

        Info<< "Interpolating fields from face points to centres" << endl;

        fields_.set(fldI, new Field<scalar> (surface_.Cf().size(), pTraits<scalar>::zero));

        forAll(surface_, fI)
        {
            forAll(surface_[fI],pI)
            {
                label pntI = surface_[fI][pI];
                fields_[fldI][fI] += pntField[pntI]/3.;
            }
        }

        Info<< "Done" << endl;
    }
}

template<class Type>
Foam::Field<Type>* Foam::fileData<Type>::mapData
(
    const pointField& tarPoints,
    const vectorField& tarNormals
)
{
    static_assert(std::is_same_v<Type, scalar>, "Function only supported for scalars");

    if (!geometryLoaded_)
    {
        Info<< "Reading thermal mesh" << endl;

        surface_ = triSurface(filename_);

        //if(debug)
        //surface_.write("exported.stl",true,false);
        //if(debug)
        //Info<< "Constructing Octree" << endl;

        const vectorField& pnts = surface_.Cf();

        //construct octree
        treeBoundBox overallBb(pnts);
        Random rndGen(123456);
        overallBb = overallBb.extend(rndGen, 1E-4);
        overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

        indexedOctree<treeDataPoint> pointTree
        (
            treeDataPoint(pnts),
            overallBb,  // overall search domain
            8,          // maxLevel
            10,         // leafsize
            3.0         // duplicity
        );

        pointTreePtr_ = pointTree.clone();

        geometryLoaded_ = true;
    }

    if (fields_.size() == 0)
    readFields();

    //if(debug)
    //Info<< "Mapping fields to CFD mesh" << endl;
    //Don't forget to add convertToMillimeter

    scalarField* mappedFieldPtr(nullptr);
    mappedFieldPtr = new scalarField(tarPoints.size(),0.);
    scalarField& mappedField = *mappedFieldPtr;

    /*
    const vectorField& pnts = surface_.Cf();

    //construct octree
    treeBoundBox overallBb(pnts);
    Random rndGen(123456);
    overallBb = overallBb.extend(rndGen, 1E-4);
    overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

    indexedOctree<treeDataPoint> pointTree
    (
        treeDataPoint(pnts),
        overallBb,  // overall search domain
        8,          // maxLevel
        10,         // leafsize
        3.0         // duplicity
    );
    */

    //loop over points, find nearest and assign mapped data
    //scalar spanSqr = Foam::sqr(pointTree.bb().mag());
    scalar spanSqr = Foam::sqr(pointTreePtr_->bb().mag());

    forAll(tarPoints, pI)
    {
        const point& p = tarPoints[pI];

        //pointIndexHit info = pointTree.findNearest
        pointIndexHit info = pointTreePtr_->findNearest
        (
            p,
            spanSqr
        );

        if (!info.hit())
        {
            //info = pointTree.findNearest
            info = pointTreePtr_->findNearest
            (
                p,
                GREAT*spanSqr
            );
        }

        label side = 0;

        if (mapBack_)
        if ((tarNormals[pI] & surface_.Sf()[info.index()]) > 0) side = 1;

        mappedField[pI] = fields_[side][info.index()];
    }

    //if(debug)
    //Info<< "Done" << endl;

    return mappedFieldPtr;
}

template<class Type>
Foam::scalar Foam::fileData<Type>::parseScalar(const string& s)
{
    size_t expSign = s.find_last_of("+-");

    if (expSign != string::npos && expSign > 0 && !isspace(s[expSign-1]))
    {
        scalar mantissa = readScalar(IStringStream(s.substr(0, expSign-1))());
        scalar exponent = readScalar(IStringStream(s.substr(expSign+1))());

        if (s[expSign] == '-')
        {
            exponent = -exponent;
        }
        return mantissa*pow(10, exponent);
    }
    else
    {
        return readScalar(IStringStream(s)());
    }
}

template<class Type>
std::string Foam::fileData<Type>::readToken
(
    const string& line,
    const size_t& width,
    size_t& index
)
{
    size_t indexStart, indexEnd;

    indexStart = index;

    indexEnd = line.find(',', indexStart);
    index = indexEnd + 1;

    if (indexEnd == std::string::npos)
    {
        indexEnd = indexStart + width;
        index = indexEnd;
    }

    return line.substr(indexStart, indexEnd - indexStart);
}

// ************************************************************************* //
