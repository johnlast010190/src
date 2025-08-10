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

\*---------------------------------------------------------------------------*/

#include "matrices/Matrix/Matrix.H"
#include "db/IOstreams/IOstreams/Istream.H"
#include "db/IOstreams/IOstreams/Ostream.H"
#include "db/IOstreams/token/token.H"
#include "primitives/contiguous/contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(Istream& is)
:
    mRows_(0),
    nCols_(0),
    v_(nullptr)
{
    operator>>(is, *this);
}


template<class Form, class Type>
Foam::Istream& Foam::operator>>(Istream& is, Matrix<Form, Type>& M)
{
    // Anull matrix
    M.clear();

    is.fatalCheck(FUNCTION_NAME);

    token firstToken(is);

    is.fatalCheck
    (
        "operator>>(Istream&, Matrix<Form, Type>&) : reading first token"
    );

    if (firstToken.isLabel())
    {
        M.mRows_ = firstToken.labelToken();
        M.nCols_ = readLabel(is);

        label mn = M.mRows_*M.nCols_;

        // Read list contents depending on data format
        if (is.format() == IOstream::ASCII || !contiguous<Type>())
        {
            // Read beginning of contents
            char listDelimiter = is.readBeginList("Matrix");

            if (mn)
            {
                M.allocate();
                Type* v = M.v_;

                if (listDelimiter == token::BEGIN_LIST)
                {
                    label k = 0;

                    // loop over rows
                    for (label i=0; i<M.m(); i++)
                    {
                        listDelimiter = is.readBeginList("MatrixRow");

                        for (label j=0; j<M.n(); j++)
                        {
                            is >> v[k++];

                            is.fatalCheck
                            (
                                "operator>>(Istream&, Matrix<Form, Type>&) : "
                                "reading entry"
                            );
                        }

                        is.readEndList("MatrixRow");
                    }
                }
                else
                {
                    Type element;
                    is >> element;

                    is.fatalCheck
                    (
                        "operator>>(Istream&, Matrix<Form, Type>&) : "
                        "reading the single entry"
                    );

                    for (label i=0; i<mn; i++)
                    {
                        v[i] = element;
                    }
                }
            }

            // Read end of contents
            is.readEndList("Matrix");
        }
        else
        {
            if (mn)
            {
                M.allocate();
                Type* v = M.v_;

                is.read(reinterpret_cast<char*>(v), mn*sizeof(Type));

                is.fatalCheck
                (
                    "operator>>(Istream&, Matrix<Form, Type>&) : "
                    "reading the binary block"
                );
            }
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "incorrect first token, expected <int>, found "
            << firstToken.info()
            << exit(FatalIOError);
    }

    return is;
}


template<class Form, class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const Matrix<Form, Type>& M)
{
    label mn = M.mRows_*M.nCols_;

    os  << M.m() << token::SPACE << M.n();

    // Write list contents depending on data format
    if (os.format() == IOstream::ASCII || !contiguous<Type>())
    {
        if (mn)
        {
            const Type* v = M.v_;

            // can the contents be considered 'uniform' (ie, identical)
            bool uniform = (mn > 1 && contiguous<Type>());
            if (uniform)
            {
                for (label i=0; i<mn; ++i)
                {
                    if (v[i] != v[0])
                    {
                        uniform = false;
                        break;
                    }
                }
            }

            if (uniform)
            {
                // Write start delimiter
                os  << token::BEGIN_BLOCK;

                // Write contents
                os << v[0];

                // Write end delimiter
                os << token::END_BLOCK;
            }
            else if (mn < 10 && contiguous<Type>())
            {
                // Write start contents delimiter
                os  << token::BEGIN_LIST;

                label k = 0;

                // loop over rows
                for (label i=0; i<M.m(); i++)
                {
                    os  << token::BEGIN_LIST;

                    // Write row
                    for (label j=0; j< M.n(); j++)
                    {
                        if (j) os << token::SPACE;
                        os << v[k++];
                    }

                    os << token::END_LIST;
                }

                // Write end of contents delimiter
                os << token::END_LIST;
            }
            else
            {
                // Write start contents delimiter
                os  << nl << token::BEGIN_LIST;

                label k = 0;

                // loop over rows
                for (label i=0; i<M.m(); i++)
                {
                    os  << nl << token::BEGIN_LIST;

                    // Write row
                    for (label j=0; j< M.n(); j++)
                    {
                        os << nl << v[k++];
                    }

                    os << nl << token::END_LIST;
                }

                // Write end of contents delimiter
                os << nl << token::END_LIST << nl;
            }
        }
        else
        {
            os  << token::BEGIN_LIST << token::END_LIST << nl;
        }
    }
    else
    {
        if (mn)
        {
            os.write(reinterpret_cast<const char*>(M.v_), mn*sizeof(Type));
        }
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
