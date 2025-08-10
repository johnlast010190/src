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
    (c) 2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

inline bool Foam::wordRe::meta(const char c)
{
    return regExp::meta(c);
}


inline bool Foam::wordRe::isPattern(const std::string& str)
{
    return string::meta<regExp>(str);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::wordRe::wordRe()
:
    word(),
    re_()
{}


inline Foam::wordRe::wordRe(const wordRe& str)
:
    word(str, false),
    re_()
{
    if (str.isPattern())
    {
        compile();
    }
}


inline Foam::wordRe::wordRe(const keyType& str)
:
    word(str, false),
    re_()
{
    if (str.isPattern())
    {
        compile();
    }
}


inline Foam::wordRe::wordRe(const char* str)
:
    word(str, false),
    re_()
{}


inline Foam::wordRe::wordRe(const std::string& str)
:
    word(str, false),
    re_()
{}


inline Foam::wordRe::wordRe(const string& str)
:
    word(str, false),
    re_()
{}


inline Foam::wordRe::wordRe(const word& str)
:
    word(str, false),
    re_()
{}


inline Foam::wordRe::wordRe(const keyType& str, const compOption opt)
:
    word(str, false),
    re_()
{
    if (str.isPattern())
    {
        compile(opt);
    }
}


inline Foam::wordRe::wordRe(const char* str, const compOption opt)
:
    wordRe(str)
{
    compile(opt);
}


inline Foam::wordRe::wordRe(const std::string& str, const compOption opt)
:
    wordRe(str)
{
    compile(opt);
}


inline Foam::wordRe::wordRe(const string& str, const compOption opt)
:
    wordRe(str)
{
    compile(opt);
}


inline Foam::wordRe::wordRe(const word& str, const compOption opt)
:
    wordRe(str)
{
    compile(opt);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::wordRe::isPattern() const
{
    return re_.exists();
}


inline bool Foam::wordRe::compile(const compOption opt) const
{
    if (opt)
    {
        bool comp = false;

        if (opt & wordRe::REGEX)
        {
            comp = true;
        }
        else if (opt & wordRe::DETECT)
        {
            comp = string::meta<regExp>(*this) || !string::valid<word>(*this);
        }
        else if (opt & wordRe::NOCASE)
        {
            comp = true;
        }

        if (comp)
        {
#if defined(WIN64) || defined(WIN32)
            re_.set(*this, (opt & wordRe::NOCASE));
            #else
            return re_.set(*this, (opt & wordRe::NOCASE));
            #endif
        }
#if defined(WIN64) || defined(WIN32)
                return re_.exists();
                #endif
    }

    // Fall-through behaviour - not a regex
    re_.clear();
    return false;
}


inline bool Foam::wordRe::compile() const
{
#if defined(WIN64) || defined(WIN32)
    re_ = *this;
    return re_.exists();
#else
    return re_.set(*this);
#endif
}


inline void Foam::wordRe::uncompile(const bool doStripInvalid) const
{
    if (re_.clear() && doStripInvalid)
    {
        // Skip stripping unless debug is active to avoid costly operations
        if (word::debug)
        {
            string::stripInvalid<word>
            (
                const_cast<word&>(static_cast<const word&>(*this))
            );
        }
    }
}


inline void Foam::wordRe::clear()
{
    word::clear();
    re_.clear();
}


inline bool Foam::wordRe::match(const std::string& text, bool literal) const
{
    if (literal || !re_.exists())
    {
        return !compare(text);  // Compare as literal string
    }
    else
    {
        return re_.match(text); // Match as regex
    }
}


inline Foam::string Foam::wordRe::quotemeta() const
{
    return string::quotemeta<regExp>(*this);
}


inline void Foam::wordRe::set(const std::string& str, const compOption opt)
{
    string::operator=(str);
    compile(opt);
}


inline void Foam::wordRe::set(const char* str, const compOption opt)
{
    string::operator=(str);
    compile(opt);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::wordRe::operator()(const std::string& text) const
{
    return match(text);
}


inline void Foam::wordRe::operator=(const wordRe& str)
{
    string::operator=(str);
    if (str.isPattern())
    {
        compile();
    }
    else
    {
        re_.clear();
    }
}


inline void Foam::wordRe::operator=(const word& str)
{
    string::operator=(str);
    re_.clear();
}


inline void Foam::wordRe::operator=(const keyType& str)
{
    string::operator=(str);
    if (str.isPattern())
    {
        compile();
    }
    else
    {
        re_.clear();
    }
}


inline void Foam::wordRe::operator=(const string& str)
{
    string::operator=(str);
    compile(wordRe::DETECT);  // Auto-detect regex
}


inline void Foam::wordRe::operator=(const std::string& str)
{
    string::operator=(str);
    compile(wordRe::DETECT);  // Auto-detect regex
}


inline void Foam::wordRe::operator=(const char* str)
{
    string::operator=(str);
    compile(wordRe::DETECT);  // Auto-detect regex
}


// ************************************************************************* //
