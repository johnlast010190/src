/*---------------------------------------------------------------------------*\
| Modified 2010-2012 Copyright (C) Esi Ltd                                  |
\*---------------------------------------------------------------------------*/
void indent(label nin, OFstream& dout)
{
    for (label i = 0; i < nin; i++)
    {
        dout << token::SPACE << token::SPACE
             << token::SPACE << token::SPACE << flush;
    }
}

void writeEntries
(
    OFstream& dout,
    const dictionary& dict,
    word kw,
    label nIndent
);


void writeDict
(
    OFstream& dout,
    const dictionary& dict,
    label nIndent,
    label newlines = 0
)
{
    wordList dicToc(dict.toc());
    forAll(dicToc, dtI)
    {
        word keyword = dicToc[dtI];

        writeEntries(dout, dict, keyword, nIndent);

        for (label i = 0; i < newlines; i++)
        {
            dout << nl;
        }
    }
}

void writeEntries
(
    OFstream& dout,
    const dictionary& dict,
    word kw,
    label nIndent
)
{
    if (dict.isDict(kw))
    {
        indent(nIndent, dout);
        dout << kw << nl;
        indent(nIndent, dout);
        dout << token::BEGIN_BLOCK << nl;
        nIndent++;

        const dictionary& sdict = dict.subDict(kw);
        wordList dicToc(sdict.toc());
        forAll(dicToc, dtI)
        {
            word keyword = dicToc[dtI];
            writeEntries(dout, sdict, keyword, nIndent);
        }

        nIndent--;
        indent(nIndent, dout);
        dout << token::END_BLOCK << nl;
    }
    else
    {
        //make an ITstream from entry.
        //Check whether it has at least 3 tokens
        //Then check whether it is a list of dictionaries

        ITstream& checkTS = dict.lookup(kw);

        List<token> tsTokens(3);

        forAll(tsTokens, tI)
        {
            if (checkTS.eof())
            {
                break;
            }
            tsTokens[tI] = token(checkTS);
        }

        //if a list, assume all entries are uniform
        //check if entries are dictionaries
        if
        (
            tsTokens[0] == token::BEGIN_LIST
            &&
            (
                tsTokens[1] == token::BEGIN_BLOCK
                || (tsTokens[1].isWord() && tsTokens[2] == token::BEGIN_BLOCK)
            )
        )
        {
            //this is a list of named/unnamed dictionaries
            indent(nIndent, dout);
            dout << kw << nl;
            indent(nIndent, dout);
            dout << token::BEGIN_LIST << flush;
            nIndent++;

            //now check if it is a named or unnamed list
            if (tsTokens[1].isWord()) // named list
            {
                //stream to output
                ITstream& its = dict.lookup(kw);

                //strip leading bracket
                token t(its);

                t = token(its);

                while  (t != token::END_LIST)
                {
                    dictionaryEntry subDict(t.wordToken(), dict, its);

                    dout << nl;
                    indent(nIndent, dout);
                    dout << t << nl;
                    indent(nIndent, dout) ;
                    dout << token::BEGIN_BLOCK << nl;
                    nIndent++;

                    writeDict(dout, subDict, nIndent);

                    nIndent--;
                    indent(nIndent, dout);
                    dout << token::END_BLOCK << endl;

                    if (its.eof())
                    {
                        break;
                    }
                    t = token(its);
                }
            }
            else //unnamed list
            {
                List<dictionary> dictList(dict.lookup(kw));

                forAll(dictList, dlI)
                {
                    const dictionary& cDict = dictList[dlI];
                    dout << nl;
                    indent(nIndent, dout) ;
                    dout << token::BEGIN_BLOCK << nl;
                    nIndent++;

                    writeDict(dout, cDict, nIndent);

                    nIndent--;
                    indent(nIndent, dout) ;
                    dout << token::END_BLOCK << endl;;
                }
            }

            nIndent--;
            indent(nIndent, dout);
            dout << token::END_LIST << token::END_STATEMENT
                 <<  endl;
        }
        else
        {
            const entry* pePtr = &dict.lookupEntry(kw,false,true);

            indent(nIndent, dout);
            pePtr->write(dout);
            dout << flush;
        }
    }
}

void writeDictFile(const IOdictionary& dict)
{
    OFstream dout(dict.IOobject::objectPath());

    dout  << "FoamFile\n{\n";
    dout  << "    version " << dout.version() << ";\n";
    dout  << "    format " << dout.format() << ";\n\n";

    dout  << "    root " << dict.rootPath() << ";\n";
    dout  << "    case " << dict.caseName() << ";\n";
    dout  << "    instance " << dict.instance() << ";\n";
    dout  << "    local " << dict.local() << ";\n\n";

    dout  << "    class " << dict.type() << ";\n";

    if (dict.note().size())
    {
        dout << "\n    note " << dict.note() << ";\n\n";
    }

    dout  << "    object " << dict.name() << ";\n";
    dout  << "}\n\n";
    dict.writeDivider(dout);
    dout  << endl;

    label nIndent = 0;

    writeDict(dout, dict, nIndent, 1);
}
