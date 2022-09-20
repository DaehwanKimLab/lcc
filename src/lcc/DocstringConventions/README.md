<!-- L++ Docstring Syntax -->
<details>
  <summary>L++ Docstring Syntax</summary>
  <ol>
    <li><a href="#Introduction">Introduction</a></li>
    <li><a href="#File-Docstrings">File Docstrings</a></li>
    <ol>
        <li><a href="#Basics">Basics</a></li>
        <li><a href="#References">References</a></li>
    </ol>
    </li> 
  </ol>
</details>

<!-- Introduction -->
## Introduction

L++ as a language comes with its own somewhat unique properties of docstring format. Perhaps most notably is the frequent reliance on referenced sources. File authorship? (almost treating files like journal articles?? -- maybe this is a bit down the road but I think eventually lpp files may be like short journal articles.)

The subsequent text will discribe the formatting of docstrings in L++.

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- File-Docstrings -->
## File Docstrings

All files should have a preferably one line description of what protein/pathway/process/etc. is being described. In short, what is the file supposed to convey?

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- References -->
## References

Types of references:
- Background: In general, L++ more than other programming languages tends to be built upon biological language and syntax. References of this type should give basic background information of the file's topic. These however should `NOT` be placed in the references section.  Instead, these should be placed in the `Background` portion of the docstring.

- Reference material:
Per datacite guidelines (https://datacite.org/cite-your-data.html), minimal information provided as follows
> Author (Year) Title Journal/Book/Publisher Identifier
Burgess R.R. (2001) Sigma Factors _Encyclopedia of Genetics_ doi: 10.1006/rwgn.2001.1192

If you wish to inline cite material (such as within docstring)

If a specific value is cited (such as a Kd), this should be specified in the reference section, if multiple values are cited from the same source...

Ex: 
```sh
""" These Proteins are involved in Transcription initiation /* 1 */

A Short description of the protein. This protein is involved in glycolysis 

Comments
--------
I have this and that to say about this process that is more opinion or experimental observation.

See Also
--------
`Glycolysis` 

Background
----------
(Glycolysis)[https://en.wikipedia.org/wiki/Glycolysis] : A Wikipedia article describing glycolysis.
(Glycolysis)[https://www.sciencedirect.com/topics/neuroscience/glycolysis] : Science Direct glycolysis topic  

Database
--------
Uniprot : (UniprotID)[link to protein in uniprot]
Ecocyc : (EcocycID)[link to Ecocyc]

References
----------
/* 1 */ 
    Some reference to my claim in docstring.
/* 2 */ : table 1, sigma factor promoter regions
    Burgess R.R.(2001) Sigma Factors _Encyclopedia of Genetics_ doi: 10.1006/rwgn.2001.1192
/* 3 */ : Kd, RpoD
    Burgess R.R.(2001) Sigma Factors _Encyclopedia of Genetics_ doi: 10.1006/rwgn.2001.1192
/* 3 */ : Kd, RpoE
    Burgess R.R.(2001) Sigma Factors _Encyclopedia of Genetics_ doi: 10.1006/rwgn.2001.1192

"""
```


<p align="right">(<a href="#top">back to top</a>)</p>
