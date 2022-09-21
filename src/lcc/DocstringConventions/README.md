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


`value` -- indicates another LPP file 
(hyperlinkTitle)[hyperlink] -- indicates text with associated hyperlink

In line reference IDs will be **in future** automatically parsed into the journal citation format of your choice.  

Ex: 
```sh
""" Sigma Factors  

There are 7 sigma factors in E. Coli which facilitate `RNAP` binding during transcription initiation /* RefID1 */. Here, each sigma factor is defined on the basis of its Kd and promoter binding region.

Discussion
----------
I have this and that to say about this process that is more opinion or experimental observation. This model does better than X model by /* RefID2 */ as shown by x, y, z.

See Also
--------
`Transcription` : LPP implementation of transcription 
`RNAP` : LPP implementation of RNA polymerase

Background
----------
(Transcription)[https://en.wikipedia.org/wiki/Transcription_(biology)] : A Wikipedia article describing transcription.
(RNAP)[https://en.wikipedia.org/wiki/RNA_polymerase] : A Wikipedia article describing RNA polymerase

Database
--------
Uniprot : (UniprotID)[link to protein in uniprot]
Ecocyc : (EcocycID)[link to Ecocyc]

About
-----
Author : 
Journal : Not published
Version : v0.1
VersionRelease : 

References
----------
/* RefID1 */ : Some reference to my claim in docstring.
/* RefID1 */ : Burgess R.R.(2001) Sigma Factors _Encyclopedia of Genetics_ doi: 10.1006/rwgn.2001.1192
    table 1, sigma factor promoter regions
/* RefID1 */ : Burgess R.R.(2001) Sigma Factors _Encyclopedia of Genetics_ doi: 10.1006/rwgn.2001.1192
    Kd, RpoD
    Kd, RpoE
    ...
/* RefID2 */ : SomeAuthor (2072) Some Title about a Model _Some Journal_ doi: someDoi

"""

# Some code
sigmafactor RpoD(ds(TTGACANNNN NNNNNNNNNN NNNTATAATN), Kd=1e-15); /* 1 */

```


<p align="right">(<a href="#top">back to top</a>)</p>
