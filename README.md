#KalignP

##Description
multiple sequence alignment supporting external supplied position specific gap penalties

###Author
Nanjiang Shu

Department of Biochemistry and Biophysics 
Stockholm University

Email: nanjiang.shu@scilifelab.se


This software is derived from 
Kalign version 2.03, Copyright (C) 2006 Timo Lassmann
http://msa.cgb.ki.se/
timolassmann@gmail.com

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

A copy of this license is in the LICENSE file.

##Installation:

    $ ./configure
    $ make
    $ sudo make install


It is recommenced to use the Intel compiler 'icc' to compile kalignP which increase the speed by about 15%.
To enable the icc compiler, change the line in the Makefile

    `CC					= gcc  `
to 

    `CC					= icc  `

and rerun the compilation by the command

    $ make -B


##Usage

        `kalignP [Options]  infile.fasta outfile.fasta`

        or:

        `kalignP [Options] -i infile.fasta -o outfile.fasta`

        or:

        `kalignP [Options] < infile.fasta > outfile.fasta`

        For more options, type

        `kalignP -h`

##How KalignP works

KalignP changes the behavior of the alignment intuitively with the externally supplied position specific gap penalties. Generally speaking, if the users want to force a gap after the first residue of a sequence, set the gap open penalty at the second residue position of that sequence to a small value. For example, for the following four sequences,

```
    >seq1
    ASNLSKLFLSDSDA
    >seq2
    ASNLDA
    >seq3
    ASNLKFFFDDDAA
    >seq4
    LLNFFSDAAAAA
```

The multiple sequence alignment (MSA) with all parameters set to default is as follows (in ClustalW format)

```
    seq1 ASNLSKLFLSDSDA-
    seq2 ASNLDA-----
    seq3 ASNLKF-FFDDDAA-
    seq4 LLN--FFSDAAAAA 
```

The tree of the MSA is (((seq1, seq2), seq3), seq4). If the users want to open a gap after the second residue in seq1, the users can set the gap open penalty at the third position to be a minus value (e.g. -5). An example setting of ESPSGP is as follows. To ensure that only one gap is opened after the second residue S in seq1 so that we can see clearly the effect of gap open, we have set the gap open penalties for seq2 and the gap open penalty at the fourth position of seq1 to be a large positive value.

```
    >seq1
    ASNLSKLFLSDSDA
    {gpo: 1 1 -5 10 1 1 1 1 1 1 1 1 1 1 }
    >seq2
    ASNLDA
    {gpo: 10 10 10 10 10 10 }
    >seq3
    ASNLKFFFDDDAA
    >seq4
    LLNFFSDAAAAA
```

The alignment will become

```
    seq1 AS-NLSKLFLSDSDA-
    seq2 -ASNLDA-----
    seq3 -ASNLKF-FFDDDAA-
    seq4 --LLN-FFSDAAAAA
```

If the users want again to open a gap after the 11th residue D in seq3, set the gap open penalty of the 11th position in seq3 to a minus value (e.g. -5). An example ESPSGP setting is as follows. Again, to ensure that gap will only be opened after the 11th residue in seq3, we have set the gap open penalty of seq4 and the neighboring residue positions in seq3 to large positive values.

```
    >seq1
    ASNLSKLFLSDSDA
    {gpo: 1 1 -5 10 1 1 1 1 1 1 1 1 1 1 }
    >seq2
    ASNLDA
    {gpo: 10 10 10 10 10 10 }
    >seq3
    ASNLKFFFDDDAA
    {gpo: 1 1 1 1 1 1 1 1 1 10 -5 20 10 }
    >seq4
    LLNFFSDAAAAA
    {gpo: 10 10 10 10 10 10 10 10 10 10 10 10}
```

The alignment will become

```
    seq1 -----AS-NLSKLFLSDSDA
    seq2 -----ASNLDA----
    seq3 ASNLKFFFDD--DAA----
    seq4 ------LLNFFSDAAAAA
```

The alignment above might not be ideal since there are many gaps at the terminals. To reduce the number of terminal gaps, one can increase the terminal gap extension penalty as shown in the example below.

```
    >seq1
    ASNLSKLFLSDSDA
    {gpo: 1 1 -5 10 1 1 1 1 1 1 1 1 1 1 }
    {tgpe: 10 10 10 10 10 10 10 10 10 10 10 10 10 10}
    >seq2
    ASNLDA
    {gpo: 10 10 10 10 10 10 }
    {tgpe:10 10 10 10 10 10}
    >seq3
    ASNLKFFFDDDAA
    {gpo: 1 1 1 1 1 1 1 1 1 10 -5 20 10}
    {tgpe: 10 10 10 10 10 10 10 10 10 10 10 10 10}
    >seq4
    LLNFFSDAAAAA
    {gpo: 10 10 10 10 10 10 10 10 10 10 10 10}
    {tgpe:10 10 10 10 10 10 10 10 10 10 10 10}
```

The alignment will become

```
    seq1 AS-NLSKLFLSDS-DA
    seq2 -ASNLDA-----
    seq3 AS-NLKFFFDD-D-AA
    seq4 --LLNFFSDAAAAA-
```

Sometimes, the gap will not be opened at the expected position if many customized gap penalties are set in multiple sequences. This is because KalignP forbid neighboring gap opens at two aligned sequences such as the following example

```
    ALDDS-D-S
    ALED-D-S-
```

*Updated 2011-03-04 by Nanjiang Shu, nanjiang@sbc.su.se
