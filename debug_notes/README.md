1. CrisprVRiants failed for some references: ep300a, etv2, h3f3d, hey2:
Cause: the reference fasta contains Windows newline character that leads psl to bam conversion fault.
Solution: when converting ffom psl to bam, remove trailing spaces:
```
ref_name=\$(awk 'NR==1 {print substr(\$0,2,length(\$0))}' $ref)
# remove all trailing spaces in the reference name because samtools is very sensitive about it, any trailing spaces will invalidate the reference check, thereby, the third column will be assigned a *, and will be treated as "unmapped" by samtools (some ref.fasta contains a special ending character)
ref_name=\$(echo \$ref_name | awk '{ gsub(/[[:space:]]+\$/,""); print }')
```

Side note:
1. use `cat -A` to display all special characters (only works on Linux cat, not MacOS)
