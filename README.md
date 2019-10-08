# pileup2sam
A partial implementation of a [pileup](https://en.wikipedia.org/wiki/Pileup_format) to 
[SAM](https://www.ncbi.nlm.nih.gov/pubmed/19505943) format conversion tool.

- This is an inefficient and naive implementation, but it seems to do the trick.
- Does not support indels! 

## Example usage

```bash
pileup2sam -r test/cases/simulated.fa test/cases/simulated.head100.pileup out.sam
```
