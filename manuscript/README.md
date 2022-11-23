# Manuscript

Rebuild with `Rake`. This will run:

```
pdflatex -draftmode ncov-escape
pdflatex -draftmode ncov-escape
bibtex ncov-escape
pdflatex ncov-escape
```

to generate the PDF `ncov-escape.pdf`. However, on subsequent builds it will skip steps if they are not required.
The PDF `ncov-escape.pdf` is intentionally not versioned.
