---
name: Bug report
about: Report something that does not work as documented
title: ""
labels: bug
assignees: ""
---

**Describe the bug**
A clear, short description of what went wrong.

**Minimal reproducible example**
Please use a small or public dataset (e.g. `data(GlobalPatterns, package = "phyloseq")`)
and, ideally, `reprex::reprex()`.

```r
library(micropubr)
data(GlobalPatterns, package = "phyloseq")

# code that triggers the problem
```

**Expected behaviour**
What you expected to happen (for a figure, attach the output or describe the
difference).

**Actual behaviour**
The full error or warning message. For errors, the output of
`rlang::last_trace()` and the condition class (`class(rlang::last_error())`)
are very helpful.

**Session information**
<details>
<summary>sessionInfo()</summary>

```r
# paste the output of sessionInfo() here
```

</details>

**Additional context**
Anything else: data source (DADA2, QIIME 2, MetaPhlAn, ...), whether
`interactive = TRUE` was used, the output file type for `mp_export()`, etc.
