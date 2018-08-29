# BrownCOBRE2018

Teaching materials for short course

# Converting Rmd scripts to ipynb

This script is convenient. It uses "notedown" which I installed using `brew install notedown`:

```
#!/bin/bash
oldfile=$1
newfile="$(basename "$oldfile" .Rmd).ipynb"
echo Creating "$newfile" from "$oldfile"
notedown $oldfile | grep -v '%%R' > $newfile
```
