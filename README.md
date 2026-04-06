# bios731_hw4_taylor

## Generic Project Directory

This is a generic project template. It might help to add a brief description.

When the project has been created, I recommend running

```usethis::use_github(private = TRUE, protocol = "https")```

to create a private remote repo corresponding to this project. If this project
is edited across multiple machines, it may be necessary to add symbolic links
to data directories using

```ln -s PATH/TO/DATA data```

## Structure

* `analysis/` is intended to contain all the source R Markdown files that
implement the analyses for the project.
* `data/` will either be a symbolic link to an external data directory, or
a subdirectory
* `drafts/` will contain relevant paper drafts
* `literature/` will contain all references
* `results/` will contain results exported by the analysis files
* `source/` will contain bare scripts (typically containing functions sourced
by the full analysis files)

project structured using the [projectr](https://github.com/julia-wrobel/projectr) package.