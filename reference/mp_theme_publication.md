# Default publication theme for micropubr figures

A clean, minimal `ggplot2` theme with no panel grid clutter, framed axes
and balanced typography. It is applied automatically by every `mp_*()`
plotting function and can be customised or replaced via the shared
`theme` argument.

## Usage

``` r
mp_theme_publication(
  font.size = getOption("micropubr.font.size", 11),
  base.family = "",
  grid = c("none", "x", "y", "both"),
  legend.position = "right"
)
```

## Arguments

- font.size:

  Base font size in points. Defaults to the `micropubr.font.size` option
  (11).

- base.family:

  Font family. Defaults to "" (device default).

- grid:

  Which major grid lines to draw: `"none"` (default), `"x"`, `"y"` or
  `"both"`.

- legend.position:

  Legend position passed to
  [`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html).

## Value

A `ggplot2` theme object.

## Examples

``` r
library(ggplot2)
ggplot(mtcars, aes(wt, mpg)) +
  geom_point() +
  mp_theme_publication()
```
