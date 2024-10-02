library(hexSticker)
library(magick)
library(sysfonts)
library(tidyverse)
library(here)

setwd(here("hex_sticker"))

ally <- image_read("mean_gator.png")

# font_add_google("Fira Code", "firacode")
font_add_google("Knewave", "knewave")

sticker(
  subplot=ally,
  package="crawl",
  p_y = 1.4,
  p_size = 32,
  p_color = "khaki2",
  p_family = "knewave",
  s_x=1,
  s_y=1,
  s_width=1.4,
  s_height=1.4,
  h_fill="skyblue2",
  h_color="darkslategray4",
  dpi=300
)
