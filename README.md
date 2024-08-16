# finished-ish
This repository contains working scripts I've written. I say finished-ish because (as of 7/2024) I am new to programming and may have made mistakes that cause the code to not work for others despite working for me. Additionally, I plan to improve my code as I become a better programmer :D

# CellProfiler to xy-plot pipeline for 96-well cell viability assays (7/19/24)
My first project is to automate the unwieldy process we have in my lab for going from images of 96-well plates to polished graphs displaying cell survival data. Normally, we would take images, count nuclei using CellProfiler, take the .csv file from CellProfiler and transform it in Excel, normalize it in Excel or GraphPad, and finally plot it in Graphpad. Needless to say, this is worthy of automation. My script "nuclei_count_to_xyplot_local_v1.0.0" works on my computer, but would need some tweaking with file directories, etc. to work on your computer. I plan on turning this into a fully packaged .exe file sometime in the near future! Hopefully it is well-annotated enough for someone with a bit of Python experience to adapt it as-is in the meantime.

# GDA synergy analysis pipeline (8/9/24)
This pipeline "GDA_synergy_local_v2.1" automates another unwieldy process in our lab, going from images of 96-well plates to a synergy analysis plot. This pipeline is designed for three, 96-well plates with DPBS in the outer wells and a 6x10 grid of two different drugs in gradient in the inner 60 wells. This pipeline, like the last one, is designed to work locally for me (James), so it would need some tweaking to get it to work on a different system. Essentially, it uses CellProfiler headless to count nuclei in 180 4x3 images of the 60 wells in triplicate. It then finds the mean for each dose combination, normalizes it to the vehicle-only mean, determines synergy using Bliss independence, then plots cell viability (z-axis), drug concentrations (x- and y-axes), and synergy score(color) on a mapped 3D surface.

# Minor updates to nuclei_count_to_xyplot_local_v1.0.0, major updates to GDA_synergy_local_v2.1 (8/16/24)
## nuclei_count_to_xyplot_local_v1.1
  Minor aesthetic and functional changes to the GUI.
  Created singular variables for IC50 values as opposed to leaving in array.
## GDA_synergy_local_v2.2
  Major functional changes, I do not recommend using v2.1 any longer.
  Made the artificial zero relative to each drug gradient. Made the plots spatially representative of the data (no longer forces a square representation).
  Added .csv outputs for data generated along the way, including pivot tables with normalized nuclei count means and Bliss independence values.
