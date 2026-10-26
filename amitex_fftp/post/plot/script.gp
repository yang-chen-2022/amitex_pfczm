# Output in an SVG image
set t svg size 1100 1100 fsize 35;
# Reduction of the margins
# Rotation of the value printed below the x axis (more readable)
set xtics rotate by 25;
set xtics offset -3.5,graph -0.11;
# Scientific notation
set format x "%.1e";
set grid;
set key off;
# plot 1: Imposed strain on polycrystal
set o "polyxCC_Exx_Sxx.svg";
set title "Imposed strain on polycrystal";
set xlabel "strain_{11}";
set ylabel "stress_{11}";
plot "../../resultats/polyxCC_def_imp/polyxCC_def_imp.std" u 8:2 w l notitle lw 3;
# plot 2: Imposed strain and relaxation of concrete
set o "concrete_Exx_Sxx_mat.svg";
set title "Imposed strain and relaxation of concrete";
set xlabel "strain_{11}";
set ylabel "stress_{11}";
plot "../../resultats/beton_relax_65/beton_relax_65.mstd" every 2::0 u 8:2 w l notitle lw 3, "../../resultats/beton_relax_65/beton_relax_65.mstd" every 2::1 u 8:2 w l notitle lw 3;
# plot 3: Creep of concrete
set o "concrete_creep_std_zone.svg";
set title "Creep of concrete";
set xlabel "t";
set ylabel "StandardDeviation(stress_{11})";
# Requires gnuplot 4.6 or later
plot for [j=1:1110] "../../resultats/beton_fluage_64/beton_fluage_64_2.zstd" every 1110::j-1 u 1:14 w l notitle lw 3;
