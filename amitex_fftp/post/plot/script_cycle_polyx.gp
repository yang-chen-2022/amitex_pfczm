set t svg size 1100 1100 fsize 35;
set xtics rotate by 25;
set xtics offset -3.5,graph -0.11;
# Scientific notation
set format x "%.1e";
set grid;
set key off;
set o "cycle_polyxCC_Exx_Sxx.svg";
set title "Cycling load on polycrystal";
set xlabel "strain_{11}";
set ylabel "stress_{11}";
plot "../../resultats/polyxCC_cycle_traction/polyxCC_cycle_traction.std" u 8:2 w l lw 4 notitle;
