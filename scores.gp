set multiplot layout 2, 2 title "Évolution Scores"
set tmargin 2

set title "ΣScores"
unset key
#plot "../out/evo_score" using 1:2 w l, "../out/evo_score" using 1:4 w l, "../out/evo_score" using 1:3 w l
plot "../out/evo_score" using 1:4 w l, "../out/evo_score" using 1:5 w l

set title "Power Degree Law"
unset key
#plot "../out/evo_pdl" using 1:2 w l, "../out/evo_pdl" using 1:4 w l, "../out/evo_pdl" using 1:3 w l
plot "../out/evo_pdl" using 1:4 w l, "../out/evo_pdl" using 1:5 w l

set title "Small World"
unset key
#plot "../out/evo_sw" using 1:2 w l, "../out/evo_sw" using 1:4 w l, "../out/evo_sw" using 1:3 w l
plot "../out/evo_sw" using 1:4 w l, "../out/evo_sw" using 1:5 w l

set title "Clique Formation"
unset key
#plot "../out/evo_cf" using 1:2 w l, "../out/evo_cf" using 1:4 w l, "../out/evo_cf" using 1:3 w l
plot "../out/evo_cf" using 1:4 w l, "../out/evo_cf" using 1:5 w l

unset multiplot

pause -1
