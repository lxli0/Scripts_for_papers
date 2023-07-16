set GMT_SESSION_NAME=97401
gmt begin ETP_BS png

gmt grdimage @earth_relief_01m -JM15c -R70/130/10/50 -Baf -BWSen -Cbroc
gmt coast -R70/130/10/50 -JM15c -Baf -W1/0.5p,black
gmt colorbar -C -DJMR+o0.6c/0 -Bx1000+l"Topography (m)" -By -S


cat > TP.txt << EOF
104 35
104 32
101.8 28
98 28
98 35
104 35
EOF
gmt plot TP.txt -W1p,BLACK

cat > BB.txt << EOF
113 35
116 35
116.3 36.3
118.7 36.5
123 42
119 39.8
117 40
114.7 38.6
114 35.7
113 35
EOF
gmt plot BB.txt -W1p,BLACK
    
gmt end show
