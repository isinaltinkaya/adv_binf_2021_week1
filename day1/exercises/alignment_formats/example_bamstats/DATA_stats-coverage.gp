
            set terminal png size 600,400 truecolor
            set output "DATA_stats-coverage.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set ylabel "Number of mapped bases"
            set xlabel "Coverage"
            set log y
            set style fill solid border -1
            set title "DATA_stats.txt"
            set xrange [:5]
            plot '-' with lines notitle
        1	141896
2	12155
3	857
4	443
5	205
6	93
7	55
8	15
9	22
10	10
14	20
15	6
16	2
17	21
18	24
end
