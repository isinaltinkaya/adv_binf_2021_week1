
        set terminal png size 600,400 truecolor
        set output "DATA_stats-indel-dist.png"
        set grid xtics ytics y2tics back lc rgb "#cccccc"
        set style line 1 linetype 1  linecolor rgb "red"
        set style line 2 linetype 2  linecolor rgb "black"
        set style line 3 linetype 3  linecolor rgb "green"
        set style increment user
        set ylabel "Indel count [log]"
        set xlabel "Indel length"
        set y2label "Insertions/Deletions ratio"
        set log y
        set y2tics nomirror
        set ytics nomirror
        set title "DATA_stats.txt"
        plot '-' w l ti 'Insertions', '-' w l ti 'Deletions', '-' axes x1y2 w l ti "Ins/Dels ratio"
    1	18
2	4
3	1
4	0
5	5
end
1	29
2	7
3	0
4	2
5	0
end
1	0.620690
2	0.571429
3	0.000000
4	0.000000
5	0.000000
end
