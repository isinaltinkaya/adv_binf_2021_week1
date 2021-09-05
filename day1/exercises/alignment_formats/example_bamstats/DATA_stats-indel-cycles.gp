
        set terminal png size 600,400 truecolor
        set output "DATA_stats-indel-cycles.png"
        set grid xtics ytics y2tics back lc rgb "#cccccc"
        set style line 1 linetype 1  linecolor rgb "red"
        set style line 2 linetype 2  linecolor rgb "black"
        set style line 3 linetype 3  linecolor rgb "green"
        set style line 4 linetype 4  linecolor rgb "blue"
        set style increment user
        set ylabel "Indel count"
        set xlabel "Read Cycle"
        set title "DATA_stats.txt"
    plot '-' w l ti 'Insertions (fwd)', '' w l ti 'Insertions (rev)', '' w l ti 'Deletions (fwd)', '' w l ti 'Deletions (rev)'
3	0
4	0
5	0
6	0
7	0
8	0
9	0
12	0
13	0
14	0
15	0
16	0
17	0
18	0
19	0
20	0
22	0
23	0
24	0
25	1
27	0
28	0
29	0
31	0
33	0
34	0
37	0
38	0
41	0
42	0
45	0
49	0
58	0
63	0
77	0
95	0
end
3	0
4	1
5	0
6	0
7	0
8	1
9	0
12	1
13	1
14	1
15	2
16	0
17	0
18	0
19	1
20	1
22	2
23	1
24	1
25	0
27	0
28	1
29	2
31	0
33	0
34	0
37	3
38	2
41	0
42	1
45	1
49	2
58	0
63	0
77	1
95	1
end
3	0
4	0
5	0
6	0
7	0
8	0
9	0
12	0
13	0
14	0
15	0
16	0
17	0
18	0
19	0
20	0
22	0
23	0
24	0
25	0
27	0
28	0
29	0
31	0
33	0
34	0
37	0
38	0
41	0
42	0
45	0
49	0
58	0
63	0
77	0
95	0
end
3	1
4	0
5	1
6	1
7	1
8	1
9	1
12	1
13	1
14	0
15	3
16	1
17	1
18	2
19	2
20	2
22	0
23	1
24	1
25	1
27	1
28	2
29	3
31	1
33	1
34	1
37	1
38	0
41	2
42	0
45	1
49	1
58	1
63	1
77	0
95	0
end
