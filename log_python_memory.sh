#!/bin/bash

LOGFILE="python_memory_$(date +%Y%m%d_%H%M%S).log"
echo "Logging to: $(realpath "$LOGFILE")"


VSZ_MAX=0
RSS_MAX=0


while true; do
    ps -u ashmat -C python3 -o comm,vsz,rss --no-headers | awk -v vsz_max="$VSZ_MAX" -v rss_max="$RSS_MAX" -v logfile="$LOGFILE" '
        $1=="python3" {
            vsz = $2; rss = $3
            if (vsz > vsz_max) vsz_max = vsz
            if (rss > rss_max) rss_max = rss
            cmd = "date +\"%Y-%m-%d %H:%M:%S\""; cmd | getline dt; close(cmd)
            printf "%s %s VSZ=%.2fGB (max=%.2fGB) RSS=%.2fGB (max=%.2fGB)\n", dt, $1, vsz/1048576, vsz_max/1048576, rss/1048576, rss_max/1048576 | "tee -a " logfile
        }
        END { print vsz_max > "/tmp/_vsz_max"; print rss_max > "/tmp/_rss_max" }
    '
    VSZ_MAX=$(cat /tmp/_vsz_max 2>/dev/null || echo $VSZ_MAX)
    RSS_MAX=$(cat /tmp/_rss_max 2>/dev/null || echo $RSS_MAX)

done
