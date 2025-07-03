for dir in dist*_sim*; do
  dx=$(echo "$dir" | sed -E 's/dist([0-9]+)_sim[0-9]+/\1/')
  sim=$(echo "$dir" | sed -E 's/dist[0-9]+_sim([0-9]+)/\1/')
  expected="sim_seqs_maxDist_${dx}_${dx}_Fasta.fa"
  if [ ! -f "$dir/$expected" ]; then
    echo "[MISSING] $expected not found in $dir (dist=$dx, sim=$sim)"
  fi
done
