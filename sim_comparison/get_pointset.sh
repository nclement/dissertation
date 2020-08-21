pts=$1

# Cols 2, 7, 8
min=($(awk 'NR==1{min[2] = $2 + 0; min[7] = $7; min[8]=$8; next} {if ($2 < min[2]) min[2] = $2; if ($7 < min[7]) min[7] = $7; if ($8 < min[8]) min[8] = $8;} END {printf "%f %f %f", min[2], min[7], min[8]}' $pts))
max=($(awk 'NR==1{max[2] = $2 + 0; max[7] = $7; max[8]=$8; next} {if ($2 > max[2]) max[2] = $2; if ($7 > max[7]) max[7] = $7; if ($8 > max[8]) max[8] = $8;} END {printf "%f %f %f", max[2], max[7], max[8]}' $pts))

awk -v min2="${min[0]}" -v min7="${min[1]}" -v min8="${min[2]}" \
  -v max2="${max[0]}" -v max7="${max[1]}" -v max8="${max[2]}" \
  '{printf "%f %f %f\n", ($2-min2)/(max2-min2), ($7-min7)/(max7-min7), ($8-min8)/(max8-min8)}' $pts | head -n 50
