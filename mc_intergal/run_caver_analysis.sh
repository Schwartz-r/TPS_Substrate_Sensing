#!/bin/bash

for ((i=1; i<=24; i++)); do
  for rep in 0 5; do
    analyze_p.py ${i} ${rep} &
  done
done

