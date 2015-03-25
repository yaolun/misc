#!/bin/bash

echo "What is the lowest model number:"
read low
echo "What is the highest model number: "
read hi

for i in $(seq $low $hi)
do
  scp "yaolun@bettyjo.as.utexas.edu:~/hyperion/bhr71/cycle2/model"$i"/*.pdf" .
  scp "yaolun@bettyjo.as.utexas.edu:~/hyperion/bhr71/cycle2/model"$i"/*.png" .
done