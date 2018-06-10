#!/usr/bin
INPUT=100-2.in
echo "Basic requirement"
./binom $INPUT
./monte $INPUT
echo "Bonus 1"
time python3 hw4_bonus1.py $INPUT
echo "Bonus 2"
time python3 hw4_bonus2.py $INPUT
