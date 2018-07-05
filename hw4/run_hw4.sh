#!/usr/bin
echo "S_max_t = 50"
echo "Basic requirement"
./binom q50.in
echo ""
./monte q50.in
echo ""
echo "Bonus 1"
time python3 hw4_bonus1.py q50.in
echo ""
echo "Bonus 2"
time python3 hw4_bonus2.py q50.in
echo "50 Done"
echo ""
echo "S_max_t = 60"
echo "Basic requirement"
./binom q60.in
echo ""
./monte q60.in
echo ""
echo "Bonus 1"
time python3 hw4_bonus1.py q60.in
echo "60 Done"
echo ""
echo "S_max_t = 70"
echo "Basic requirement"
./binom q70.in
echo ""
./monte q70.in
echo ""
echo "Bonus 1"
time python3 hw4_bonus1.py q70.in
echo "70 Done"
