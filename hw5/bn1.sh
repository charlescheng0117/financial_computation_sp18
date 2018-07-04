FILE=graph.in

echo "compiling source code"
make
echo "run binomial model with different M"
./hw5_bn1 < $FILE
echo "plotting option value v.s. M"
python3 hw5_plot.py
echo "opening plots"
open *.png
echo "Done"
